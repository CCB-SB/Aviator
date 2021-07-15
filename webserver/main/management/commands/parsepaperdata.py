from collections import defaultdict

from django.core.management.base import BaseCommand, CommandError
from django.db.models import Count

import pandas as pd

from main.models import Website
from main.models import Publication


class Command(BaseCommand):
    help = 'Creates the paper model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')
        parser.add_argument('url_filter', default=None)
        parser.add_argument('biotoolsid2pmid', default=None)

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
        except:
            raise CommandError('Please provide an input folder')

        filter_orig_urls = set()
        if options["url_filter"] is not None and \
            options["url_filter"] not in {"None", "null", ""}:
            filter_tbl = pd.read_csv(options['url_filter'], sep='\t')
            filter_orig_urls = set(filter_tbl["Original URL"])
            self.stdout.write("Filter applied")
        else:
            self.stdout.write("No Filter applied")

        pmid2biotoolsid = None
        if options["biotoolsid2pmid"] is not None and \
            options["biotoolsid2pmid"] not in {"None", "null", ""}:
            df = pd.read_csv(options['biotoolsid2pmid'], sep='\t', names = ["biotoolsID","pmid"])
            pmid2biotoolsid = dict(zip(list(df.pmid), list(df.biotoolsID)))
            self.stdout.write("Applied biotoolsid2pmid")
        else:
            self.stdout.write("No biotoolsid2pmid")

        pub_tbl = pd.read_csv(csv_file, sep='\t', index_col="PMID")
        pub_tbl.index = pub_tbl.index.astype(str)
        pub_dict = pub_tbl.to_dict()

        # remove publications that are not in our list anymore
        pubs_to_remove = Publication.objects.exclude(pubmed_id__in=pub_tbl.index)
        pubs_to_remove.delete()
        self.stdout.write(self.style.SUCCESS(
            f"Removed {pubs_to_remove.count()} publications that are not in our table anymore. This can happen when duplicate records are deleted from PubMed."))

        papers_to_update = list(Publication.objects.filter(pubmed_id__in=pub_tbl.index))

        def split_w_nan(s, delimiter):
            if pd.isnull(s):
                return []
            else:
                return s.split(delimiter)

        for p in papers_to_update:
            p.title = pub_dict["title"][p.pubmed_id]
            p.abstract = pub_dict["abstract"][p.pubmed_id]
            p.year = pub_dict["year"][p.pubmed_id]
            p.authors = str(pub_dict["authors"][p.pubmed_id]).split(", ")
            p.journal = pub_dict["journal"][p.pubmed_id]
            urls = pub_dict["URL"][p.pubmed_id].split('; ')
            if filter_orig_urls:
                urls = [u for u in urls if u in filter_orig_urls]
            p.url = urls
            p.user_kwds = split_w_nan(pub_dict["keywords_all"][p.pubmed_id], ';')
            p.mesh_terms = split_w_nan(pub_dict["mesh_terms_all"][p.pubmed_id], ';')
            p.contact_mail = split_w_nan(pub_dict["Email"][p.pubmed_id], ';')
            if pmid2biotoolsid:
                if p.pubmed_id in pmid2biotoolsid:
                    p.biotools_id = pmid2biotoolsid[p.pubmed_id]
                else:
                    p.biotools_id = ""

        Publication.objects.bulk_update(papers_to_update,
                                        ["title", "abstract", "year", "authors", "journal", "url",
                                         "user_kwds", "mesh_terms", "contact_mail"])

        papers_to_update_ids = set(p.pubmed_id for p in papers_to_update)

        def get_paper_models(tbl):
            for i, row in tbl.iterrows():
                urls = row["URL"].split('; ')
                if filter_orig_urls:
                    urls = [u for u in urls if u in filter_orig_urls]
                biotools_id = ""
                if pmid2biotoolsid:
                    if i in pmid2biotoolsid:
                        biotools_id = pmid2biotoolsid[i]
                if len(urls) > 0:
                    yield Publication(pubmed_id=i,
                                      biotools_id = biotools_id,
                                      authors=str(row["authors"]).split(", "),
                                      title=row["title"],
                                      abstract=row["abstract"],
                                      year=row["year"],
                                      journal=row["journal"],
                                      url=urls,
                                      user_kwds=split_w_nan(row["keywords_all"], ';'),
                                      mesh_terms=split_w_nan(row["mesh_terms_all"], ';'),
                                      contact_mail=split_w_nan(row["Email"], ';')
                                      )

        new_pubs = pub_tbl[~pub_tbl.index.isin(papers_to_update_ids)]
        new_pub_mdls = list(get_paper_models(new_pubs))
        Publication.objects.bulk_create(new_pub_mdls)

        pmid2pub = {p.pubmed_id: p for p in Publication.objects.filter(pubmed_id__in=pub_tbl.index)}
        url2pub = defaultdict(set)
        for p, url_str in pub_dict["URL"].items():
            for u in url_str.split('; '):
                if filter_orig_urls and u in filter_orig_urls:
                    url2pub[u].add(pmid2pub[p])
                elif not filter_orig_urls:
                    url2pub[u].add(pmid2pub[p])

        updated_websites = 0
        websites_to_update = Website.objects.filter(original_url__in=url2pub.keys())
        for website in websites_to_update:
            website.papers.add(*url2pub[website.original_url])
            updated_websites += 1

        self.stdout.write(self.style.SUCCESS(
            f'Successfully added {len(new_pub_mdls)} new publications and updated {len(papers_to_update)} publications with {updated_websites} connections to webpages'))

        # remove websites not linked to any papers
        websites_wo_papers = Website.objects.annotate(p_count=Count("papers")).filter(p_count=0)
        num_websites_wo_papers = websites_wo_papers.count()
        self.stdout.write(
            f"Found {num_websites_wo_papers} websites without associated publication.")
        websites_wo_papers.delete()
        self.stdout.write(self.style.SUCCESS(
            f'Successfully removed {num_websites_wo_papers} websites without associated publication.'))
