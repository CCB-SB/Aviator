from django.core.management.base import BaseCommand, CommandError
from main.models import Website
from main.models import CuratedWebsite
from main.models import Tag
import pandas as pd


class Command(BaseCommand):
    help = 'Creates the CuratedWebsite model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
        except:
            raise CommandError('Please provide an input folder')

        cw_tbl = pd.read_csv(csv_file, sep='\t', index_col="PMID")
        cw_tbl.index = cw_tbl.index.astype(str)
        cw_dict = cw_tbl.to_dict()
        cw_to_update = list(CuratedWebsite.objects.filter(pubmed_id__in=cw_tbl.index))

        def split_w_nan(s, delimiter):
            if pd.isnull(s):
                return []
            else:
                return s.split(delimiter)

        for cw in cw_to_update:
            cw.title = cw_dict["title"][cw.pubmed_id]
            cw.description = cw_dict["description"][cw.pubmed_id]
            cw.year = cw_dict["year"][cw.pubmed_id]
            cw.authors = cw_dict["authors"][cw.pubmed_id].split(", ")
            cw.journal = cw_dict["journal"][cw.pubmed_id]
            cw.api_url = cw_dict["api_url"][cw.pubmed_id]
            cw.contact_mail = cw_dict["contact_mail"][cw.pubmed_id]
            cw.url = cw_dict["url"][cw.pubmed_id]
            if Website.objects.filter(original_url=cw.url).count() > 0:
                cw.website = Website.objects.filter(original_url=cw.url)[0]
            if Website.objects.filter(derived_url=cw.url).count() > 0:
                cw.website = Website.objects.filter(derived_url=cw.url)[0]
            kwds = split_w_nan(cw_dict["keywords"][cw.pubmed_id], ';')
            for kwd in kwds:
                if Tag.objects.filter(name=kwd).count() > 0:
                    cw.tags.add(Tag.objects.filter(name=kwd)[0])
            #Bulk Update With many to many relation doesn't seem to work
            cw.save()

        #CuratedWebsite.objects.bulk_update(cw_to_update, ["title", "description", "year", "authors", "journal", "api_url", "url", "tags"])#
        cw_to_update_ids = set(cw.pubmed_id for cw in cw_to_update)

        def get_cw_models(tbl):
            for i, row in tbl.iterrows():
                website = None
                if Website.objects.filter(original_url=row["url"]).count() > 0:
                    website = Website.objects.filter(original_url=row["url"])[0]
                if Website.objects.filter(derived_url=row["url"]).count() > 0:
                    website = Website.objects.filter(derived_url=row["url"])[0]
                self.stdout.write(str(website))
                if website is not None:
                    cw = CuratedWebsite(pubmed_id=i,
                                      authors=str(row["authors"]).split(", "),
                                      title=row["title"],
                                      description=row["description"],
                                      year=row["year"],
                                      journal=row["journal"],
                                      url=row["url"],
                                      api_url=row["api_url"],
                                      contact_mail=row["contact_mail"],
                                      website=website)
                    kwds = row["keywords"].split(';')
                    cw.save()
                    for kwd in kwds:
                        if Tag.objects.filter(name=kwd).count() > 0:
                            cw.tags.add(Tag.objects.filter(name=kwd)[0])
                    cw.save()
                    yield cw

        new_cws = cw_tbl[~cw_tbl.index.isin(cw_to_update_ids)]
        new_cw_mdls = list(get_cw_models(new_cws))
        #CuratedWebsite.objects.bulk_create(new_cw_mdls)

        self.stdout.write(self.style.SUCCESS(
            f'Successfully added {len(new_cw_mdls)} new curated websites and updated {len(cw_to_update)} curated websites'))
