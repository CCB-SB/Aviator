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

        cw_tbl = pd.read_csv(csv_file, sep='\t', index_col="PubMed ID")
        cw_tbl.index = cw_tbl.index.astype(str)
        cw_dict = cw_tbl.to_dict()
        cw_to_update = list(CuratedWebsite.objects.filter(pubmed_id__in=cw_tbl.index))

        def split_w_nan(s, delimiter):
            if pd.isnull(s):
                return []
            else:
                return s.split(delimiter)

        for cw in cw_to_update:
            cw.title = cw_dict["Tool"][cw.pubmed_id]
            cw.description = cw_dict["Description"][cw.pubmed_id]
            cw.year = cw_dict["year"][cw.pubmed_id]
            cw.authors = cw_dict["Authors"][cw.pubmed_id].split(", ")
            cw.journal = cw_dict["journal"][cw.pubmed_id]
            cw.api_url = cw_dict["API"][cw.pubmed_id]
            cw.contact_mail = cw_dict["Email"][cw.pubmed_id]
            cw.url = cw_dict["URL"][cw.pubmed_id]
            cw.days_reminder = cw_dict["Days before reminder"][cw.pubmed_id]
            other_url = cw.url.replace("http:", "https:") if cw.url.startswith(
                "http:") else cw.url.replace("https:", "http:")
            if Website.objects.filter(original_url=cw.url).exists():
                cw.website = Website.objects.get(original_url=cw.url)
            elif Website.objects.filter(original_url=other_url).exists():
                cw.website = Website.objects.get(original_url=other_url)
            elif Website.objects.filter(derived_url=cw.url).exists():
                cw.website = Website.objects.get(derived_url=cw.url)
            elif Website.objects.filter(derived_url=other_url).exists():
                cw.website = Website.objects.get(derived_url=other_url)
            kwds = split_w_nan(cw_dict["Tags"][cw.pubmed_id], ', ')
            for kwd in kwds:
                if Tag.objects.filter(name=kwd).exists():
                    cw.tags.add(Tag.objects.get(name=kwd))
            # Bulk Update With many to many relation doesn't seem to work
            cw.save()

        cw_to_update_ids = set(cw.pubmed_id for cw in cw_to_update)

        def get_cw_models(tbl):
            for i, row in tbl.iterrows():
                website = None
                other_url = row["URL"].replace("http:", "https:") if row["URL"].startswith(
                    "http:") else row["URL"].replace("https:", "http:")
                if Website.objects.filter(original_url=row["URL"]).exists():
                    website = Website.objects.get(original_url=row["URL"])
                if Website.objects.filter(derived_url=row["URL"]).exists():
                    website = Website.objects.get(derived_url=row["URL"])
                if Website.objects.filter(original_url=other_url).exists():
                    website = Website.objects.get(original_url=other_url)
                if Website.objects.filter(derived_url=other_url).exists():
                    website = Website.objects.get(derived_url=other_url)

                if website is not None:
                    cw = CuratedWebsite(pubmed_id=i,
                                        authors=str(row["Authors"]).split(", "),
                                        title=row["Tool"],
                                        description=row["Abstract"],
                                        year=row["year"],
                                        journal=row["journal"],
                                        url=row["URL"],
                                        api_url=row["API"],
                                        contact_mail=row["Email"],
                                        days_reminder=row["Days before reminder"],
                                        website=website)
                    kwds = row["Tags"].split(';')
                    cw.save()
                    for kwd in kwds:
                        if Tag.objects.filter(name=kwd).count() > 0:
                            cw.tags.add(Tag.objects.filter(name=kwd)[0])
                    cw.save()
                    yield cw

        new_cws = cw_tbl[~cw_tbl.index.isin(cw_to_update_ids)]
        new_cw_mdls = list(get_cw_models(new_cws))
        # CuratedWebsite.objects.bulk_create(new_cw_mdls)

        self.stdout.write(self.style.SUCCESS(
            f'Successfully added {len(new_cw_mdls)} new curated websites and updated {len(cw_to_update)} curated websites'))
