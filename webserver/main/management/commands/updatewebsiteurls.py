from django.core.management.base import BaseCommand, CommandError
from main.models import Website
from main.models import Publication
from collections import defaultdict
import pandas as pd


class Command(BaseCommand):
    help = 'Creates the paper model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
        except:
            raise CommandError('Please provide an input folder')

        orig2new = pd.read_csv(csv_file, sep='\t', index_col="orig").to_dict()["new"]

        orig_urls = {w.original_url for w in Website.objects.all()}
        websites_to_check = Website.objects.filter(original_url__in=orig2new.keys())

        # first https pages, then the others
        websites_to_check = sorted(websites_to_check,
                                   key=lambda w: (w.original_url.startswith("https"), w.original_url),
                                   reverse=True)
        websites_to_update = []
        websites_w_papers_removed = 0
        for website in websites_to_check:
            new_url = orig2new[website.original_url]
            # url was changed
            if website.original_url != new_url:
                # already have a website for the new_url
                if new_url in orig_urls:
                    # remove the associated papers
                    # since they will be linked to the other website entry
                    website.papers.clear()
                    websites_w_papers_removed += 1
                else:
                    print(f"{website.original_url} -> {new_url}")
                    website.original_url = new_url
                    orig_urls.add(new_url)
                    websites_to_update.append(website)

        Website.objects.bulk_update(websites_to_update,
                                    ["original_url"])

        # remove websites not linked to any papers
        websites_wo_papers = Website.objects.annotate(p_count=Count("papers")).filter(p_count=0)
        num_websites_wo_papers = websites_wo_papers.count()
        websites_wo_papers.delete()

        self.stdout.write(self.style.SUCCESS(
            f'Updated {len(websites_to_update)} websites and removed papers from {websites_w_papers_removed} websites. '
            f'Removed {num_websites_wo_papers} websites.'))
