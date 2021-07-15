from django.core.management.base import BaseCommand, CommandError
from django.db.models import Count

from main.models import Website

import pandas as pd
import json

class Command(BaseCommand):
    help = 'Creates the paper model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')
        parser.add_argument('new2orig_prev')

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
            new2orig_file = options['new2orig_prev']
        except:
            raise CommandError('Please provide an input file')

        orig2new = pd.read_csv(csv_file, sep='\t', index_col="orig").to_dict()["new"]

        with open(new2orig_file) as f:
            new2orig_prev = json.load(f)
        orig_urls = {w.original_url:w for w in Website.objects.all()}

        # add previously updated urls to their original representation
        for orig, w in list(orig_urls.items()):
            if orig in new2orig_prev:
                for prev in new2orig_prev[orig]:
                    orig_urls[prev] = w

        websites_to_check = [(o, w) for o, w in orig_urls.items() if o in orig2new]

        # first https pages, then the others
        websites_to_check = sorted(websites_to_check,
                                   key=lambda w: (w[0].startswith("https"), w[0]),
                                   reverse=True)
        websites_to_update = []
        websites_w_papers_removed = 0
        for o, website in websites_to_check:
            new_url = orig2new[o]
            # url was changed
            if website.original_url != new_url:
                # already have a website for the new_url
                if new_url in orig_urls and orig_urls[new_url].original_url == new_url:
                    # remove the associated papers
                    # since they will be linked to the other website entry
                    if website.id != orig_urls[new_url].id:
                        website.papers.clear()
                        websites_w_papers_removed += 1
                else:
                    print(f"{website.original_url} -> {new_url}")
                    website.original_url = new_url
                    orig_urls[new_url] = website
                    websites_to_update.append(website)

        Website.objects.bulk_update(websites_to_update,
                                    ["original_url"])

        # remove websites not linked to any papers
        websites_wo_papers = Website.objects.annotate(p_count=Count("papers")).filter(p_count=0)
        num_websites_wo_papers = websites_wo_papers.count()
        self.stdout.write(f"Found {num_websites_wo_papers} websites without associated publication.")
        websites_wo_papers.delete()

        self.stdout.write(self.style.SUCCESS(
            f'Updated {len(websites_to_update)} websites and removed papers from {websites_w_papers_removed} websites. '
            f'Removed {num_websites_wo_papers} websites.'))
