from django.core.management.base import BaseCommand, CommandError

from main.models import Publication
from collections import defaultdict
import pandas as pd


class Command(BaseCommand):
    help = 'Adds the bio.tools to the publications model'

    def add_arguments(self, parser):
        parser.add_argument('biotoolsid2pmid', default=None)

    def handle(self, *args, **options):
        try:
            biotoolsid2pmid = options['biotoolsid2pmid']
        except:
            raise CommandError('Please provide an biotoolsid2pmid csv')


        df = pd.read_csv(biotoolsid2pmid, sep='\t', names = ["biotoolsID","pmid"])
        pmid2biotoolsid = dict(zip(list(df.pmid), list(df.biotoolsID)))

        to_update = []
        for pub in Publication.objects.all():
            if pmid2biotoolsid:
                if pub.pubmed_id in pmid2biotoolsid:
                    pub.biotools_id = pmid2biotoolsid[pub.pubmed_id]
                    to_update.append(pub)

        Publication.objects.bulk_update(to_update, ["biotools_id"])

        self.stdout.write(self.style.SUCCESS(f'Updated {len(to_update)} publications'))