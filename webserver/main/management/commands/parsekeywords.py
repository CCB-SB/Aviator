from django.core.management.base import BaseCommand, CommandError
from main.models import Tag
import pandas as pd


class Command(BaseCommand):
    help = 'Creates the Tag model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
        except:
            raise CommandError('Please provide an input file')

        tag_tbl = pd.read_csv(csv_file, delimiter='\t')
        new_tags = list()
        for line, values in tag_tbl.iterrows():
            tag = values[0]
            if Tag.objects.filter(name=tag).count() == 0:
                new_tags.append(Tag(name=tag))
        Tag.objects.bulk_create(new_tags)
        self.stdout.write(self.style.SUCCESS(f'Successfully added {len(new_tags)} new tags'))
