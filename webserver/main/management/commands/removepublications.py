from django.core.management.base import BaseCommand, CommandError
from main.models import Website, WebsiteCall, Publication

class Command(BaseCommand):
    help = 'Removes all data'

    def handle(self, *args, **options):
        Publication.objects.all().delete()
        self.stdout.write(self.style.SUCCESS('Successfully removed all publications'))
