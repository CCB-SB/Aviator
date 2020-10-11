from django.core.management.base import BaseCommand, CommandError
from main.models import WebsiteCall
from datetime import datetime, timedelta

class Command(BaseCommand):
    help = 'Removes data older than 14 days'

    def handle(self, *args, **options):
        max_days = 14
        WebsiteCall.objects.filter(datetime__lte=datetime.now()-timedelta(days=max_days)).delete()
        self.stdout.write(self.style.SUCCESS('Successfully removed old Website Calls'))
