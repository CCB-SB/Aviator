from django.core.management.base import BaseCommand, CommandError
from main.models import Website, WebsiteCall
from datetime import date, datetime, timedelta
import math
import pytz
from tqdm import tqdm

class Command(BaseCommand):
    help = 'Updates the website statistics for the last 30 days'
    def handle(self, *args, **options):
        max_days = 30
        websites = Website.objects.all().prefetch_related("calls")
        website_updates = []
        now = datetime.utcnow()
        for website in tqdm(websites):
            states = list()
            offline = 0
            online = 0
            latest_state = None
            for day_delta in range(1, max_days+1):
                calls = website.calls.filter(datetime__date=now-timedelta(days=(max_days - day_delta)))
                found = False
                state = False
                for call in calls:
                    state = state | (call.ok & (call.error == "") & (call.code == 200))
                    found = True
                states.append(state)
                if found:
                    if state:
                        online += 1
                        latest_state = True
                    else:
                        offline += 1
                        latest_state = False
            website.states = states
            website.percentage = -1
            if online > 0 or offline > 0:
                website.percentage = math.ceil(100*(online/(online+offline)))
            if latest_state is False and online > 1:
                latest_state = None
            website.status = latest_state
            website_updates.append(website)
        Website.objects.bulk_update(website_updates, ["status", "states", "percentage"])
        self.stdout.write(self.style.SUCCESS("Successfully calculated website statistics"))

