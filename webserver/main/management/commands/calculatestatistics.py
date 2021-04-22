from django.core.management.base import BaseCommand
from django.conf import settings
from django.contrib.postgres.aggregates import BoolOr
from django.core.mail import send_mail
from main.models import Website, WebsiteCall, WebsiteStatus, CuratedWebsite

from datetime import timedelta, date, datetime

from collections import defaultdict
import random
from urllib.request import Request, urlopen

from tqdm import tqdm


class Command(BaseCommand):
    help = 'Updates the website statistics for the last 30 days'

    def handle(self, *args, **options):
        max_days = settings.TEMPORAL_INFO_DAYS
        website_updates = []
        latest_time = WebsiteCall.objects.latest("datetime").datetime
        website_ok = WebsiteCall.objects.filter(datetime__date__gte=latest_time-timedelta(days=max_days)).values("website", "datetime").annotate(ok=BoolOr("ok"))

        # update website status according to the last 20 hours (i.e. 2 runs)
        latest_time = WebsiteCall.objects.latest("datetime").datetime

        website_statuses = WebsiteCall.objects.filter(
            datetime__gt=latest_time - timedelta(hours=20)).values("website").annotate(
            final_ok=BoolOr("ok"))

        website2status = {e["website"]: (WebsiteStatus.ONLINE if e["final_ok"] else WebsiteStatus.OFFLINE) for e in website_statuses}

        website_states = defaultdict(dict)
        for w in website_ok:
            website_states[w['website']][w["datetime"].date()] = w["ok"]

        for website in tqdm(Website.objects.all()):
            states = []
            offline = 0
            online = 0
            for day_delta in range(max_days, -1, -1):
                d_date = (latest_time-timedelta(days=day_delta)).date()
                state = website_states[website.id].get(d_date, None)
                states.append(state)
                if state is not None:
                    if state:
                        online += 1
                    else:
                        offline += 1
            if "total_heap_size" in website.calls.latest("datetime").json_data:
                website.last_heap_size = website.calls.latest("datetime").json_data["total_heap_size"]#used_heap_size
            else:
                website.last_heap_size = 0
            website.states = states
            website.percentage = None
            if online > 0 or offline > 0:
                website.percentage = 100 * (online / (online + offline))
            status = website2status.get(website.id, WebsiteStatus.UNKNOWN)
            if status == WebsiteStatus.OFFLINE and any(states[-settings.TEMP_OFFLINE_DAYS-1:]):
                status = WebsiteStatus.TEMP_OFFLINE
            website.status = status
            website_updates.append(website)

        Website.objects.bulk_update(website_updates, ["status", "states", "percentage", "last_heap_size"])

        def sum_digits(n):
            s = 0
            while n:
                s += n % 10
                n //= 10
            return s

        #Automated Email message (offline 3 days and online before)
        def handleAutomatedMessage(website):
            if website.days_reminder > 0:
                remind = True
                for i in range(1, website.days_reminder + 1):
                    if (website.states[len(website.states) - i] is None) or (website.states[len(website.states) - i]):
                        remind = False
                        break
                if remind and (website.states[len(website.states) - (website.days_reminder + 1)]):
                    try:
                        send_mail(
                            f'Your website has been offline for {website.days_reminder} days',
                            f'This is an automated email from aviator to inform you, that your website {website.url} has been offline for {website.days_reminder} days',
                            'no-reply@aviator.ccb.uni-saarland.de',
                            [website.contact_mail],
                            fail_silently=False,
                        )
                    except Exception as e:
                        self.stdout.write(str(e.message))

        today = date.today()
        today_dt = datetime.now()
        website_updates = []
        for website in tqdm(CuratedWebsite.objects.all()):
            first_check = False
            input = random.randint(10, 999)
            result = str(sum_digits(input))
            try:
                req = Request(f"{website.api_url}?input={input}", headers={'User-Agent': 'Aviator/1.0'})
                response = urlopen(req, timeout=30).read().decode()
            except Exception:
                response = "0"
            ok = response == result
            if not (len(website.dates) > 0 and website.dates[len(website.dates) - 1].date() == today):
                first_check = True
                website.dates.append(today_dt)
                if ok:
                    website.states.append(True)
                    website.status = WebsiteStatus.ONLINE
                else:
                    website.states.append(False)
                    website.status = WebsiteStatus.OFFLINE
            else:
                website.dates[len(website.dates) - 1] = today_dt
                website.states[len(website.states) - 1] = (ok or website.states[len(website.states) - 1])
                if website.states[len(website.states) - 1]:
                    website.status = WebsiteStatus.ONLINE
                else:
                    website.status = WebsiteStatus.OFFLINE
            while len(website.states) > max_days + 1:
                website.states.pop(0)
                website.dates.pop(0)
            while len(website.states) < max_days + 1:
                website.states.insert(0, None)
                website.dates.insert(0, (website.dates[len(website.dates) - 1] - timedelta(days=1)).date())
            website_updates.append(website)
            if first_check:
                handleAutomatedMessage(website)
        CuratedWebsite.objects.bulk_update(website_updates, ["dates", "states", "status"])

        website_updates = []
        for website in tqdm(CuratedWebsite.objects.all()):
            states = []
            offline = 0
            online = 0
            for day_delta in range(max_days, -1, -1):
                d_date = (datetime.now() - timedelta(days=day_delta)).date()
                state = None
                for d in range(len(website.dates)):
                    if d_date == website.dates[d].date():
                        state = website.dates[d]
                states.append(state)
                if state is not None:
                    if state:
                        online += 1
                    else:
                        offline += 1
            website.percentage = None
            if online > 0 or offline > 0:
                website.percentage = 100 * (online / (online + offline))
            status = website.status
            if status == WebsiteStatus.OFFLINE and any(states[-settings.TEMP_OFFLINE_DAYS-1:]):
                status = WebsiteStatus.TEMP_OFFLINE
            website.status = status
            website_updates.append(website)
        CuratedWebsite.objects.bulk_update(website_updates, ["status", "percentage"])

        self.stdout.write(self.style.SUCCESS("Successfully calculated website statistics"))

