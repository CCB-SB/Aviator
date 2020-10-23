from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Publication
from collections import Counter
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr
from django.db.models.functions import TruncDate, Cast
from django.db.models import Count
from django.db import models
from cache_memoize import cache_memoize
from django.conf import settings

# cache for 6h
@cache_memoize(settings.CACHE_TIMEOUT)
def get_index_stats():
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Publication.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    latest_time = WebsiteCall.objects.latest("datetime")
    tmp_offline = WebsiteCall.objects.filter(
        datetime__gt=latest_time.datetime - timedelta(days=30)).values("website").annotate(
        final_ok=BoolOr("ok")).filter(final_ok=True, website__status=False)
    context['temp_offline_count'] = tmp_offline.count()
    context['offline_count'] = Website.objects.count() - context['online_count'] - context[
        'temp_offline_count']
    return context

@cache_memoize(settings.CACHE_TIMEOUT)
def get_all_statistics():
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Publication.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    latest_time = WebsiteCall.objects.latest("datetime")
    tmp_offline = WebsiteCall.objects.filter(
        datetime__gt=latest_time.datetime - timedelta(days=30)).values("website").annotate(
        final_ok=BoolOr("ok")).filter(final_ok=True, website__status=False)
    tmp_offline_websites = set(e["website"] for e in tmp_offline)
    context['temp_offline_count'] = tmp_offline.count()
    context['offline_count'] = Website.objects.count() - context['online_count'] - context[
        'temp_offline_count']

    current_date = latest_time.datetime - timedelta(days=14)
    stat_names = ""
    stat_online = ""
    stat_offline = ""
    stat_tmp_offline = ""
    for i in range(0, 14):
        current_date = current_date + timedelta(days=1)
        if i > 0:
            stat_online += ", "
            stat_offline += ", "
            stat_tmp_offline += ", "
            stat_names += ", "
        stat_names += '"' + str(current_date.day) + "." + str(current_date.month) + "." + str(
            current_date.year) + '"'
        ws_online = WebsiteCall.objects.filter(datetime__date=current_date).values(
            "website").annotate(final_ok=BoolOr("ok")).filter(final_ok=True).count()
        stat_online += f'"{ws_online}"'
        ws_tmp_offline = WebsiteCall.objects.filter(
            datetime__range=(current_date - timedelta(days=30), current_date)).values(
            "website").annotate(
            final_ok=BoolOr("ok")).filter(final_ok=True, website__status=False).count()
        stat_tmp_offline += f'"{ws_tmp_offline}"'
        ws_offline = WebsiteCall.objects.filter(
            datetime__range=(current_date - timedelta(days=30), current_date)).values(
            "website").annotate(
            final_ok=BoolOr("ok")).filter(final_ok=False, website__status=False).count()
        stat_offline += f'"{ws_offline}"'
    context['stat1_names'] = stat_names
    context['stat1_online'] = stat_online
    context['stat1_offline'] = stat_offline
    context['stat1_tmp_offline'] = stat_tmp_offline

    context["top10_journals"] = Publication.objects.values("journal").annotate(
        count=Count("journal")).order_by("-count")[:10]
    journals2count = {e["journal"]: e["count"] for e in context["top10_journals"]}
    top_journals_online = Counter(e["journal"] for e in Publication.objects.filter(
        journal__in=journals2count.keys()).annotate(status=BoolOr("websites__status")).filter(
        status=True).values("journal"))
    tmp_offline_all_journals = Counter(e for p in
                                       Website.objects.filter(id__in=tmp_offline_websites).annotate(
                                           paper=ArrayAgg("papers__journal")).values("paper") for e
                                       in p["paper"])
    tmp_offline_top_journals = {j: tmp_offline_all_journals.get(j, 0) for j in top_journals_online}
    offline_top_journals = {j: (c - top_journals_online[j] - tmp_offline_top_journals[j]) for j, c
                            in journals2count.items()}
    journals = [j for j, _ in sorted(journals2count.items(), key=lambda e: e[1], reverse=True)]
    context["top10_journals_names"] = json.dumps(journals)
    context["top10_journals_online"] = json.dumps([top_journals_online[j] for j in journals])
    context["top10_journals_tmp_offline"] = json.dumps(
        [tmp_offline_top_journals[j] for j in journals])
    context["top10_journals_offline"] = json.dumps([offline_top_journals[j] for j in journals])

    context["publications_per_year"] = Publication.objects.values("year").annotate(
        count=Count("year")).order_by("-count")
    year2count = {e["year"]: e["count"] for e in context["publications_per_year"]}
    years_online = Counter(e["year"] for e in Publication.objects.all().annotate(
        status=BoolOr("websites__status")).filter(
        status=True).values("year"))
    tmp_offline_years = Counter(e for p in
                                Website.objects.filter(id__in=tmp_offline_websites).annotate(
                                    paper=ArrayAgg("papers__year")).values("paper") for e
                                in p["paper"])
    offline_years = {j: (c - years_online[j] - tmp_offline_years[j]) for j, c
                     in year2count.items()}
    years = sorted(year2count.keys())
    context["pubs_per_year_names"] = json.dumps(years)
    context["pubs_per_year_online"] = json.dumps([years_online[y] for y in years])
    context["pubs_per_year_tmp_offline"] = json.dumps([tmp_offline_years[y] for y in years])
    context["pubs_per_year_offline"] = json.dumps([offline_years[y] for y in years])
    return context
