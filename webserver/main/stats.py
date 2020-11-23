from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Publication, WebsiteStatus
from collections import Counter
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr
from django.db.models.functions import TruncDate, Cast
from django.db.models import Count, Q, BooleanField
from django.db import models
from cache_memoize import cache_memoize
from django.conf import settings


# cache for 6h
@cache_memoize(settings.CACHE_TIMEOUT)
def get_index_stats():
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Publication.objects.count()
    context['online_count'] = Website.objects.filter(status=WebsiteStatus.ONLINE).count()
    latest_time = WebsiteCall.objects.latest("datetime")
    context['temp_offline_count'] = Website.objects.filter(
        status=WebsiteStatus.TEMP_OFFLINE).count()
    context['offline_count'] = Website.objects.filter(status=WebsiteStatus.OFFLINE).count()
    return context


def get_all_statistics(pub_queryset):
    # website count
    # paper count
    # online count
    # temp offline count
    # offline count
    # temporal online, offline, tmp offline
    # top 10 journals, online, offline, tmp offline
    # per year online, offline, tmp offline
    website_ids = set(e[0] for e in pub_queryset.values_list("websites").distinct())
    websites = [w for w in Website.objects.all().values_list("id", "status", "states") if
                w[0] in website_ids]
    context = dict()
    context['website_count'] = len(websites)
    context['paper_count'] = pub_queryset.count()
    latest_time = WebsiteCall.objects.latest("datetime")
    temp_info_num_days = 15
    stat1_names = []
    for i in range(temp_info_num_days):
        c = latest_time.datetime - timedelta(days=temp_info_num_days - i - 1)
        stat1_names.append("{}.{}.{}".format(c.day, c.month, c.year))
    stat1_online = [0 for _ in range(temp_info_num_days)]
    stat1_offline = [0 for _ in range(temp_info_num_days)]
    stat1_tmp_offline = [0 for _ in range(temp_info_num_days)]
    online_count = 0
    offline_count = 0
    tmp_offline_count = 0
    tmp_offline_websites = set()
    for w_id, w_status, w_states in websites:
        for i in range(temp_info_num_days):
            if temp_info_num_days - i <= len(w_states):
                pos = len(w_states) - temp_info_num_days + i
                if w_states[pos]:
                    stat1_online[i] += 1
                else:
                    srange = w_states[(pos - settings.TEMP_OFFLINE_DAYS):pos + 1]
                    online_in_range = any(e is True for e in srange)
                    if online_in_range and (
                        w_status == WebsiteStatus.OFFLINE or w_status == WebsiteStatus.TEMP_OFFLINE):
                        stat1_tmp_offline[i] += 1
                    elif any(e is False for e in srange):
                        stat1_offline[i] += 1
        if w_status == WebsiteStatus.ONLINE:
            online_count += 1
        elif w_status == WebsiteStatus.TEMP_OFFLINE:
            tmp_offline_count += 1
            tmp_offline_websites.add(w_id)
        elif w_status == WebsiteStatus.OFFLINE:
            offline_count += 1
    context['online_count'] = online_count
    context['offline_count'] = offline_count
    context['temp_offline_count'] = tmp_offline_count
    context['stat1_names'] = json.dumps(stat1_names)
    context['stat1_online'] = json.dumps(stat1_online)
    context['stat1_offline'] = json.dumps(stat1_offline)
    context['stat1_tmp_offline'] = json.dumps(stat1_tmp_offline)
    context["top10_journals"] = list(pub_queryset.values("journal").annotate(
        count=Count("journal")).order_by("-count")[:10])
    journals2count = {e["journal"]: e["count"] for e in context["top10_journals"]}
    top_journals_online = Counter(e["journal"] for e in pub_queryset.filter(
        journal__in=journals2count.keys()).annotate(
        status=BoolOr(Q(websites__status=WebsiteStatus.ONLINE),
                      output_field=BooleanField())).filter(
        status=True).values("journal"))
    tmp_offline_all_journals = Counter(e for p in
                                       Website.objects.filter(id__in=tmp_offline_websites).annotate(
                                           paper=ArrayAgg("papers__journal")).values("paper") for e
                                       in p["paper"])
    tmp_offline_top_journals = {j: tmp_offline_all_journals.get(j, 0) for j in top_journals_online}
    #################
    offline_top_journals = {}
    for j, c in journals2count.items():
        if j in top_journals_online and j in tmp_offline_top_journals:
            offline_top_journals[j] = c - top_journals_online[j] - tmp_offline_top_journals[j]
    #offline_top_journals = {j: (c - top_journals_online[j] - tmp_offline_top_journals[j]) for j, c
     #                       in journals2count.items()}
    #################
    journals = []
    for j, _ in sorted(journals2count.items(), key=lambda e: e[1], reverse=True):
        if j in offline_top_journals:
            journals.append(j)
    #journals = [j for j, _ in sorted(journals2count.items(), key=lambda e: e[1], reverse=True)]
    #################
    context["top10_journals_names"] = json.dumps(journals)
    context["top10_journals_online"] = json.dumps([top_journals_online[j] for j in journals])
    context["top10_journals_tmp_offline"] = json.dumps([tmp_offline_top_journals[j] for j in journals])
    context["top10_journals_offline"] = json.dumps([offline_top_journals[j] for j in journals])
    context["publications_per_year"] = list(pub_queryset.values("year").annotate(
        count=Count("year")).order_by("-count"))
    year2count = {e["year"]: e["count"] for e in context["publications_per_year"]}
    years_online = Counter(e["year"] for e in pub_queryset.annotate(
        status=BoolOr(Q(websites__status=WebsiteStatus.ONLINE),
                      output_field=BooleanField())).filter(
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
