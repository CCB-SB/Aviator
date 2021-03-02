from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Publication, WebsiteStatus, CuratedWebsite
from collections import Counter, namedtuple, defaultdict
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
    context['temp_offline_count'] = Website.objects.filter(
        status=WebsiteStatus.TEMP_OFFLINE).count()
    context['offline_count'] = Website.objects.filter(status=WebsiteStatus.OFFLINE).count()
    return context


def get_all_statistics(pub_queryset, curated=False):
    # website count
    # paper count
    # online count
    # temp offline count
    # offline count
    # temporal online, offline, tmp offline
    # top 10 journals, online, offline, tmp offline
    # per year online, offline, tmp offline
    if curated:
        p_websites = pub_queryset.values_list("journal", "year", "website__id", "states", "status",
                                              "pubmed_id")
        num_pubs = len(p_websites)
        website_papers = defaultdict(list)
        website2states = dict()
        for e in p_websites:
            website_papers[e[2]].append((e[0], e[1]))
            website2states[e[2]] = (e[3], e[4], e[5])
    else:
        p_websites = pub_queryset.annotate(website_states=ArrayAgg("websites__states")).values_list(
            "journal", "year", "website_pks", "website_states", "status", "pubmed_id")
        num_pubs = len(p_websites)
        website_papers = defaultdict(list)
        website2states = dict()
        for e in p_websites:
            for w_id, w_states, w_status in zip(e[2], e[3], e[4]):
                website_papers[w_id].append((e[0], e[1]))
                website2states[w_id] = (w_states, w_status, e[5])

    context = dict()
    context['website_count'] = len(website2states)
    context['paper_count'] = num_pubs
    latest_time = WebsiteCall.objects.latest("datetime")
    temp_info_num_days = 15
    stat1_names = []
    for i in range(temp_info_num_days):
        c = latest_time.datetime - timedelta(days=temp_info_num_days - i - 1)
        stat1_names.append("{}.{}.{}".format(c.day, c.month, c.year))
    stat1_online = [0 for _ in range(temp_info_num_days)]
    stat1_offline = [0 for _ in range(temp_info_num_days)]
    stat1_tmp_offline = [0 for _ in range(temp_info_num_days)]

    tmp_offline_websites = set()
    online_websites = set()
    offline_websites = set()
    for w_id, (w_states, w_status, _) in website2states.items():
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
            online_websites.add(w_id)
        elif w_status == WebsiteStatus.TEMP_OFFLINE:
            tmp_offline_websites.add(w_id)
        elif w_status == WebsiteStatus.OFFLINE:
            offline_websites.add(w_id)

    latest_date = latest_time.datetime
    if curated:
        website_states = [{"pubmed_id": e[5], "websites__states": e[3]} for e in p_websites]
        if CuratedWebsite.objects.exists():
            dates = CuratedWebsite.objects.first().dates
            latest_date = dates[len(dates) - 1]
    else:
        website_states = [{"pubmed_id": e[5], "websites__states": w_states} for e in p_websites for
                          w_states in e[3]]

    state_dates = [(latest_date - timedelta(days=day_delta)).date().strftime("%Y-%m-%d")
                   for
                   day_delta in range(settings.TEMPORAL_INFO_DAYS, -1, -1)]

    context["website_states"] = website_states
    context["state_dates"] = state_dates

    context['online_count'] = len(online_websites)
    context['offline_count'] = len(offline_websites)
    context['temp_offline_count'] = len(tmp_offline_websites)
    context['stat1_names'] = json.dumps(stat1_names)
    context['stat1_online'] = json.dumps(stat1_online)
    context['stat1_offline'] = json.dumps(stat1_offline)
    context['stat1_tmp_offline'] = json.dumps(stat1_tmp_offline)
    journal_counts = Counter(e[0] for papers in website_papers.values() for e in papers)
    context["top10_journals"] = [{"journal": e[0], "count": e[1]} for e in
                                 journal_counts.most_common(10)]
    journals2count = {e["journal"]: e["count"] for e in context["top10_journals"]}

    top_journals_online = Counter(
        e[0] for w in online_websites for e in website_papers[w] if e[0] in journals2count)
    tmp_offline_top_journals = Counter(
        e[0] for w in tmp_offline_websites for e in website_papers[w] if e[0] in journals2count)
    offline_top_journals = Counter(
        e[0] for w in offline_websites for e in website_papers[w] if e[0] in journals2count)

    journals = [j for j, _ in sorted(journals2count.items(), key=lambda e: e[1], reverse=True)]
    context["top10_journals_names"] = json.dumps(journals)
    context["top10_journals_online"] = json.dumps([top_journals_online.get(j, 0) for j in journals])
    context["top10_journals_tmp_offline"] = json.dumps(
        [tmp_offline_top_journals.get(j, 0) for j in journals])
    context["top10_journals_offline"] = json.dumps(
        [offline_top_journals.get(j, 0) for j in journals])

    year_counts = Counter(e[1] for papers in website_papers.values() for e in papers)
    context["publications_per_year"] = [{"year": e[0], "count": e[1]} for e in
                                        sorted(year_counts.items(), key=lambda e: e[0])]
    year2count = {e["year"]: e["count"] for e in context["publications_per_year"]}
    years_online = Counter(e[1] for w in online_websites for e in website_papers[w])
    tmp_offline_years = Counter(e[1] for w in tmp_offline_websites for e in website_papers[w])
    offline_years = Counter(e[1] for w in offline_websites for e in website_papers[w])
    years = sorted(year2count.keys())

    context["pubs_per_year_names"] = json.dumps(years)
    context["pubs_per_year_online"] = json.dumps([years_online[y] for y in years])
    context["pubs_per_year_tmp_offline"] = json.dumps([tmp_offline_years[y] for y in years])
    context["pubs_per_year_offline"] = json.dumps([offline_years[y] for y in years])
    return context
