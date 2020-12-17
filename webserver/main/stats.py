from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Publication, WebsiteStatus
from collections import Counter, namedtuple
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
    publication_ids = set(pub_queryset.values_list("id", flat=True))
    if curated:
        CWebsite = namedtuple("CWebsite", ["id", "status", "states", "journal", "year"])
        websites = [CWebsite(*e) for e in
                    pub_queryset.values_list("website__id", "status", "states", "journal", "year")]
        website_papers = {w: [(w.journal, w.year)] for w in websites}
    else:
        websites = [w for w in
                    Website.objects.filter(papers__in=pub_queryset).distinct().prefetch_related(
                        "papers")]
        website_papers = {
            w: [(p.journal, p.year) for p in w.papers.all() if p.id in publication_ids] for w in
            websites}

    context = dict()
    context['website_count'] = len(websites)
    context['paper_count'] = len(publication_ids)
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
    for w in websites:
        for i in range(temp_info_num_days):
            if temp_info_num_days - i <= len(w.states):
                pos = len(w.states) - temp_info_num_days + i
                if w.states[pos]:
                    stat1_online[i] += 1
                else:
                    srange = w.states[(pos - settings.TEMP_OFFLINE_DAYS):pos + 1]
                    online_in_range = any(e is True for e in srange)
                    if online_in_range and (
                        w.status == WebsiteStatus.OFFLINE or w.status == WebsiteStatus.TEMP_OFFLINE):
                        stat1_tmp_offline[i] += 1
                    elif any(e is False for e in srange):
                        stat1_offline[i] += 1
        if w.status == WebsiteStatus.ONLINE:
            online_websites.add(w)
        elif w.status == WebsiteStatus.TEMP_OFFLINE:
            tmp_offline_websites.add(w)
        elif w.status == WebsiteStatus.OFFLINE:
            offline_websites.add(w)

    # website_states = list(qs.values('pubmed_id', 'websites__states'))
    # latest_date = WebsiteCall.objects.latest("datetime").datetime
    #
    # state_dates = [(latest_date - timedelta(days=day_delta)).date().strftime("%Y-%m-%d") for
    #                day_delta in range(settings.TEMPORAL_INFO_DAYS, -1, -1)]

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
