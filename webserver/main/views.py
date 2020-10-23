from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Publication
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr
from django.db.models.functions import TruncDate, Cast
from django.db.models import Count
from django.db import models

# Create your views here.
def index(request):
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Publication.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    context['offline_count'] = Website.objects.count() - online.count()
    return render(request, 'index.html', context)

def overview(request):
    context = {'search_column':0, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'overview.html', context)

def publications(request):
    context = {'search_column': -1, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    context["websites"] = json.dumps({x['pk']:x for x in list(Website.objects.all().values('pk', 'status', 'original_url', 'derived_url').annotate(calls=ArrayAgg('calls')))})
    context["calls"] = json.dumps({x['pk']:x for x in list(WebsiteCall.objects.filter(datetime__gt=(date.today() - timedelta(days=140))).values('pk', 'website', 'ok', 'error', 'code').annotate(datetime=Cast(TruncDate('datetime'), models.CharField())))})
    return render(request, 'publications.html', context)

def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    context['calls'] = website.calls.all().order_by('datetime')
    context['website'] = website
    return render(request, 'details.html', context)

def publication(request, pk):
    context = {}
    paper = get_object_or_404(Publication, pk=pk)
    context['paper'] = paper
    context['websites'] = paper.websites.all()
    return render(request, 'publication.html', context)

def author(request):
    context = {}
    context['websites'] = Website.objects.all()
    return render(request, 'author.html', context)

def websiteData(request):
    return JsonResponse({"data": list(Website.objects.all().values('original_url', 'derived_url', 'status', 'created_at', 'updated_at', 'pk', 'papers'))})

def paperData(request):
    data_papers = list(Publication.objects.all().values('pk', 'title', 'url', 'authors', 'abstract', 'year', 'journal', 'pubmed_id', 'contact_mail', 'user_kwds').annotate(websites=ArrayAgg('websites')))
    return JsonResponse({"data": data_papers})


def statistics(request):
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Publication.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    latest_time = WebsiteCall.objects.latest("datetime")
    context['temp_offline_count'] = WebsiteCall.objects.filter(
        datetime__gt=latest_time.datetime - timedelta(days=30)).values("website").annotate(
        final_ok=BoolOr("ok")).filter(final_ok=True, website__status=False).count()
    context['offline_count'] = Website.objects.count() - context['online_count'] - context['temp_offline_count']

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
        stat_names += '"' + str(current_date.day) + "." + str(current_date.month) + "." + str(current_date.year) + '"'
        ws_online = WebsiteCall.objects.filter(datetime__date=current_date).values("website").annotate(final_ok=BoolOr("ok")).filter(final_ok=True).count()
        stat_online += f'"{ws_online}"'
        ws_tmp_offline = WebsiteCall.objects.filter(
            datetime__range=(current_date - timedelta(days=30), current_date)).values("website").annotate(
            final_ok=BoolOr("ok")).filter(final_ok=True, website__status=False).count()
        stat_tmp_offline += f'"{ws_tmp_offline}"'
        ws_offline = WebsiteCall.objects.filter(
            datetime__range=(current_date - timedelta(days=30), current_date)).values("website").annotate(
            final_ok=BoolOr("ok")).filter(final_ok=False, website__status=False).count()
        stat_offline += f'"{ws_offline}"'
    context['stat1_names'] = stat_names
    context['stat1_online'] = stat_online
    context['stat1_offline'] = stat_offline
    context['stat1_tmp_offline'] = stat_tmp_offline

    context["top10_journals"] = Publication.objects.values("journal").annotate(count=Count("journal")).order_by("-count")[:10]
    journals = [e["journal"] for e in context["top10_journals"]]



    return render(request, 'statistics.html', context)
