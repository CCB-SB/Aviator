from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.cache import cache_page

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

from .stats import get_index_stats, get_all_statistics

def index(request):
    context = get_index_stats()
    return render(request, 'index.html', context)

def overview(request):
    context = {'search_column':0, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'overview.html', context)

@cache_memoize(settings.CACHE_TIMEOUT)
def get_publication_datatable_info():
    context = {}
    context["websites"] = json.dumps({x['pk']: x for x in list(
        Website.objects.all().values('pk', 'status', 'original_url', 'derived_url').annotate(
            calls=ArrayAgg('calls')))})
    context["calls"] = json.dumps({x['pk']: x for x in list(
        WebsiteCall.objects.filter(datetime__gt=(date.today() - timedelta(days=140))).values('pk',
                                                                                             'website',
                                                                                             'ok',
                                                                                             'error',
                                                                                             'code').annotate(
            datetime=Cast(TruncDate('datetime'), models.CharField())))})
    return context

def publications(request):
    context = {'search_column': -1, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    context.update(get_publication_datatable_info())
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

# cache for 6 hours
@cache_page(60 * 60 * 6)
def websiteData(request):
    return JsonResponse({"data": list(Website.objects.all().values('original_url', 'derived_url', 'status', 'created_at', 'updated_at', 'pk', 'papers'))})

# cache for 6 hours
@cache_page(60 * 60 * 6)
def paperData(request):
    data_papers = list(Publication.objects.all().values('pk', 'title', 'url', 'authors', 'abstract', 'year', 'journal', 'pubmed_id', 'contact_mail', 'user_kwds').annotate(websites=ArrayAgg('websites')))
    return JsonResponse({"data": data_papers})


def statistics(request):
    context = get_all_statistics()
    return render(request, 'statistics.html', context)
