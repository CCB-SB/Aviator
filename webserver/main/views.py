from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Paper
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime

# Create your views here.
def index(request):
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Paper.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    context['offline_count'] = Website.objects.count() - online.count()
    return render(request, 'index.html', context)

def overview(request):
    context = {'search':''}
    if request.method == 'POST':
        context['search'] = request.POST['search']
    return render(request, 'overview.html', context)

def publications(request):
    context = {}
    return render(request, 'publications.html', context)

def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    context['calls'] = website.calls.all()
    context['website'] = website
    return render(request, 'details.html', context)

def publication(request, pk):
    context = {}
    paper = get_object_or_404(Paper, pk=pk)
    context['paper'] = paper
    context['websites'] = paper.websites.all()
    return render(request, 'publication.html', context)

def author(request):
    context = {}
    context['websites'] = Website.objects.all()
    return render(request, 'author.html', context)

def websiteData(request):
    return JsonResponse({"data": list(Website.objects.all().values('url', 'status', 'created_at', 'updated_at', 'pk', 'papers'))})

def paperData(request):
    return JsonResponse({"data": list(Paper.objects.all().values('title', 'authors', 'abstract', 'year', 'journal', 'pubmed_id', 'websites'))})

def statistics(request):
    context = {}
    context['website_count'] = Website.objects.count()
    context['paper_count'] = Paper.objects.count()
    online = Website.objects.filter(status=True)
    context['online_count'] = online.count()
    context['offline_count'] = Website.objects.count() - online.count()

    current_date = date.today() - timedelta(days=14)
    #current_date = current_date - timedelta(days=40)
    stat_names = ""
    stat_online = ""
    stat_offline = ""
    for i in range(0, 14):
        current_date = current_date + timedelta(days=1)
        if i > 0:
            stat_online += ", "
            stat_offline += ", "
            stat_names += ", "
        stat_names += '"' + str(current_date.day) + "." + str(current_date.month) + "." + str(current_date.year) + '"'
        ws = WebsiteCall.objects.filter(datetime__date=current_date, ok=True, error="", code=200)
        stat_online += '"'+str(ws.count())+'"'
        ws = WebsiteCall.objects.filter(datetime__date=current_date).exclude(ok=True, error="", code=200)
        stat_offline += '"'+str(+ws.count())+'"'
    context['stat1_names'] = stat_names
    context['stat1_online'] = stat_online
    context['stat1_offline'] = stat_offline
    return render(request, 'statistics.html', context)
