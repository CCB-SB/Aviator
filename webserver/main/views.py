from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Paper
from django.core import serializers
from django.core.paginator import Paginator

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
    context = {}
    context['websites'] = Website.objects.all().prefetch_related('papers')
    return render(request, 'overview.html', context)

def publications(request):
    context = {}
    #context['papers'] = serializers.serialize('json', Paper.objects.all().prefetch_related('websites'))
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
    return JsonResponse({"data": list(Website.objects.all().values('url', 'status', 'created_at', 'updated_at', 'pk', 'papers'))})#, 'ip', 'certificate_secure'

def paperData(request):
    search_values = []
    #fields = ['title', 'authors', 'year', 'journal', 'pubmed_id']
    #papers = Paper.objects.all()#filter(reduce(AND, (Q(**{fields[i]+'__icontains': value} ) for i, value in enumerate(search_values)))).values('title', 'authors', 'year', 'journal', 'pubmed_id')

    return JsonResponse({"data": list(Paper.objects.all().values('title', 'authors', 'abstract', 'year', 'journal', 'pubmed_id', 'websites'))})
    #paginator = Paginator(Paper.objects.all().values('title', 'authors', 'year', 'journal', 'pubmed_id'), 25)
    #page = request.GET.get('page')
    #context = {"data": list(paginator.get_page(page).object_list)}
    #context["draw"] = paginator.num_pages
    #context["recordsTotal"] = paginator.count
    #context["recordsFiltered"] = paginator.count
    #return JsonResponse(context)
