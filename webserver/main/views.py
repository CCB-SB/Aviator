from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall, Paper

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
    context['websites'] = Website.objects.all()
    return render(request, 'overview.html', context)

def papers(request):
    context = {}
    context['papers'] = Paper.objects.all()
    return render(request, 'papers.html', context)

def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    context['calls'] = website.calls.all()
    return render(request, 'details.html', context)

def author(request):
    context = {}
    context['websites'] = Website.objects.all()
    return render(request, 'author.html', context)
