from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from .models import Website, WebsiteCall

# Create your views here.
def index(request):
    context = {}
    websites = Website.objects.all()
    context['websites'] = websites
    return render(request, 'index.html', context)

def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    context['calls'] = website.calls.all()
    return render(request, 'details.html', context)

