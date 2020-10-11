from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render
from .models import Website

# Create your views here.
def index(request):
    context = {}
    websites = Website.objects.all()
    context['websites'] = websites
    return render(request, 'index.html', context)

