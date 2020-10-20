from django.urls import path
from . import views

app_name = 'main'
urlpatterns = [
    path('', views.index, name='index'),
    path('details/<int:pk>/', views.details, name="details"),
    path('overview', views.overview, name='overview'),
    path('publications', views.publications, name='publications'),
    path('publication/<int:pk>/', views.publication, name='publication'),
    path('author', views.author, name='author'),
    path('paperData', views.paperData, name='paperData'),
    path('websiteData', views.websiteData, name='websiteData'),
]
