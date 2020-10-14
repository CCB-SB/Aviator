from django.urls import path
from . import views

app_name = 'main'
urlpatterns = [
    path('', views.index, name='index'),
    path('details/<int:pk>/', views.details, name="details"),
    path('overview', views.overview, name='overview'),
    path('papers', views.papers, name='papers'),
    path('author', views.author, name='author'),
]
