from django.urls import include, path
from . import views

app_name = 'main'
urlpatterns = [
    path('', views.index, name='index'),
    path('help', views.help, name='help'),
    path('details/<int:pk>/', views.details, name="details"),
    path('overview', views.overview, name='overview'),
    path('publications', views.publications, name='publications'),
    path('publication/<int:pk>/', views.publication, name='publication'),
    path('author', views.author, name='author'),
    path('paperData', views.paperData, name='paperData'),
    path('websiteData', views.websiteData, name='websiteData'),
    path('statistics', views.statistics, name='statistics'),
    path("table_data", view=views.Table.as_view(), name='table_data'),
    path("curated_data", view=views.CuratedTable.as_view(), name='curated_data'),
    path('autocomplete', views.autocomplete, name='autocomplete'),
    path('api', views.api, name='api'),
    path('curated', views.curated, name='curated'),
    path('curated_autocomplete', views.curated_autocomplete, name='curated_autocomplete'),
    path('aviator_api', views.aviator_api, name='aviator_api'),
]
