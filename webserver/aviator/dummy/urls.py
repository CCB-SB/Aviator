from django.urls import path

from aviator.dummy.views import index

app_name = "dummy"
urlpatterns = [
    path("", view=index, name="index"),
]
