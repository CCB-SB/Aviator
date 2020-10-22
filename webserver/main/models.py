from django.db import models
from django.contrib.postgres.fields import JSONField, ArrayField


# Create your models here.
class Publication(models.Model):
    title = models.TextField()
    abstract = models.TextField()
    authors = models.TextField()
    year = models.SmallIntegerField()
    journal = models.TextField()
    pubmed_id = models.CharField(max_length=50, unique=True, db_index=True)
    url = ArrayField(models.TextField())


class Website(models.Model):
    original_url = models.TextField(db_index=True, unique=True)
    derived_url = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    status = models.NullBooleanField(null=True, default=None)
    papers = models.ManyToManyField(Publication, related_name='websites')
    ip = models.GenericIPAddressField(null=True)
    server = models.CharField(max_length=200)
    analytics = models.CharField(max_length=50)
    timezone = models.CharField(max_length=50)
    framework = models.CharField(max_length=50)
    script = models.CharField(max_length=50)
    certificate_secure = models.BooleanField()
    security_issuer = models.CharField(max_length=100)


class WebsiteCall(models.Model):
    website = models.ForeignKey(Website, related_name='calls', on_delete=models.CASCADE)
    datetime = models.DateTimeField(db_index=True)
    ok = models.BooleanField()
    error = models.CharField(max_length=200)
    msg = models.TextField()
    code = models.IntegerField()
    json_data = JSONField()


class Tag(models.Model):
    name = models.CharField(max_length=50)


class AuthorInfo(models.Model):
    tool_name = models.CharField(max_length=200)
    description = models.TextField()
    webserver = models.CharField(max_length=100)
    website = models.OneToOneField(Website, related_name='authorinfo', on_delete=models.CASCADE)
    tags = models.ManyToManyField(Tag, related_name='authorinfos')
