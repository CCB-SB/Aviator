from django.db import models
from django.contrib.postgres.fields import JSONField, ArrayField
from djchoices import DjangoChoices, ChoiceItem

#Single Row table
class GlobalStatistics(models.Model):
    data_size = models.BigIntegerField()
    num_calls = models.BigIntegerField()
    weekdays_online = ArrayField(models.IntegerField(), default=list)
    weekdays_offline = ArrayField(models.IntegerField(), default=list)
    recovery_rate = ArrayField(models.IntegerField(), default=list)

class Publication(models.Model):
    title = models.TextField(db_index=True)
    abstract = models.TextField()
    authors = ArrayField(models.TextField(), blank=True, default=list)
    year = models.SmallIntegerField(db_index=True)
    journal = models.TextField(db_index=True)
    pubmed_id = models.CharField(max_length=50, unique=True, db_index=True)
    biotools_id = models.CharField(max_length=250, default="")
    user_kwds = ArrayField(models.TextField(null=True, blank=True), blank=True, default=list)
    mesh_terms = ArrayField(models.TextField(null=True, blank=True), blank=True, default=list)
    contact_mail = ArrayField(models.EmailField(null=True, blank=True), blank=True, default=list)
    url = ArrayField(models.TextField())


class WebsiteStatus(DjangoChoices):
    ONLINE = ChoiceItem("O")
    OFFLINE = ChoiceItem("F")
    TEMP_OFFLINE = ChoiceItem("T")
    UNKNOWN = ChoiceItem("U")


class Website(models.Model):
    original_url = models.TextField(db_index=True, unique=True)
    derived_url = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    status = models.CharField(max_length=1, choices=WebsiteStatus.choices,
                              default=WebsiteStatus.UNKNOWN)
    papers = models.ManyToManyField(Publication, related_name='websites')
    ip = models.GenericIPAddressField(null=True)
    server = models.CharField(max_length=200)
    analytics = models.CharField(max_length=50)
    timezone = models.CharField(max_length=50)
    framework = models.CharField(max_length=50)
    script = models.CharField(max_length=50)
    certificate_secure = models.BooleanField()
    security_issuer = models.CharField(max_length=100)
    # store offline/online/NA per day in ascending date order
    states = ArrayField(models.BooleanField(null=True), default=list)
    percentage = models.FloatField(null=True, default=None)
    last_heap_size = models.IntegerField(default=0)


class WebsiteCall(models.Model):
    website = models.ForeignKey(Website, related_name='calls', on_delete=models.CASCADE)
    datetime = models.DateTimeField(db_index=True)
    ok = models.BooleanField()
    error = models.CharField(max_length=200)
    msg = models.TextField()
    code = models.IntegerField()
    json_data = JSONField()

    class Meta:
        indexes = [
            models.Index(fields=['datetime', 'website']),
        ]


class Tag(models.Model):
    name = models.CharField(max_length=50)


class CuratedWebsite(models.Model):
    title = models.CharField(max_length=200)
    pubmed_id = models.CharField(max_length=50, unique=True, db_index=True)
    year = models.SmallIntegerField(db_index=True)
    journal = models.TextField(db_index=True)
    description = models.TextField()
    authors = ArrayField(models.TextField(), blank=True, default=list)
    url = models.CharField(max_length=100)
    api_url = models.CharField(max_length=100)
    tags = models.ManyToManyField(Tag, related_name='curated')
    website = models.OneToOneField(Website, related_name='curated', on_delete=models.CASCADE)
    status = models.CharField(max_length=1, choices=WebsiteStatus.choices,
                              default=WebsiteStatus.UNKNOWN)
    dates = ArrayField(models.DateTimeField(), default=list)
    # store offline/online/NA per day in ascending date order
    states = ArrayField(models.BooleanField(null=True), default=list)
    percentage = models.FloatField(null=True, default=None)
    contact_mail = models.EmailField(null=True, blank=True)
    days_reminder = models.IntegerField(default=0)
