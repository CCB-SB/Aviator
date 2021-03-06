# Generated by Django 2.2.16 on 2020-12-07 15:05

import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0002_website_last_heap_size'),
    ]

    operations = [
        migrations.CreateModel(
            name='CuratedWebsite',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tool_name', models.CharField(max_length=200)),
                ('pubmed_id', models.CharField(db_index=True, max_length=50, unique=True)),
                ('description', models.TextField()),
                ('authors', models.TextField()),
                ('url', models.CharField(max_length=100)),
                ('api_url', models.CharField(max_length=100)),
                ('status', models.CharField(choices=[('O', 'ONLINE'), ('F', 'OFFLINE'), ('T', 'TEMP OFFLINE'), ('U', 'UNKNOWN')], default='U', max_length=1)),
                ('states', django.contrib.postgres.fields.ArrayField(base_field=models.NullBooleanField(), default=list, size=None)),
                ('percentage', models.FloatField(default=None, null=True)),
                ('tags', models.ManyToManyField(related_name='authorinfos', to='main.Tag')),
                ('website', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='authorinfo', to='main.Website')),
            ],
        ),
        migrations.DeleteModel(
            name='AuthorInfo',
        ),
    ]
