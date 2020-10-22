# Generated by Django 2.2.16 on 2020-10-22 20:39

import django.contrib.postgres.fields
import django.contrib.postgres.fields.jsonb
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Publication',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.TextField()),
                ('abstract', models.TextField()),
                ('authors', models.TextField()),
                ('year', models.SmallIntegerField()),
                ('journal', models.TextField()),
                ('pubmed_id', models.CharField(db_index=True, max_length=50, unique=True)),
                ('url', django.contrib.postgres.fields.ArrayField(base_field=models.TextField(), size=None)),
            ],
        ),
        migrations.CreateModel(
            name='Tag',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50)),
            ],
        ),
        migrations.CreateModel(
            name='Website',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('original_url', models.TextField(db_index=True, unique=True)),
                ('derived_url', models.TextField()),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('status', models.NullBooleanField(default=None)),
                ('ip', models.GenericIPAddressField(null=True)),
                ('server', models.CharField(max_length=200)),
                ('analytics', models.CharField(max_length=50)),
                ('timezone', models.CharField(max_length=50)),
                ('framework', models.CharField(max_length=50)),
                ('script', models.CharField(max_length=50)),
                ('certificate_secure', models.BooleanField()),
                ('security_issuer', models.CharField(max_length=100)),
                ('papers', models.ManyToManyField(related_name='websites', to='main.Publication')),
            ],
        ),
        migrations.CreateModel(
            name='WebsiteCall',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('datetime', models.DateTimeField(db_index=True)),
                ('ok', models.BooleanField()),
                ('error', models.CharField(max_length=200)),
                ('msg', models.TextField()),
                ('code', models.IntegerField()),
                ('json_data', django.contrib.postgres.fields.jsonb.JSONField()),
                ('website', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='calls', to='main.Website')),
            ],
        ),
        migrations.CreateModel(
            name='AuthorInfo',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tool_name', models.CharField(max_length=200)),
                ('description', models.TextField()),
                ('webserver', models.CharField(max_length=100)),
                ('tags', models.ManyToManyField(related_name='authorinfos', to='main.Tag')),
                ('website', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='authorinfo', to='main.Website')),
            ],
        ),
    ]
