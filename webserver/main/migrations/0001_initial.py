# Generated by Django 2.1.15 on 2020-10-19 08:41

import django.contrib.postgres.fields.jsonb
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='AuthorInfo',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tool_name', models.CharField(max_length=200)),
                ('description', models.TextField()),
                ('webserver', models.CharField(max_length=100)),
            ],
        ),
        migrations.CreateModel(
            name='Paper',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.TextField()),
                ('abstract', models.TextField()),
                ('authors', models.TextField()),
                ('year', models.CharField(max_length=4)),
                ('journal', models.TextField()),
                ('pubmed_id', models.CharField(max_length=50)),
                ('url', models.TextField()),
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
                ('url', models.TextField()),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('updated_at', models.DateTimeField(auto_now=True)),
                ('status', models.NullBooleanField(default=None)),
                ('ip', models.CharField(max_length=45)),
                ('server', models.CharField(max_length=200)),
                ('analytics', models.CharField(max_length=50)),
                ('timezone', models.CharField(max_length=50)),
                ('framework', models.CharField(max_length=50)),
                ('script', models.CharField(max_length=50)),
                ('certificate_secure', models.BooleanField()),
                ('security_issuer', models.CharField(max_length=100)),
                ('papers', models.ManyToManyField(related_name='websites', to='main.Paper')),
            ],
        ),
        migrations.CreateModel(
            name='WebsiteCall',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('datetime', models.DateTimeField()),
                ('ok', models.BooleanField()),
                ('error', models.CharField(max_length=200)),
                ('msg', models.TextField()),
                ('code', models.IntegerField()),
                ('json_data', django.contrib.postgres.fields.jsonb.JSONField()),
                ('website', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='calls', to='main.Website')),
            ],
        ),
        migrations.AddField(
            model_name='authorinfo',
            name='tags',
            field=models.ManyToManyField(related_name='authorinfos', to='main.Tag'),
        ),
        migrations.AddField(
            model_name='authorinfo',
            name='website',
            field=models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='authorinfo', to='main.Website'),
        ),
    ]
