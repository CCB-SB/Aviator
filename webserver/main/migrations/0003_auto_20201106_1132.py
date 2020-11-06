# Generated by Django 2.2.16 on 2020-11-06 10:32

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0002_auto_20201023_0810'),
    ]

    operations = [
        migrations.AddField(
            model_name='website',
            name='percentage',
            field=models.IntegerField(default=-1),
        ),
        migrations.AddField(
            model_name='website',
            name='states',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.NullBooleanField(), default=[], size=None),
        ),
    ]
