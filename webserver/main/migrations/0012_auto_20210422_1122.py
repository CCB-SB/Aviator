# Generated by Django 3.2 on 2021-04-22 09:22

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0011_auto_20210210_1851'),
    ]

    operations = [
        migrations.AlterField(
            model_name='curatedwebsite',
            name='states',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.BooleanField(null=True), default=list, size=None),
        ),
        migrations.AlterField(
            model_name='website',
            name='states',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.BooleanField(null=True), default=list, size=None),
        ),
    ]
