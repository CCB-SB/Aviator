# Generated by Django 2.2.16 on 2020-12-07 15:42

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0004_auto_20201207_1639'),
    ]

    operations = [
        migrations.AlterField(
            model_name='curatedwebsite',
            name='authors',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.TextField(), blank=True, default=list, size=None),
        ),
    ]
