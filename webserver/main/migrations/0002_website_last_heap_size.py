# Generated by Django 2.2.16 on 2020-12-06 17:59

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='website',
            name='last_heap_size',
            field=models.IntegerField(default=0),
        ),
    ]
