# Generated by Django 2.2.16 on 2021-04-24 14:45

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0016_publication_biotools_id'),
    ]

    operations = [
        migrations.AlterField(
            model_name='publication',
            name='biotools_id',
            field=models.CharField(default='', max_length=50),
        ),
    ]
