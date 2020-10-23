# Generated by Django 2.2.16 on 2020-10-23 06:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='publication',
            name='journal',
            field=models.TextField(db_index=True),
        ),
        migrations.AlterField(
            model_name='publication',
            name='year',
            field=models.SmallIntegerField(db_index=True),
        ),
    ]
