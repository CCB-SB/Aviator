# Generated by Django 2.2.16 on 2020-12-10 00:03

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0008_curatedwebsite_dates'),
    ]

    operations = [
        migrations.AddField(
            model_name='curatedwebsite',
            name='contact_mail',
            field=models.EmailField(blank=True, max_length=254, null=True),
        ),
    ]
