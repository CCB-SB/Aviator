# Generated by Django 2.2.16 on 2020-12-07 15:39

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0003_auto_20201207_1605'),
    ]

    operations = [
        migrations.RenameField(
            model_name='curatedwebsite',
            old_name='tool_name',
            new_name='title',
        ),
    ]
