from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

import environ

def create_superuser():
    env = environ.Env()
    username = env("ADMIN_USER")
    password = env("ADMIN_PASS")
    email = env("ADMIN_EMAIL")

    if User.objects.filter(username=username).count() == 0:
        User.objects.create_superuser(username, email, password)
        print('Superuser {} created.'.format(username))
    else:
        print('Superuser creation skipped. {} exists already.'.format(username))


class Command(BaseCommand):
    def handle(self, *args, **options):
        create_superuser()
