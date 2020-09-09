import subprocess

from django.core.management.base import BaseCommand
from django.utils import autoreload


def restart_celery():
    subprocess.run(["pkill", "celery"])
    subprocess.run(["/start-celeryworker"])


class Command(BaseCommand):
    def handle(self, *args, **options):
        print('Starting celery worker with autoreload...')
        autoreload.main(restart_celery)
