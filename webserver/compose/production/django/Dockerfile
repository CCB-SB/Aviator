FROM registry.ccb-gitlab.cs.uni-saarland.de/webserver/miniconda3/debian/stretch:4.7.10

ENV LANG C.UTF-8

ENV PYTHONUNBUFFERED 1

RUN groupadd -r django \
    && useradd -r -g django django

# Requirements have to be pulled and installed here, otherwise caching won't work
COPY ./requirements /requirements

RUN pip install conda-merge==0.1.4 &&\
	conda-merge /requirements/base.yml /requirements/production.yml > /requirements/merged.yml &&\
	pip uninstall --yes conda-merge

RUN conda env update -n root -f /requirements/merged.yml &&\
    conda clean --all --yes &&\
    rm -rf /requirements

COPY --chown=django ./compose/production/django/entrypoint /entrypoint
RUN sed -i 's/\r//' /entrypoint
RUN chmod +x /entrypoint

COPY --chown=django ./compose/production/django/start /start
RUN sed -i 's/\r//' /start
RUN chmod +x /start

COPY --chown=django ./compose/production/django/create_superuser /create_superuser
RUN sed -i 's/\r//' /create_superuser
RUN chmod +x /create_superuser

COPY --chown=django ./compose/production/django/celery/worker/start /start-celeryworker
RUN sed -i 's/\r//' /start-celeryworker
RUN chmod +x /start-celeryworker

COPY --chown=django ./compose/production/django/celery/beat/start /start-celerybeat
RUN sed -i 's/\r//' /start-celerybeat
RUN chmod +x /start-celerybeat

COPY --chown=django ./compose/production/django/celery/flower/start /start-flower
RUN sed -i 's/\r//' /start-flower
RUN chmod +x /start-flower


RUN mkdir -p /app && chown django /app

COPY --chown=django . /app

RUN mkdir /app/staticfiles && chown -R django /app/staticfiles
RUN mkdir /app/aviator/media && chown -R django /app/aviator/media

USER django

RUN cd /app &&\
    CELERY_BROKER_URL=dummy\
    REDIS_URL=dummy DATABASE_URL=dummy DJANGO_SECRET_KEY=dummy DJANGO_ADMIN_URL=dummy\
    python /app/manage.py compress --settings config.settings.production &&\
    CELERY_BROKER_URL=dummy\
    REDIS_URL=dummy DATABASE_URL=dummy DJANGO_SECRET_KEY=dummy DJANGO_ADMIN_URL=dummy\
    python /app/manage.py collectstatic --settings config.settings.production --noinput


WORKDIR /app

ENTRYPOINT ["/entrypoint"]
