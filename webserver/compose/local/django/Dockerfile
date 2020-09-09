FROM registry.ccb-gitlab.cs.uni-saarland.de/webserver/miniconda3/debian/stretch:4.7.10

ENV LANG C.UTF-8

ENV PYTHONUNBUFFERED 1


# install pkill for celery autoreload
RUN apt-get update && apt-get install -y procps


COPY ./requirements /requirements

RUN pip install conda-merge==0.1.4 &&\
	conda-merge /requirements/base.yml /requirements/local.yml > /requirements/merged.yml &&\
	pip uninstall --yes conda-merge

RUN conda env update -n root -f /requirements/merged.yml &&\
    conda clean --all --yes

COPY ./compose/production/django/entrypoint /entrypoint
RUN sed -i 's/\r//' /entrypoint
RUN chmod +x /entrypoint

COPY ./compose/local/django/start /start
RUN sed -i 's/\r//' /start
RUN chmod +x /start

COPY ./compose/production/django/create_superuser /create_superuser
RUN sed -i 's/\r//' /create_superuser
RUN chmod +x /create_superuser

COPY ./compose/local/django/celery/worker/start /start-celeryworker
RUN sed -i 's/\r//' /start-celeryworker
RUN chmod +x /start-celeryworker

COPY ./compose/local/django/celery/beat/start /start-celerybeat
RUN sed -i 's/\r//' /start-celerybeat
RUN chmod +x /start-celerybeat

COPY ./compose/local/django/celery/flower/start /start-flower
RUN sed -i 's/\r//' /start-flower
RUN chmod +x /start-flower


WORKDIR /app

ENTRYPOINT ["/entrypoint"]
