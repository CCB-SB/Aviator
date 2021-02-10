FROM selenium/node-chrome:85.0.4183.83-chromedriver-85.0.4183.83

USER root

RUN groupadd -g 1000 -r crawler \
    && useradd -u 1000 -r -g crawler crawler

RUN mkdir /app && chown crawler /app

RUN apt-get update && apt-get install -y python3-pip

COPY requirements.txt /tmp
RUN pip3 install -r /tmp/requirements.txt

COPY download_website.py /usr/local/bin

USER crawler

ENTRYPOINT [ "python3", "/usr/local/bin/download_website.py" ]
