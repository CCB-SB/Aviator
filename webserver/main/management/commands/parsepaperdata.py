from django.core.management.base import BaseCommand, CommandError
from main.models import Website
from main.models import Paper
import os, csv, json, gzip

class Command(BaseCommand):
    help = 'Creates the paper model from a given csv'

    def add_arguments(self, parser):
        parser.add_argument('csv_file')
        parser.add_argument('url_filter')

    def handle(self, *args, **options):
        try:
            csv_file = options['csv_file']
        except:
            raise CommandError('Please provide an input folder')
        header_line = True
        num = 0
        num2 = 0
        num_updated = 0
        header = dict()

        #handle filter
        filter = False
        filter_ids = []
        filter_original = []
        filter_derived = []
        try:
            filename = options['url_filter']
            with open(filename, 'r', encoding='utf-8') as csvfile:
                csv_reader = csv.reader(csvfile, delimiter='\t')
                filter = True
                for row in csv_reader:
                    counter = 0
                    id_url = ""
                    original_url = ""
                    for entry in row:
                        if counter == 0:
                            id_url = str(entry).encode('unicode-escape').decode('utf-8')
                        elif counter == 1:
                            original_url = str(entry).encode('unicode-escape').decode('utf-8')
                        elif counter == 2:
                            derived_url = str(entry).encode('unicode-escape').decode('utf-8')
                        elif counter > 2:
                            break
                        counter += 1
                    filter_ids.append(id_url)
                    filter_original.append(original_url)
                    filter_derived.append(derived_url)
            self.stdout.write("Filter applied")
        except:
            self.stdout.write("No Filter applied")

        with open(csv_file, 'r', encoding='utf-8') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter='\t')
            for row in csv_reader:
                if header_line:
                    header_line = False
                    counter = 0
                    for entry in row:
                        header[counter] = entry
                        counter += 1
                    continue
                counter = 0
                ok = True
                urls = []
                title = ""
                pubmed_id = ""
                abstract = ""
                authors = ""
                year = ""
                journal = ""
                paper = Paper()
                for entry in row:
                    if counter in header:
                        if header[counter] == 'title':
                            title = str(entry).encode('unicode-escape').decode('utf-8')
                            paper.title = title
                        if header[counter] == 'PMID':
                            pubmed_id = str(entry).encode('unicode-escape').decode('utf-8')
                            paper.pubmed_id = pubmed_id
                        if header[counter] == 'abstract':
                            abstract = str(entry).encode('unicode-escape').decode('utf-8')
                            paper.abstract = abstract
                        if header[counter] == 'authors':
                            authors = str(entry).encode('unicode-escape').decode('utf-8')
                            paper.authors = authors
                        if header[counter] == 'year':
                            year = entry[0:4]
                            paper.year = year
                        if header[counter] == 'journal':
                            journal = str(entry).encode('unicode-escape').decode('utf-8')
                            paper.journal = journal
                        #if header[counter] == 'comment' and entry == 'delete':
                        #    ok = False
                        #    break
                        if header[counter] == 'URL':
                            urls = str(entry).encode('unicode-escape').decode('utf-8').split(';')
                            if(filter):
                                for url in urls:
                                    if url not in filter_ids and url not in filter_original and url not in filter_derived:
                                        urls.remove(url)
                                if len(urls) == 0:
                                    ok = False
                                    break
                            paper.url = url
                    counter += 1
                if ok:
                    new_url = ""
                    for url in urls:
                        new_url += (url + ", ")
                    new_url[:-2]
                    papers = Paper.objects.filter(url=new_url, title=title)
                    if papers.count() > 0:
                        for old_paper in papers:
                            old_paper.pubmed_id = pubmed_id
                            old_paper.abstract = abstract
                            old_paper.authors = authors
                            old_paper.year = year
                            old_paper.journal = journal
                            old_paper.save()
                            num_updated += 1
                            continue
                    num += 1
                    paper.save()
                    #websites = Website.objects.filter(url=url)
                    websites = Website.objects.filter(url=new_url)
                    for website in websites:
                        website.papers.add(paper)
                        website.save()
                        num2 += 1
            self.stdout.write(self.style.SUCCESS('Successfully added ' + str(num) + ' new Papers and updated ' + str(num_updated) + ' Papers with ' + str(num2) + ' connections to webpages'))
