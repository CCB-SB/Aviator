from django.core.management.base import BaseCommand, CommandError
from main.models import Website
from main.models import WebsiteCall
from main.models import Paper
import os, csv, json, gzip
from datetime import datetime
from tqdm import tqdm
import pytz

class Command(BaseCommand):
    help = 'Creates the website / website calls model and a csv from a given source folder'

    def add_arguments(self, parser):
        parser.add_argument('folder')
        #parser.add_argument('csv_target')

    def handle(self, *args, **options):
        # variables
        info_identifier = '.info.json'
        html_identifier = '.html'
        logs_identifier = '.logs.json.gz'
        csv_identifier = '.logs.json.gz'

        folder_header = 'folder'
        html_title_header = 'html-title'
        analytics_header = 'analytics'
        error_header = 'error'
        datetime_header = 'datetime'
        result = {}
        original_urls = {}
        derived_urls = {}
        columns = [folder_header, html_title_header, analytics_header, error_header, datetime_header]
        analytics_names = ['google-analytics', 'matomo', 'woopra', 'gosquared', 'go-squared', 'foxmetrics',
                           'fox-metrics', 'mixpanel', 'heap', 'statcounter', 'stat-counter', 'chartbeat', 'clicky',
                           'leadfeeder']

        error_phrases = ['Not Found', 'Server unavailable', 'Maintainance in progress',
                         'Out server is down temporarily', 'Out server is down', '503 Service Unavailable',
                         'Service Unavailable', 'No server is available to handle this request',
                         'has been discontinued', '502 Bad Gateway', 'The resource could not be found', '404 Not Found',
                         'decimalLongitude', 'This is the default web page for this server',
                         'The web server software is running but no content has been added, yet',
                         'The requested URL was not found on this server',
                         'The AlloPred web server is currently encountering some technical difficulties',
                         'Sorry, Page Not Found',
                         '403 Forbidden', 'Error 406 - Not Acceptable', '<Error><Code>NoSuchBucket</Code>',
                         'Server under maintenance', '503 Service Temporarily Unavailable',
                         'Temporarily Unavailable', '502 Proxy Error', 'Proxy Erro', 'currently unavailable',
                         'invalid request', 'Service unavailable', 'The service is down due to technical issues',
                         'Site Maintenance', 'Service temporarily unavailable', 'Object not found!',
                         '404 Page Not Found', 'The page you requested was not found', 'Error 404',
                         'ERROR 404: Seite nicht gefunden', 'Errore 404', 'ERROR 404', 'No longer supported webserver',
                         'This transfer is blocked',
                         'The transfer has triggered a Web Application Firewall', 'The resource cannot be found',
                         '403 - Forbidden: Access is denied',
                         'IP address could not be found', 'DNS_PROBE_FINISHED_NXDOMAIN',
                         'Your connection was interrupted', 'A network change was detected', 'ERR_NETWORK_CHANGED',
                         'This site can’t provide a secure connection', 'sent an invalid response',
                         'ERR_SSL_PROTOCOL_ERROR',
                         'This site can’t be reached', 'refused to connect', 'ERR_CONNECTION_REFUSED']

        url_replacements = {
            "∼": "~"
        }

        def winapi_path(dos_path, encoding=None):
            if (not isinstance(dos_path, str) and encoding is not None):
                dos_path = dos_path.decode(encoding)
            path = os.path.abspath(dos_path)
            if path.startswith(u"\\\\"):
                return u"\\\\?\\UNC\\" + path[2:]
            return u"\\\\?\\" + path

        def handleInfo(filename, target):
            target[datetime_header] = filename[filename.rfind("\\") + 1:][0:19].replace("_", ":") #e.g. 2020-09-06T16_01_22 => 2020-09-06T16:01:22
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'r') as file:
                    data = json.load(file)
                    iterateInfoJSON(data, target, '')
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    data = json.load(file)
                    iterateInfoJSON(data, target, '')

        def iterateInfoJSON(source, target, name_prefix):
            for key, value in source.items():
                if isinstance(value, dict):
                    iterateInfoJSON(value, target, name_prefix + key + '_')
                else:
                    name = name_prefix + str(key)
                    target[name] = str(value).replace("\r\n", "").replace("\n", "")
                    if name == 'code' and value != 200 and str(value).isnumeric():
                        if error_header in target:
                            if len(target[error_header]) > 0:
                                target[error_header] = target[error_header] + ', '
                            target[error_header] = target[error_header] + str(value)
                    if (not name in columns):
                        columns.append(name)

        def handleHTML(filename, target):
            text = ""
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'r') as file:
                    text = file.read().replace('\n', '')
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    text = file.read().replace('\n', '')

            index = text.find("</title>")
            if (index > 0):
                substr = text[:index]
                title = substr[substr.rfind(">") + 1:]
                target[html_title_header] = title
            target[error_header] = ''
            for error_phrase in error_phrases:
                if error_phrase in text and not error_phrase in target[error_header]:
                    if len(target[error_header]) > 0:
                        target[error_header] = target[error_header] + ', '
                    target[error_header] = target[error_header] + error_phrase

        def handleLogs(filename, target):
            target[analytics_header] = ''
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'r') as file:
                    data = json.load(file)
                    iterateLogJSON(data, target)
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    data = json.load(file)
                    iterateLogJSON(data, target)

        def iterateLogJSON(source, target):
            if isinstance(source, dict):
                for key, value in source.items():
                    if isinstance(value, dict) or isinstance(value, list):
                        iterateLogJSON(value, target)
                    elif isinstance(value, str):
                        for name in analytics_names:
                            if name in value and not name in target[analytics_header]:
                                if len(target[analytics_header]) > 0:
                                    target[analytics_header] = target[analytics_header] + ', '
                                target[analytics_header] = target[analytics_header] + name
            elif isinstance(source, list):
                for value in source:
                    iterateLogJSON(value, target)

        def handleCSV(filename, target1, target2):
            with open(filename, 'r', encoding='utf-8') as csvfile:
                csv_reader = csv.reader(csvfile, delimiter='\t')
                for row in csv_reader:
                    counter = 0
                    id_url = ""
                    original_url = ""
                    derived_url = ""
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
                    target1[id_url] = original_url
                    target2[id_url] = derived_url

        #Iterate New data
        try:
            inputfolder = options['folder']
        except:
            raise CommandError('Please provide an input folder')
        self.stdout.write("Starting to iterate and parse data:")
        for subdir, dirs, files in tqdm(os.walk(inputfolder)):
            for file in files:
                folder_name = subdir[subdir.rfind("\\") + 1:]
                if len(folder_name) > 0 and not folder_name in result.keys():
                    result[folder_name] = {}
                if info_identifier in file:
                    handleInfo((os.path.join(subdir, file)), result[folder_name])
                    #handleInfo(winapi_path(os.path.join(subdir, file)), result[folder_name])
                elif html_identifier in file:
                    handleHTML((os.path.join(subdir, file)), result[folder_name])
                    #handleHTML((os.path.join(subdir, file)), result[folder_name])
                elif logs_identifier in file:
                    handleLogs((os.path.join(subdir, file)), result[folder_name])
                    #handleLogs(winapi_path(os.path.join(subdir, file)), result[folder_name])
                elif csv_identifier in file:
                    handleCSV((os.path.join(subdir, file)), original_urls, derived_urls)
                    #handleCSV(winapi_path(os.path.join(subdir, file)), original_urls, derived_urls)
            for key, value in result.items():
                result[key][folder_header] = key
                # replacements
                if "url" in result[key]:
                    for replace in url_replacements:
                        result[key]["url"] = result[key]["url"].replace(replace, url_replacements[replace])

        new_websites = 0
        new_websitecalls = 0
        new_paper_connections = 0
        create_websites = []
        update_websites = []
        #Create new Websites
        self.stdout.write("Creating new websites:")
        for key, value in tqdm(result.items()):
            if "url" not in result[key]:
                continue
            wp_url = result[key]["url"]
            derived_url = derived_urls[wp_url] if wp_url in derived_urls else ""
            wp_url = original_urls[wp_url] if wp_url in original_urls else wp_url
            #Check and create Webpage if there is none
            website = None
            status = None
            if "code" in result[key] and result[key]["code"] == "200":
                status = True
            elif "code" in result[key]:
                status = False
            if "ok" in result[key] and result[key]["ok"] == "Pass":
                status = True
            if "ok" in result[key] and result[key]["ok"] != "Pass":
                status = False
            if "error" in result[key] and len(result[key]["error"]) > 0:
                status = False
            websites = Website.objects.filter(url=wp_url)
            if websites.count() == 0:
                website = Website(url=wp_url)
                website.ip = result[key]["ip"] if "ip" in result[key] else 0
                website.server = result[key]["response_headers_server"] if "response_headers_server" in result[key] else ""
                website.analytics = result[key]["analytics"] if "analytics" in result[key] else ""
                website.timezone = result[key]["response_headers_Pragma"] if "response_headers_Pragma" in result[key] else ""
                website.certificate_secure = result[key]["certificateSecurityState"] == "secure" if "certificateSecurityState" in result[key] else False
                website.script = ""
                if "response_headers_Set-Cookie" in result[key]:
                    if "PHPSESSID" in result[key]["response_headers_Set-Cookie"]:
                        website.script += "PHP"
                    if "JSESSIONID" in result[key]["response_headers_Set-Cookie"]:
                        website.script += "Java"
                website.derived_url = derived_url
                create_websites.append(website)
                new_websites += 1

                papers = Paper.objects.filter(url=wp_url)
                for paper in papers:
                    website.papers.add(paper)
                    new_paper_connections += 1
                website.status = status
            else:
                website = Website.objects.filter(url=wp_url)[0]
                website.status = status
                update_websites.append(website)

        Website.objects.bulk_update(update_websites, ["status"])
        Website.objects.bulk_create(create_websites)
        #Create new WebsiteCalls
        self.stdout.write("Creating new website calls:")
        create_websitecalls = []
        for key, value in tqdm(result.items()):
            if "url" not in result[key]:
                continue
            wp_url = result[key]["url"]
            wp_url = original_urls[wp_url] if wp_url in original_urls else wp_url
            website = Website.objects.filter(url=wp_url)[0]
            #Create Website Call
            call = WebsiteCall(website=website)
            if datetime_header in result[key]:
                dt_string = result[key][datetime_header] #e.g. 2020-09-06T16:01:22
                utc = pytz.timezone('UTC')
                dt = datetime.strptime(dt_string, "%Y-%m-%dT%H:%M:%S")
                dt = utc.localize(dt, is_dst=True)
                call.datetime = dt
            else:
                call.datetime = datetime.now(tzinfo=pytz.UTC)
            call.ok = result[key]["ok"] == "Pass" if "ok" in result[key] else False
            call.error = result[key]["error"] if "error" in result[key] else ""
            call.msg = result[key]["msg"] if "msg" in result[key] else ""
            call.code = result[key]["code"] if "code" in result[key] and result[key]["code"] != "NA" else 0
            call.json_data = json.dumps(result[key])
            create_websitecalls.append(call)
            new_websitecalls += 1
            #Set status of website
            website.status = status
            website.save()
        WebsiteCall.objects.bulk_create(create_websitecalls)
        #Save as csv
        #try:
        #    outputfile = options['csv_target']
        #    with open(outputfile, 'w', newline='', encoding="utf-8") as csvfile:
        #        writer = csv.DictWriter(csvfile, fieldnames=columns, delimiter=';')
        #        writer.writeheader()
        #        for key, value in result.items():
        #            writer.writerow(value)
        #except:
        #    self.stdout.write("CSV was not saved, no target file given or an IO Error occured")
        self.stdout.write(self.style.SUCCESS('Successfully added ' + str(new_websitecalls) + ' website calls with ' + str(new_websites) + ' new websites with ' + str(new_paper_connections) + ' new connections to papers '))
