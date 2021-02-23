from django.core.management.base import BaseCommand, CommandError
from main.models import Website
from main.models import WebsiteCall
from main.models import WebsiteStatus
from main.models import Publication
import os, gzip
import orjson as json
from datetime import datetime, timedelta
from tqdm import tqdm
import pytz
import pandas as pd
from glob import glob
from os.path import join, basename
from collections import defaultdict
from django.contrib.postgres.aggregates import BoolOr
import ahocorasick




class Command(BaseCommand):
    help = 'Creates the website / website calls model and a csv from a given source folder'

    def add_arguments(self, parser):
        parser.add_argument('folder')
        parser.add_argument('origurl2folder')
        parser.add_argument('url_filter')
        # parser.add_argument('csv_target')

    def handle(self, *args, **options):
        # variables
        info_identifier = '.info.json'
        html_identifier = '.html'
        logs_identifier = '.logs.json.gz'
        csv_identifier = '.csv'

        folder_header = 'folder'
        html_title_header = 'html-title'
        analytics_header = 'analytics'
        error_header = 'error'
        datetime_header = 'datetime'

        columns = [folder_header, html_title_header, analytics_header, error_header,
                   datetime_header]
        analytics_names = ['google-analytics', 'matomo', 'woopra', 'gosquared', 'go-squared',
                           'foxmetrics',
                           'fox-metrics', 'mixpanel', 'heap', 'statcounter', 'stat-counter',
                           'chartbeat', 'clicky',
                           'leadfeeder', 'piwik']

        def make_automaton(words, idx_only=False):
            result = ahocorasick.Automaton()
            for i, w in enumerate(words):
                if idx_only:
                    result.add_word(w, i)
                else:
                    result.add_word(w, w)
            result.make_automaton()
            return result

        analytics_automaton = make_automaton(analytics_names)

        error_phrases = ['Not Found', 'Server unavailable', 'Maintainance in progress',
                         'Out server is down temporarily', 'Out server is down',
                         '503 Service Unavailable',
                         'Service Unavailable', 'No server is available to handle this request',
                         'has been discontinued', '502 Bad Gateway',
                         'The resource could not be found', '404 Not Found',
                         'decimalLongitude', 'This is the default web page for this server',
                         'The web server software is running but no content has been added, yet',
                         'The requested URL was not found on this server',
                         'The AlloPred web server is currently encountering some technical difficulties',
                         'Sorry, Page Not Found',
                         '403 Forbidden', 'Error 406 - Not Acceptable',
                         '<Error><Code>NoSuchBucket</Code>',
                         'Server under maintenance', '503 Service Temporarily Unavailable',
                         'Temporarily Unavailable', '502 Proxy Error', 'Proxy Erro',
                         'currently unavailable',
                         'invalid request', 'Service unavailable',
                         'The service is down due to technical issues',
                         'Site Maintenance', 'Service temporarily unavailable', 'Object not found!',
                         '404 Page Not Found', 'The page you requested was not found', 'Error 404',
                         'ERROR 404: Seite nicht gefunden', 'Errore 404', 'ERROR 404',
                         'No longer supported webserver',
                         'This transfer is blocked',
                         'The transfer has triggered a Web Application Firewall',
                         'The resource cannot be found',
                         '403 - Forbidden: Access is denied',
                         'IP address could not be found', 'DNS_PROBE_FINISHED_NXDOMAIN',
                         'Your connection was interrupted', 'A network change was detected',
                         'ERR_NETWORK_CHANGED',
                         'This site can’t provide a secure connection', 'sent an invalid response',
                         'ERR_SSL_PROTOCOL_ERROR',
                         'This site can’t be reached', 'refused to connect',
                         'ERR_CONNECTION_REFUSED']

        error_phrases_automaton = make_automaton(error_phrases)
        url_replacements = {
            "∼": "~"
        }

        def handleInfo(filename, target):
            target[datetime_header] = os.path.basename(filename)[0:19].replace("_",
                                                                               ":")  # e.g. 2020-09-06T16_01_22 => 2020-09-06T16:01:22
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'r') as file:
                    data = json.loads(file.read())
                    target.update(data)
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    data = json.loads(file.read())
                    target.update(data)
                    code_val = target.get("code", None)
                    if code_val is not None and code_val != 200 and str(code_val).isnumeric():
                        if error_header in target:
                            if len(target[error_header]) > 0:
                                target[error_header] = f"{target[error_header]}, {code_val}"
                            else:
                                target[error_header] = f"{code_val}"

        def handleHTML(filename, target):
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'r') as file:
                    text = file.read().replace('\n', '')
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    text = file.read().replace('\n', '')

            index = text.find("</title>")
            if index > 0:
                substr = text[:index]
                title = substr[substr.rfind(">") + 1:]
                target[html_title_header] = title

            shiny_kws = ["shiny-singletons", "shiny-server-client", "shiny.min.js"]
            for k in shiny_kws:
                if text.find(k) > 0:
                    target["programming"] = "Shiny"
                    break

            dash_kws = ["dash-component-suites", "_dash-renderer"]
            for k in dash_kws:
                if text.find(k) > 0:
                    target["programming"] = "Dash"
                    break

            error_phrase_hits = set()
            for end_index, error_phrase in error_phrases_automaton.iter(text):
                error_phrase_hits.add(error_phrase)
            target[error_header] = ", ".join(sorted(error_phrase_hits))

        def handleLogs(filename, target):
            target[analytics_header] = ''
            if filename[-3:] == '.gz':
                with gzip.open(filename, 'rt') as file:
                    analytics_hits = set()
                    # for l in file:
                    for _, a in analytics_automaton.iter(file.read()):
                        analytics_hits.add(a)
                    target[analytics_header] = ", ".join(sorted(analytics_hits))
                    # data = json.loads(file.read())
                    # iterateLogJSON(data, target)
            else:
                with open(filename, 'r', encoding="utf-8") as file:
                    analytics_hits = set()
                    for l in file:
                        for _, a in analytics_automaton.iter(l):
                            analytics_hits.add(a)
                    target[analytics_header] = ", ".join(sorted(analytics_hits))
                    # data = json.loads(file.read())
                    # iterateLogJSON(data, target)

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

        filter_orig_urls = set()
        if options["url_filter"] is not None  and \
            options["url_filter"] not in {"None", "null", ""}:
            filter_tbl = pd.read_csv(options['url_filter'], sep='\t')
            filter_orig_urls = set(filter_tbl["Original URL"])
            self.stdout.write("Filter applied")
        else:
            self.stdout.write("No Filter applied")

        # Iterate New data
        try:
            inputfolder = options['folder']
        except:
            raise CommandError('Please provide an input folder')
        self.stdout.write("Starting to iterate and parse data:")

        origurl2folder = pd.read_csv(options["origurl2folder"], sep='\t')
        folder2orig_urls = origurl2folder.groupby("ID")["URL"].apply(set).to_dict()

        # either we get multiple telemetry folders or just one
        telemetry_folders = glob(join(inputfolder, "20*-*-*T*"))
        if len(telemetry_folders) == 0:
            telemetry_folders = [inputfolder]
        self.stdout.write(f"Found {len(telemetry_folders)} telemetry folders")

        new_websites = 0
        new_websitecalls = 0
        new_paper_connections = 0
        for tf in tqdm(telemetry_folders):
            result = {}
            all_tbls = pd.concat((pd.read_csv(f, sep='\t') for f in glob(join(tf, "*.csv"))))
            all_tbls.set_index("Original URL", inplace=True)
            orig_urls = set(all_tbls.index)
            all_dict = all_tbls.to_dict()
            # update websites
            websites_to_update = list(Website.objects.filter(original_url__in=orig_urls))
            websites_to_update_final = []
            for w in websites_to_update:
                new_url = all_dict["Derived URL"][w.original_url]
                if w.derived_url != all_dict["Derived URL"][w.original_url]:
                    w.derived_url = new_url
                websites_to_update_final.append(w)

            Website.objects.bulk_update(websites_to_update_final, fields=["derived_url"])

            for subdir, dirs, files in tqdm(os.walk(tf)):
                for file in files:
                    folder_name = subdir[subdir.rfind("\\") + 1:]
                    if len(folder_name) > 0 and folder_name not in result.keys():
                        result[folder_name] = {}
                    if info_identifier in file:
                        handleInfo((os.path.join(subdir, file)), result[folder_name])
                    elif html_identifier in file:
                        handleHTML((os.path.join(subdir, file)), result[folder_name])
                    elif logs_identifier in file:
                        handleLogs((os.path.join(subdir, file)), result[folder_name])
                for key, value in result.items():
                    result[key][folder_header] = key
                    # replacements
                    if "url" in result[key]:
                        for replace in url_replacements:
                            result[key]["url"] = result[key]["url"].replace(replace,
                                                                            url_replacements[
                                                                                replace])

            create_websites = []
            update_websites = []
            # Create new Websites
            self.stdout.write("Creating new websites:")
            already_done = set()

            orig_url_2_website = {w.original_url: w for w in Website.objects.all()}
            for key, value in tqdm(result.items()):
                if "url" not in value:
                    continue
                # get original url(s) for this folder
                orig_urls = folder2orig_urls.get(basename(value[folder_header]), None)
                if orig_urls is None:
                    self.stderr.write(
                        f"No folder 2 original url mapping found for folder {basename(value[folder_header])}")
                    continue

                if filter_orig_urls:
                    orig_urls = set(u for u in orig_urls if u in filter_orig_urls)
                if not orig_urls:
                    continue

                # Check and create Webpage if there is none
                status = WebsiteStatus.OFFLINE
                if "code" in value and value["code"] == 200:
                    status = WebsiteStatus.ONLINE
                elif "code" in value:
                    status = WebsiteStatus.OFFLINE
                if "ok" in value and value["ok"] != "Pass":
                    status = WebsiteStatus.OFFLINE
                if "error" in value and len(value["error"]) > 0:
                    status = WebsiteStatus.OFFLINE

                for o in orig_urls:
                    if o not in already_done:
                        already_done.add(o)
                        ip = value.get("ip", None)
                        if ip == "NA":
                            ip = None
                        if o not in orig_url_2_website:
                            website = Website(original_url=o,
                                              status=status,
                                              ip=ip,
                                              derived_url=value["url"]
                                              )
                        else:
                            website = orig_url_2_website[o]
                        # update if we can reach the website, or if it's the first time we query the website
                        if status == WebsiteStatus.ONLINE or o not in orig_url_2_website:
                            website.ip = ip
                            website.server = value.get("response_headers_server", "")
                            website.analytics = value.get("analytics", "")
                            website.timezone = value.get("response_headers_Pragma", "")
                            website.certificate_secure = value.get("certificateSecurityState",
                                                                   "unknown") == "secure"
                            website.script = value.get("programming", "")
                            if "response" in value and value.get("response") is not None:
                                if "headers" in value.get("response"):
                                    if "server" in value.get("response").get("headers"):
                                        website.server = value.get("response").get("headers").get("server", "")
                                    if "Pragma" in value.get("response").get("headers"):
                                        website.timezone = value.get("response").get("headers").get("Pragma", "")
                                    if "Set-Cookie" in value.get("response").get("headers"):
                                        c2lang = {"ASPSESSIONID": "ASP", "JSESSIONID": "Java", "PHPSESSID": "PHP", "CGISESSID": "Perl"}
                                        for k, v in c2lang.items():
                                            if k in value.get("response").get("headers").get("Set-Cookie"):
                                                website.script = "{}/{}".format(website.script, v)
                        if o not in orig_url_2_website:
                            new_websites += 1
                            create_websites.append(website)
                        else:
                            update_websites.append(website)


            create_websites = Website.objects.bulk_create(create_websites)
            Website.objects.bulk_update(update_websites, ["server", "analytics", "timezone", "certificate_secure", "script", "ip"])

            # connect new websites to publications
            origurl2paper = defaultdict(list)
            for p in Publication.objects.all():
                for u in p.url:
                    origurl2paper[u].append(p)

            # Create connections to paper
            for website in create_websites:
                pubs = origurl2paper.get(website.original_url, None)
                if pubs is None:
                    self.stderr.write(
                        f"No publications could be found for URL {website.original_url}!")
                else:
                    website.papers.set(pubs)
                    new_paper_connections += len(pubs)

            # update mapping to latest
            orig_url_2_website = {w.original_url: w for w in Website.objects.all()}

            # Create new WebsiteCalls
            self.stdout.write("Creating new website calls:")
            create_websitecalls = []
            for key, value in tqdm(result.items()):
                if "url" not in value:
                    continue

                # get original url(s) for this folder
                orig_urls = folder2orig_urls.get(basename(value[folder_header]), None)
                if orig_urls is None:
                    self.stderr.write(
                        f"No folder 2 original url mapping found for folder {basename(value[folder_header])}")
                    continue

                if filter_orig_urls:
                    orig_urls = set(u for u in orig_urls if u in filter_orig_urls)
                if not orig_urls:
                    continue

                websites = {orig_url_2_website[u] for u in orig_urls}

                # Create Website Call(s)
                for w in websites:
                    call = WebsiteCall(website=w)
                    if datetime_header in value:
                        dt_string = value[datetime_header]  # e.g. 2020-09-06T16:01:22
                        utc = pytz.timezone('UTC')
                        dt = datetime.strptime(dt_string, "%Y-%m-%dT%H:%M:%S")
                        dt = utc.localize(dt, is_dst=True)
                        call.datetime = dt
                    else:
                        call.datetime = datetime.now(tz=pytz.UTC)
                    call.ok = value.get("ok", "unk") == "Pass"
                    call.error = value.get("error", "")
                    call.msg = value.get("msg", "")
                    call.code = value["code"] if "code" in value and not pd.isnull(
                        value["code"]) and not value["code"] == "NA" else 0
                    call.ok = call.ok and call.error == "" and call.code == 200
                    call.json_data = value
                    create_websitecalls.append(call)
                    new_websitecalls += 1

            WebsiteCall.objects.bulk_create(create_websitecalls)



        # Save as csv
        # try:
        #    outputfile = options['csv_target']
        #    with open(outputfile, 'w', newline='', encoding="utf-8") as csvfile:
        #        writer = csv.DictWriter(csvfile, fieldnames=columns, delimiter=';')
        #        writer.writeheader()
        #        for key, value in result.items():
        #            writer.writerow(value)
        # except:
        #    self.stdout.write("CSV was not saved, no target file given or an IO Error occured")
        self.stdout.write(self.style.SUCCESS(
            f'Successfully added {new_websitecalls} website calls with {new_websites} new websites with {new_paper_connections} new connections to papers '))
