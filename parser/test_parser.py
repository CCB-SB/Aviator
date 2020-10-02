import os
import sys
import json
import csv
import gzip

#arguments
if len(sys.argv) < 3:
    print('usage: script.py <inputfile> <outputfile>')
    print('or: script.py <inputfile> <outputfile> <info_identifier> <html_identifier> <logs_identifier>')
    sys.exit(2)

inputfolder = sys.argv[1]
outputfile = sys.argv[2]
info_identifier = '.info.json' if len(sys.argv) < 4 else sys.argv[3]
html_identifier = '.html' if len(sys.argv) < 5 else sys.argv[4]
logs_identifier = '.logs.json.gz' if len(sys.argv) < 6 else sys.argv[5]

print("input folder: "+inputfolder)
print("output file: "+outputfile)
print("info file identifier: "+info_identifier)
print("html file identifier: "+html_identifier)
print("logs file identifier: "+logs_identifier)

#variables
folder_header = 'folder'
html_title_header = 'html-title'
analytics_header = 'analytics'
error_header = 'error'
result = {}
columns = [folder_header, html_title_header, analytics_header, error_header]
analytics_names = ['google-analytics', 'matomo', 'woopra', 'gosquared', 'go-squared', 'foxmetrics', 'fox-metrics', 'mixpanel', 'heap', 'statcounter', 'stat-counter', 'chartbeat', 'clicky', 'leadfeeder']

error_phrases = ['Not Found', 'Server unavailable', 'Maintainance in progress', 'Out server is down temporarily', 'Out server is down', '503 Service Unavailable', 
'Service Unavailable', 'No server is available to handle this request', 'has been discontinued', '502 Bad Gateway', 'The resource could not be found', '404 Not Found', 
'decimalLongitude', 'This is the default web page for this server', 'The web server software is running but no content has been added, yet', 
'The requested URL was not found on this server', 'The AlloPred web server is currently encountering some technical difficulties', 'Sorry, Page Not Found', 
'403 Forbidden', 'Error 406 - Not Acceptable', '<Error><Code>NoSuchBucket</Code>', 'Server under maintenance', '503 Service Temporarily Unavailable',
'Temporarily Unavailable', '502 Proxy Error', 'Proxy Erro', 'currently unavailable', 'invalid request', 'Service unavailable', 'The service is down due to technical issues', 
'Site Maintenance', 'Service temporarily unavailable', 'Object not found!', '404 Page Not Found', 'The page you requested was not found', 'Error 404', 
'ERROR 404: Seite nicht gefunden', 'Errore 404', 'ERROR 404', 'No longer supported webserver', 'This transfer is blocked', 
'The transfer has triggered a Web Application Firewall', 'The resource cannot be found', '403 - Forbidden: Access is denied',
'IP address could not be found', 'DNS_PROBE_FINISHED_NXDOMAIN', 'Your connection was interrupted', 'A network change was detected', 'ERR_NETWORK_CHANGED',
'This site can’t provide a secure connection', 'sent an invalid response', 'ERR_SSL_PROTOCOL_ERROR',
'This site can’t be reached', 'refused to connect', 'ERR_CONNECTION_REFUSED']

url_replacements = {
  "∼": "~"
}

def handleInfo(filename, target):
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
            if name == 'code' and value is not 200 and str(value).isnumeric():
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

for subdir, dirs, files in os.walk(inputfolder):
    for file in files:
        folder_name = subdir[subdir.rfind("\\") + 1:]
        if len(folder_name) > 0 and not folder_name in result.keys():
            result[folder_name] = {}
        if info_identifier in file:
            handleInfo(os.path.join(subdir, file), result[folder_name])
        elif html_identifier in file:
            handleHTML(os.path.join(subdir, file), result[folder_name])
        elif logs_identifier in file:
            handleLogs(os.path.join(subdir, file), result[folder_name])
    for key, value in result.items():
        result[key][folder_header] = key
        # replacements
        if "url" in result[key]:
            for replace in url_replacements:
                result[key]["url"] = result[key]["url"].replace(replace, url_replacements[replace])

try:
    with open(outputfile, 'w', newline='', encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns, delimiter=';')
        writer.writeheader()
        for key, value in result.items():
            writer.writerow(value)
except IOError:
    print("I/O error")