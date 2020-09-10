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
google_analytics_header = 'google-analytics'
result = {}
columns = [folder_header, html_title_header, google_analytics_header]

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


def handleLogs(filename, target):
    target[google_analytics_header] = 'NA'
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
                if 'google-analytics' in value:
                    target['google-analytics'] = 'true'
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
try:
    with open(outputfile, 'w', newline='', encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns, delimiter=';')
        writer.writeheader()
        for key, value in result.items():
            writer.writerow(value)
except IOError:
    print("I/O error")