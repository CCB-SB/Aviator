import requests
from urllib.parse import quote 
from collections import defaultdict
from time import sleep
from tqdm import tqdm
import pickle
import pandas as pd
import json

terms = ["Web service", "Web API", "Web application"]
optional_terms = ["Bioinformatics portal", "Database portal"]

base_url = 'https://bio.tools/api/t/?format=json&toolType="{}"'

result = defaultdict(list)
 
for t in tqdm(terms + optional_terms):
    res = requests.get(base_url.format(t))
    total_items = res.json()["count"]
    print("Found {} entries for '{}'".format(total_items, t))
    if res.json()["next"] is not None:
        url = base_url.format(t)+res.json()["next"].replace("?", "&")
    else:
        url = None
    result[t].extend(res.json()["list"])
    with tqdm(total=total_items) as pbar:
        pbar.update(len(res.json()["list"]))
        while url is not None:
            sleep(0.1)
            res = requests.get(url)
            if res.json()["next"] is not None:
                url = base_url.format(t)+res.json()["next"].replace("?", "&")
            else:
                url = None
            result[t].extend(res.json()["list"])
            pbar.update(len(res.json()["list"]))

# code for "all" tools
#base_url = 'https://bio.tools/api/t/?format=json'
# res = requests.get(base_url)
# total_items = res.json()["count"]
# print("Found {} entries".format(total_items))
# if res.json()["next"] is not None:
#     url = base_url+res.json()["next"].replace("?", "&")
# else:
#     url = None
#
#result["all"].extend(res.json()["list"])
# with tqdm(total=total_items) as pbar:
#     pbar.update(len(res.json()["list"]))
#     while url is not None:
#         sleep(0.1)
#         res = requests.get(url)
#         if res.json()["next"] is not None:
#             url = base_url+res.json()["next"].replace("?", "&")
#         else:
#             url = None
#         result["all"].extend(res.json()["list"])
#         pbar.update(len(res.json()["list"]))


with open(snakemake.output.pkl, 'wb') as f:
    pickle.dump(result, f)

with open(snakemake.output.json, 'w') as f:
    json.dump(result, f)
