import pandas as pd

from tqdm import tqdm
from urlextract import URLExtract

extractor = URLExtract()
extractor.update_when_older(7)

tbl = pd.read_csv(snakemake.input.query, sep='\t')
java_urls = pd.read_csv(snakemake.input.java_urls, sep='\t')
tbl = tbl.merge(java_urls, on="PMID")

def reformat_abstract(e):
    val = e.replace(';', ' ').replace('\xa0', ' ')
    if "http:/" in e and not "http://" in e:
        val = val.replace('http:/', 'http://')
    if "https:/" in e and not "https://" in e:
        val = val.replace('https:/', 'https://')
    if val[-1] == '.':
        return val[:-1]
    else:
        return val

tbl["URL_Extractor"] = [";".join(set(extractor.find_urls("" if pd.isnull(e) else reformat_abstract(e)))) for e in tqdm(tbl["abstract"])]
tbl["has_url_kw"] = ["http://" in e or "https://" in e or "www." in e if not pd.isnull(e) else False for e in tqdm(tbl["abstract"])]

tbl.to_csv(snakemake.output.csv, sep='\t', index=False)
