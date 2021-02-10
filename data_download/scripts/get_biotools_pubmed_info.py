from tqdm import tqdm
import time
import pandas as pd

from Bio import Entrez

from collections import defaultdict

Entrez.email = "tobias.fehlmann@gmail.com"
Entrez.api_key = "2a26cfc6c7e92ad9eee8a2b29615e8a14209"

tbl = pd.read_csv(snakemake.input.biotools, sep='\t')

pmid_list = list(set(e for p in tbl["pmid"] for e in p.split(';')))
pmid2url = defaultdict(list)

for p in pmid_list:
    urls = list(tbl[tbl["pmid"].str.contains(p)]["homepage"])
    pmid2url[p].extend(urls)

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

step_size = 1000

def get_author_name(author):
    if "LastName" in author and "Initials" in author:
        return f"{author['LastName']} {author['Initials']}"
    if "LastName" in author:
        return author['LastName']
    if "CollectiveName" in author:
        return author["CollectiveName"]
    raise RuntimeError("Could not determine author name {}".format(author))

pmid_info = dict()
for i in tqdm(range(0, len(pmid_list), step_size)):
    result = fetch_details(pmid_list[i:min(i+step_size, len(pmid_list))])
    for e in result["PubmedArticle"]:
        pmid = str(e["MedlineCitation"]["PMID"])
        article = e["MedlineCitation"]["Article"]
        abstract = article.get("Abstract", {"AbstractText": [""]})["AbstractText"]
        abstract_text = "\n".join(str(e) for e in abstract)
        title = article["ArticleTitle"]
        pubdate = article["Journal"]["JournalIssue"]["PubDate"]
        year = pubdate.get("Year", pubdate.get("MedlineDate", "")).split(' ')[0]
        journal_info = article["Journal"]
        if "AuthorList" in article:
            authors = ", ".join(get_author_name(a) for a in article["AuthorList"])
        else:
            authors = ""
        pmid_info[pmid] = {
            "abstract": abstract_text,
            "title": title,
            "year": year,
            "journal": journal_info.get("ISOAbbreviation", journal_info["Title"]),
            "authors": authors
        }


pmid_tbl = pd.DataFrame.from_dict(pmid_info, orient='index')
pmid_tbl.index.name = "PMID"
pmid_tbl["biotools_url"] = [";".join(pmid2url[p]) for p in pmid_tbl.index]
pmid_tbl.to_csv(snakemake.output.csv, sep='\t')
