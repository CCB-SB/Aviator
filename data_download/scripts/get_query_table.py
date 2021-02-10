from tqdm import tqdm
import time
import pandas as pd

from Bio import Entrez

query = """(webserver [Title/Abstract] OR "web server" [Title/Abstract] OR web-server [Title/Abstract] OR "web service" [Title/Abstract] OR web-service [Title/Abstract] OR "web portal" [Title/Abstract] OR web-portal [Title/Abstract] OR "online tool" [Title/Abstract] OR website [Title/Abstract] OR "web based" [Title/Abstract] OR web-based [Title/Abstract] OR "web interface" [Title/Abstract] OR web-interface [Title/Abstract] OR "web tool" [Title/Abstract] OR web-tool [Title/Abstract] OR "network service" [Title/Abstract] OR network-service [Title/Abstract] OR server[Title/Abstract]) AND (("2010/01/01"[Date - Publication] : "3000"[Date - Publication]))"""

Entrez.email = "tobias.fehlmann@gmail.com"
Entrez.api_key = "2a26cfc6c7e92ad9eee8a2b29615e8a14209"
handle = Entrez.esearch(db='pubmed',
                        retmax='1000000',
                        retmode='xml',
                        term=query
)
pmids = Entrez.read(handle)

pmid_list = pmids["IdList"]

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
pmid_tbl.to_csv(snakemake.output.csv, sep='\t')
