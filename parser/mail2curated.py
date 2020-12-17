import pandas as pd

from Bio import Entrez

import argparse
import os.path

from glob import glob



parser = argparse.ArgumentParser()
parser.add_argument("--email-dir", required=True)
parser.add_argument("--output", required=True)

args = parser.parse_args(["--email-dir=/home/tobias/git_projects/aviator_pascal/webserver/data/curated_mails", "--output=/home/tobias/git_projects/aviator_pascal/webserver/data/curated.csv"])

Entrez.email = "tobias.fehlmann@gmail.com"
Entrez.api_key = "2a26cfc6c7e92ad9eee8a2b29615e8a14209"


mails = glob(os.path.join(args.email_dir, "*"))

def parse_mail(mail):
    result = dict()
    kws = ["Name", "Email", "Tool", "PubMed ID", "URL",
           "API", "Authors", "Abstract", "Tags",
           "I want to be informed, if the website is down",
           "Days before reminder",
           "Comments"
    ]
    with open(mail) as f:
        current_kw = ""
        for l in f:
            l = l.replace("+", " ")
            for kw in kws:
                if l.startswith(f"{kw}:"):
                    kwline = True
                    current_kw = kw
                    break
            else:
                kwline = False
            
            if kwline:
                result[current_kw] = l[len(current_kw)+1:].strip()
            else:
                result[current_kw] += '\n' + l.strip()
    return result


def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def get_author_name(author):
    if "LastName" in author and "Initials" in author:
        return f"{author['LastName']} {author['Initials']}"
    if "LastName" in author:
        return author['LastName']
    if "CollectiveName" in author:
        return author["CollectiveName"]
    raise RuntimeError("Could not determine author name {}".format(author))


def get_pmid_info(pmid_list):
    pmid_info = dict()
    step_size = 1000
    for i in range(0, len(pmid_list), step_size):
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

    return pd.DataFrame.from_dict(pmid_info).transpose()


mail_info = pd.DataFrame(parse_mail(m) for m in mails)

pmid_info = get_pmid_info(mail_info["PubMed ID"])

mail_info["year"] = list(pmid_info.loc[mail_info["PubMed ID"], "year"])
mail_info["journal"] = list(pmid_info.loc[mail_info["PubMed ID"], "journal"])

mail_info.to_csv(args.output, sep='\t', index=False)
