from tqdm import tqdm
import time
import pandas as pd

from Bio import Entrez

from collections import defaultdict
import sys


Entrez.email = "tobias.fehlmann@gmail.com"
Entrez.api_key = "2a26cfc6c7e92ad9eee8a2b29615e8a14209"

tbl = pd.read_csv(sys.argv[1], sep='\t')
tbl["PMID"] = tbl["PMID"].astype(str)
pmid_list = [e for e in tbl["PMID"]]

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

step_size = 500

pmid_info = dict()
for i in tqdm(range(0, len(pmid_list), step_size)):
    result = fetch_details(pmid_list[i:min(i+step_size, len(pmid_list))])
    for e in result["PubmedArticle"]:
        pmid = str(e["MedlineCitation"]["PMID"])
        meshlist = e["MedlineCitation"].get("MeshHeadingList", [])
        mesh_terms_all = [str(h["DescriptorName"]) for h in meshlist]
        mesh_terms_major_only = [str(h["DescriptorName"]) for h in meshlist if h["DescriptorName"].attributes["MajorTopicYN"] == "Y"]

        kwlist = e["MedlineCitation"].get("KeywordList", [])
        if len(kwlist) == 1:
            kwlist = kwlist[0]
        kw_all = [str(k) for k in kwlist]
        kw_major_only = [str(k) for k in kwlist if k.attributes.get("MajorTopicYN", "NA") == "Y"]

        if not meshlist and not kwlist:
            print(pmid)

        pmid_info[pmid] = {
            "mesh_terms_all": ";".join(mesh_terms_all),
            "mesh_terms_major_only": ";".join(mesh_terms_major_only),
            "keywords_all": ";".join(kw_all),
            "keywords_major_only": ";".join(kw_major_only)
        }


pmid_tbl = pd.DataFrame.from_dict(pmid_info, orient='index')
pmid_tbl.index.name = "PMID"
pmid_tbl.to_csv(sys.argv[2], sep='\t')
