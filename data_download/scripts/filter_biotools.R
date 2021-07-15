library(data.table)
library(curl)
library(rjson)
library(yaml)
library(reshape2)
library(pbapply)

pboptions(type="timer")

biotools = rjson::fromJSON(file=snakemake@input$biotools)

tbl = rbindlist(lapply(c("Web service", "Web API", "Web application", "Bioinformatics portal", "Database portal"), function(s) {
  pmids = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$pmid)), collapse = ';')))
  dois = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$doi)), collapse = ';')))
  homepage = unlist(lapply(biotools[[s]], function(x) x$homepage))
  data.table(pmid=pmids, doi=dois, homepage=homepage)
}))

tbl = tbl[pmid != "" | doi != ""]
tbl = tbl[!duplicated(tbl)]

doi2pmid_df = fread(snakemake@input$doi2pmid)

tbl[pmid == "" & doi != "", pmid:=unlist(lapply(doi, function(x) {
  dois = unlist(strsplit(x, ';')[[1]])
  pmids = doi2pmid_df[match(dois, doi)]$pmid
  pmids = pmids[!is.na(pmids)]
  paste(pmids, collapse=';')
}))]

conv_result = pblapply(unique(unlist(strsplit(tbl[pmid == "" & doi != ""]$doi, ';'))), function(x) lapply(x, function(e) system(sprintf('NCBI_API_KEY="2a26cfc6c7e92ad9eee8a2b29615e8a14209" /home/tobias/miniconda3/envs/metapub/bin/convert doi2pmid \'%s\'; sleep 1', e), intern=TRUE)), cl=5)
names(conv_result) = unique(unlist(strsplit(tbl[pmid == "" & doi != ""]$doi, ';')))

# doi2pmid = unlist(lapply(names(conv_result), function(doi){
#   r = conv_result[[doi]][[1]]
#   if(length(r) == 1){
#     return("")
#   } else{
#     pmid = suppressWarnings(as.numeric(stringr::str_extract(r[2], "\\d+")))
#     ifelse(is.na(pmid), "", as.character(pmid))
#   }
# }))
# 
# doi2pmid_df = data.table(doi=names(conv_result),
#                          pmid=doi2pmid)
# fwrite(doi2pmid_df[pmid != ""], "data/doi2pmid.csv")

rescued_pmids = unlist(lapply(tbl[pmid == "" & doi != ""]$doi, function(x) {
  ids = unlist(lapply(strsplit(x, ';')[[1]], function(e){
    r = conv_result[[e]][[1]]
    if(length(r) == 1){
      return("NA")
    } else{
      stringr::str_extract(r[2], "\\d+")
    }
  }))
  ids = suppressWarnings(as.numeric(ids))
  ids = ids[!is.na(ids)]
  paste(ids, collapse=';')
  })
)

tbl[pmid == "" & doi != "", pmid:=rescued_pmids]

tbl = tbl[!duplicated(tbl) & pmid != ""]
fwrite(tbl, snakemake@output[[1]], sep='\t')
