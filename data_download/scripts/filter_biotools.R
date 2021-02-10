library(data.table)
library(curl)
library(rjson)
library(yaml)
library(reshape2)
library(pbapply)

biotools = rjson::fromJSON(file=snakemake@input$biotools)

tbl = rbindlist(lapply(c("Web service", "Web API", "Web application"), function(s) {
  pmids = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$pmid)), collapse = ';')))
  dois = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$doi)), collapse = ';')))
  homepage = unlist(lapply(biotools[[s]], function(x) x$homepage))
  data.table(pmid=pmids, doi=dois, homepage=homepage)
}))

tbl = tbl[pmid != "" | doi != ""]
tbl = tbl[!duplicated(tbl)]

conv_result = pblapply(unique(unlist(strsplit(tbl[pmid == "" & doi != ""]$doi, ';'))), function(x) lapply(x, function(e) system(sprintf('NCBI_API_KEY="2a26cfc6c7e92ad9eee8a2b29615e8a14209" /home/tobias/miniconda3/envs/metapub/bin/convert doi2pmid \'%s\'', e), intern=TRUE)), cl=3)
names(conv_result) = unique(unlist(strsplit(tbl[pmid == "" & doi != ""]$doi, ';')))

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
#save.image("conv_result.RData")

#View(tbl[pmid != ""])
