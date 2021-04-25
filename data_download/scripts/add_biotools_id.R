library(data.table)
library(rjson)

biotools = rjson::fromJSON(file=snakemake@input$json)

tbl = rbindlist(lapply(c("Web service", "Web API", "Web application"), function(s) {
  pmids = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$pmid)), collapse = ';')))
  dois = unlist(lapply(biotools[[s]], function(x) paste(unlist(lapply(x$publication, function(y) y$doi)), collapse = ';')))
  homepage = unlist(lapply(biotools[[s]], function(x) x$homepage))
  biotoolsID = unlist(lapply(biotools[[s]], function(x) x$biotoolsID))
  data.table(pmid=pmids, doi=dois, homepage=homepage, biotoolsID=biotoolsID)
}))

filtered_tbl = fread(snakemake@input$filtered)

filtered_tbl[homepage %in% tbl$homepage, biotoolsID:=tbl$biotoolsID[match(homepage, tbl$homepage)]]

pmid2biotoolsID = filtered_tbl[, c("pmid", "biotoolsID")][, list(pmid=unlist(strsplit(pmid, ';'))), by="biotoolsID"]

fwrite(unique(pmid2biotoolsID), snakemake@output[[1]], sep='\t')
