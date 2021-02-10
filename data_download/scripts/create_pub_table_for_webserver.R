library(data.table)

pub_tbl = fread(snakemake@input$tbl)
pub_tbl = pub_tbl[, c("PMID", "abstract", "title", "year", "journal", "authors", "URL")]

mails = fread(snakemake@input$emails)
mails = mails[!duplicated(mails)]
pubkw = fread(snakemake@input$mesh_terms)
pubkw = pubkw[, c("PMID", "keywords_all", "mesh_terms_all")]
pubkw = pubkw[!duplicated(pubkw)]

pub_tbl = merge(pub_tbl, mails, all.x = TRUE)

pub_tbl = merge(pub_tbl, pubkw)
pub_tbl[, Email:=mails$Email[match(PMID, mails$PMID)]]

fwrite(pub_tbl, snakemake@output[[1]], sep='\t')
