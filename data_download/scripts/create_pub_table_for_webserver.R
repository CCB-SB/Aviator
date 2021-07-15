library(data.table)

pub_tbl = fread(snakemake@input$tbl, colClasses=c(PMID="character"))
pub_tbl = pub_tbl[, c("PMID", "abstract", "title", "year", "journal", "authors", "URL")]

mails = fread(snakemake@input$emails, colClasses=c(PMID="character"))
mails = mails[!duplicated(mails)]
mails = mails[!duplicated(PMID)]
pubkw = fread(snakemake@input$mesh_terms, colClasses=c(PMID="character"))
pubkw = pubkw[, c("PMID", "keywords_all", "mesh_terms_all")]
pubkw = pubkw[!duplicated(pubkw)]

pub_tbl = merge(pub_tbl, mails, all.x = TRUE)

pub_tbl = merge(pub_tbl, pubkw)
pub_tbl[, Email:=mails$Email[match(PMID, mails$PMID)]]

# if year contains -, e.g. 2011-2012, take first year
pub_tbl[grep("-", year), year:=gsub("-.*", "", year)]

fwrite(pub_tbl, snakemake@output[[1]], sep='\t')
