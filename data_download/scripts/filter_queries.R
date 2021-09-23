library(data.table)
library(stringi)
library(pbapply)

# for debugging
#save.image()
#stop()

tbl = fread(snakemake@input$pubmed, colClasses=c(PMID="character"))
biotools = fread(snakemake@input$biotools, colClasses = c(PMID="character"))
#biotools = biotools[year>=2010]
nar = fread(snakemake@input$nar, colClasses = c(PMID="character"))
nar_manual = fread(snakemake@input$nar_manual_cur, colClasses = c(PMID="character"))
other_manual = fread(snakemake@input$manual_cur, colClasses = c(PMID="character"))
other_pminfo = fread(snakemake@input$manual_cur_pubmed, colClasses = c(PMID="character"))

manual_urls = rbind(nar_manual, other_manual)[URL != ""]
setnames(manual_urls, "URL", "URL_Manual")

tbl = rbind(tbl, biotools, nar, other_pminfo, fill=TRUE)
tbl[, biotools_url:=biotools$biotools_url[match(PMID, biotools$PMID)]]

# merge duplicate pmids
tbl = tbl[, list(URL_BRACKET=paste(URL_BRACKET, collapse='; '), URL_XML=paste(URL_XML, collapse='; '), URL_Extractor=paste(URL_Extractor, collapse='; '), biotools_url=paste(URL_Extractor, collapse='; ')), by=c("PMID", "authors", "year", "title", "abstract", "journal")]

tbl = merge(tbl, manual_urls, all.x = TRUE)

# remove URLs for publications where we have manually curated entries
tbl[!is.na(URL_Manual), `:=`(URL_BRACKET="", URL_Extractor="", URL_XML="", biotools_url="")]
tbl[is.na(URL_Manual), URL_Manual:=""]

# exclude PMIDs that are not ws
no_ws_pmids = readLines(snakemake@input$no_ws_list)
tbl = tbl[!PMID %in% no_ws_pmids]

fixes = c("www.AntibioticScout. ch"="www.AntibioticScout.ch",
          "www.donnainformata-mammografia"="www.donnainformata-mammografia.it",
          "http://pgdbj. jp/en/"="http://pgdbj.jp/en/",
          "http://myeloid-risk.\nedu/"="http://myeloid-risk.edu/",
          "https://bioserver.iiita.ac.in/DeEPn/.Communicated"="https://bioserver.iiita.ac.in/DeEPn/",
          "NeuroMorpho.Org"="neuromorpho.org",
          " https://www.cdc.gov/vaccines/hcp/acip-recs/vacc-specific/flu.html.These"="https://www.cdc.gov/vaccines/hcp/acip-recs/vacc-specific/flu.html",
          "DermPatientEd.Com"="DermPatientEd.com",
          "www.envihomolog.eawag.ch.Graphical"="www.envihomolog.eawag.ch",
          "http://bioinfo.lifl.fr/norine/smiles2monomers.jsp.Graphical"="http://bioinfo.lifl.fr/norine/smiles2monomers.jsp",
          "AKUTNE.CZ"="akutne.cz",
          "http://bl210.caspur.it/MODEL-DB/MODEL-DB_web/MODindex.php.Operating"="http://bl210.caspur.it/MODEL-DB/MODEL-DB_web/MODindex.php",
          "http://www.ensembl.org.The"="http://www.ensembl.org",
          " http://www.ensemblgenomes.org.In"=" http://www.ensemblgenomes.org",
          "GlycanStructure.Org"="glycanstructure.org",
          "http://www.\u2028neals.org"="http://www.neals.org",
          "http://sit.mfu.ac.th/lcgdb/index_FoxM1.php.Communicated"="http://sit.mfu.ac.th/lcgdb/index_FoxM1.php")

for(e in names(fixes)){
  pos = grep(e, tbl$abstract)
  set(tbl, pos, "URL_BRACKET", fixes[e])
  set(tbl, pos, "URL_XML", fixes[e])
  set(tbl, pos, "URL_Extractor", fixes[e])
}

# remove clinical trials, removes many false positives, but also a tiny fraction of real web servers
tbl = tbl[!grepl("clinical trial", abstract, ignore.case = T) | (PMID %in% manual_urls$PMID)]

# filter URLs according to regex
url_regex = "^(https?://(?:www\\.|(?!www))[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|www\\.[a-zA-Z0-9][a-zA-Z0-9-]+[a-zA-Z0-9]\\.[^\\s]{2,}|https?://(?:www\\.|(?!www))[a-zA-Z0-9]+\\.[^\\s]{2,}|www\\.[a-zA-Z0-9]+\\.[^\\s]{2,})"

tbl[, URL:=unlist(pblapply(strsplit(paste(URL_BRACKET, URL_XML, URL_Extractor, biotools_url, URL_Manual, sep = ';'), ';', fixed=T), function(x) {
  res = unique(stringi::stri_trim_both(gsub("^\\s+", "", trimws(x, whitespace="[ \t\r\n.)(,;>\\[\\]]"))))
  
  # replace tilde
  res = gsub("∼", "~", res)
  
  # replace non-braking space with space
  res = gsub(intToUtf8(160), " ", res)
  # remove starting //
  res = gsub("^//", "", res)
  # remove number-number
  res = gsub("^\\d+\\.\\d+-\\d+\\.\\d+$", "", res)
  # remove alphanum-dot-alphanum+
  res = gsub("^\\w\\.\\w+$", "", res)
  # remove number minus
  res = gsub("^\\d+-.*", "", res)
  # remove starting without alnum
  res = gsub("^\\W.*", "", res)
  # remove number dot number space
  res = gsub("^\\d+\\.\\d+\\s+", "", res)
  # remove all
  # one dot starting with less than 3 alphanum 
  # one dot ending with numbers only
  # one dot starting with digits
  # without dot
  to_discard = (stringr::str_count(res, stringr::fixed(".")) == 0) | (stringr::str_count(res, stringr::fixed(".")) == 1 &
                                                                        (grepl("^\\w\\.", res) | grepl("^\\w\\w\\.", res) | grepl("\\d+$", res) | grepl("^\\d+", res)))
  # ends with dot uppercase alnum
  to_discard = to_discard | grepl("\\.[A-Z]\\w+$", res)
  # whitespace following dot
  to_discard = to_discard | grepl("\\.\\s", res)
  # email - @
  to_discard = to_discard | grepl("@", res)
  
  # non-ascii character urls
  to_discard = to_discard | unlist(lapply(strsplit(res, ""), function(x) any(unlist(lapply(x, utf8ToInt)) > 255)))
  
  # ending with '
  to_discard = to_discard | grepl("'$", res)
  
  # containing ""
  to_discard = to_discard | grepl('""', res)
  
  # containing </ 
  to_discard = to_discard | grepl("</", res)
  
  # containing \\
  to_discard = to_discard | grepl("\\", res, fixed=TRUE)
  
  # ending with .pdf, .zip, .gz
  to_discard = to_discard | grepl("\\.pdf|\\.zip|\\.gz", res)
  
  # ending with ):
  to_discard = to_discard | grepl("):", res, fixed=TRUE)
  
  # ending with }
  to_discard = to_discard | grepl("\\}$|\\{$", res)
  
  # urls to exclude
  to_exclude = c("youtube.com", "github.com", "bitbucket.org", "gitlab.com", "hub.docker.com", "code.google.com",
                 "bitbucket.com", "code.google", "cran.r-project", "sourceforge",
                 "pypi.org", "raw.githubusercontent.com", "bioconductor.org", "www.crd.york.ac.uk",
                 "ClinicalTrials.gov", "ClinicalTrial.gov", "readthedocs",
                 "isrctn.com", "monkeysurvey.com", "anzctr.org.au", "trialregister.nl",
                 "ncbi.nlm.nih.gov/genome/", "irct.ir", "clinicaltrialsregister.eu",
                 "Preprints.org", "clinicaltrials.gouv", "euroqol.org", "randomization.com",
                 "Trials.gov", "randomize.net", "RR2-10", "amazon.co", "yelp.com",
                 "apps.who.int", "controlled-trials.com", "webcitation.org", "linkedin.com", "facebook.com",
                 "google", "trial", "guidetopharmacology", "onlinelibrary", "doi.org",
                 "www.ncbi.nlm.nih.gov", "www.ncbi.nim.nih.gov", "http://www.w3.org", '""""',
                 "wiki.nci.nih.gov", "journals.sagepub", "www.mdpi.com",
                 "goo.gl", "youtu.be", "twitter", "survey", "trialregister",
                 " or", "-the$", "www.cdc.gov", "supplementary",
                 "europepmc", "grdr.ncats.nih.gov", "professional.diabetes", "www.ema.europa.eu",
                 "journals.na.lww.com", "orthoguidelines", "guideline", "mammalsociety", "researchregistry", "commondataelements.ninds.nih.gov",
                 "www.interpol.int", "academic.oup.com", "oxfordjournals.org", "sciencedirect", "rhinologyjournal",
                 "scirp.org", "health.gov.on.ca", "dentistryresearch.org", "tripod-statement.org",
                 "paineurope.com", "www.who.int", "www.fda.gov", "umin.ac.jp", "random.org",
                 "www.drks.de", "evidenceupdate-tatarstan.ru", "genefriends.org",
                 "springer.com", "eucast.org",
                 "prepareforyourcare.org", "ivf-worldwide.com", "www.nasn.org",
                 "iwantthekit.org", "menopausematters.co", "cure4kids", "schoology",
                 "www.acr.org", "nobelprize.org", "acousticalsociety.org",
                 "contraceptionchoices.org",
                 "bing.com; yahoo.com; duckduckgo",
                 "www.gov.uk", "webofknowledge.com", "ijpc.com",
                 "cartilage.org", "guide"
                 )
  to_discard = to_discard | grepl(paste(to_exclude, collapse = "|"), res, ignore.case = T)
  res[to_discard] = ""
  # remove alphanum=
  res = gsub("^\\w*=.*", "", res)
  res = res[res != ""]
  res = res[is.na(suppressWarnings(as.numeric(res)))]
  
  new_urls = unique(grep(url_regex, res, value=TRUE, perl=TRUE))
  new_urls = new_urls[!duplicated(fs::path_sanitize(new_urls, replacement="."))]

  paste(new_urls, collapse='; ')
}, cl = snakemake@threads))]

# ensure manually added URLs are kept
tbl[URL == "" & URL_Manual != "", URL:=URL_Manual]

tbl = tbl[URL != ""]
tbl = tbl[!duplicated(tbl)]

##### normalize URLs and write mapping original URL -> normalized URL
orig_urls = unlist(strsplit(tbl$URL, '; '))
new_urls = unlist(lapply(orig_urls, function(u) {
  new_u = u
  if(!startsWith(u, "http")){
    new_u = gsub("^", "http://", new_u)
  }
  new_u
}))

orig_new_urls = data.table(orig=orig_urls, new=new_urls, new_fs_norm=fs::path_sanitize(new_urls, replacement = "."),
                           final_new="")

for(i in 1:nrow(orig_new_urls)) {
  u = orig_new_urls$new_fs_norm[i]
  if(startsWith(u, "http..")){
    u_http = u
    u_https = gsub("http..", "https..", u)
  } else { # https
    u_http = gsub("https..", "http..", u)
    u_https = u
  }
  
  # prioritize https, then http
  # multiple hits -> take shortest url
  which_https = which(u_https == orig_new_urls$new_fs_norm)
  if(length(which_https) > 0){
      orig_new_urls$final_new[i] = orig_new_urls$new[which_https[which.min(nchar(orig_new_urls$new[which_https]))]]
      next
  }
  which_http = which(u_http == orig_new_urls$new_fs_norm)
  if(length(which_http) > 0){
      orig_new_urls$final_new[i] = orig_new_urls$new[which_http[which.min(nchar(orig_new_urls$new[which_http]))]]
  }
}
orig_new_df = copy(orig_new_urls[, c("orig", "final_new")])
setnames(orig_new_df, "final_new", "new")
orig_new_df[, new_ID:=fs::path_sanitize(new, replacement = ".")]
orig_new_df[, orig_ID:=fs::path_sanitize(orig, replacement = ".")]

fwrite(unique(orig_new_df), snakemake@output$orig_url2new_id, sep='\t')

# update table urls
tbl[, URL:=unlist(lapply(URL, function(x) {
  urls = strsplit(x, '; ')[[1]]
  new_urls = orig_new_df$new[match(urls, orig_new_df$orig)]
  new_urls = unique(new_urls)
  paste(new_urls, collapse='; ')
}))]

fwrite(tbl[order(as.numeric(PMID))], snakemake@output$filtered_tbl, sep='\t')

urls = unique(unlist(strsplit(tbl$URL, '; ')))
ids = fs::path_sanitize(urls, replacement = ".")
urls2query = data.table(URL=urls, ID=ids)
fwrite(urls2query, snakemake@output$urls2query, sep='\t')

# split in X batches
nbatches = snakemake@params$batches
chunks = split(urls2query, f=as.factor(rep(seq(1, nbatches), len=nrow(urls2query))))

for(i in 1:length(chunks)){
  fwrite(chunks[[i]], snakemake@output$urls2query_batches[i], sep='\t')
}
