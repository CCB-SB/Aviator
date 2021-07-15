## Create/update the list of web servers

1. ```bash
   snakemake copy_old_ws_data
   snakemake --use-conda
   snakemake copy_to_app
   ```
   
2. Zip results 
   
   ```bash
   zip -r results.$(date -I).zip results
   ```

3. Push changes
   
4. On blackbox: pull repository. Update files in crawler directory if necessary.

5. To pull the docker crawler image via singularity:

   ```bash
   SINGULARITY_DOCKER_USERNAME=YOUR_USERNAME SINGULARITY_DOCKER_PASSWORD=YOUR_PASSWORD singularity pull docker://registry.ccb-gitlab.cs.uni-saarland.de/ccb-staff/aviator/crawler
   ```

6. Update gitlab-ci variables (update the **timestamp**):

```bash
   BIOTOOLSID_PATH=data/biotoolsid2pmid.2021-07-15.csv
   FOLDERID_URL_PATH=data/urls2query.2021-07-15.csv
   PUBLICATIONS_PATH=data/publication_table_for_webservice.2021-07-15.csv
   NEW2ORIG_PREV_PATH=data/new2orig_urls_all_previous.2021-07-15.json
   ORIG2NEW_PATH=data/orig2new_urls_and_id.2021-07-15.csv
```

7. To update the database run **before** cron job website call the **updatedataset** job via the gitlab-ci

## To manually add new publication and URL
Update data/curated_urls.csv
