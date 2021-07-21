#!/bin/bash

export PATH="$PATH":/usr/local/bin

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"
export SINGULARITY_BINDPATH="/local,/home"
export TIMESTAMP=$(date -Is)
cat cmds | parallel -j 13 > /tmp/dl.${TIMESTAMP}.log 2>&1
rsync -a /local/downloads/aviator/downloads/${TIMESTAMP} ccbadmin@ccb-compute2.cs.uni-saarland.de:/data/aviator_telemetry/
/home/tobias/miniconda3/bin/python update_dataset.py --token cRbtB795TX-Xa6Q8_yXf
/home/tobias/miniconda3/bin/python trigger_db_update.py --token cRbtB795TX-Xa6Q8_yXf
ssh ccbadmin@ccb-compute2.cs.uni-saarland.de rm -rf /data/aviator_telemetry/${TIMESTAMP}
rsync -a /local/downloads/aviator/downloads/${TIMESTAMP} main:git_projects/aviator_dl/downloads/ && rm -rf /local/downloads/aviator/downloads/${TIMESTAMP}

