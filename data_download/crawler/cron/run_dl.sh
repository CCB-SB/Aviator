#!/bin/bash

export PATH="$PATH":/usr/local/bin

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"
export SINGULARITY_BINDPATH="/local"
export TIMESTAMP=$(date -Is)
cat cmds | parallel -j 13 > /tmp/dl.${TIMESTAMP}.log 2>&1
rsync -a /local/downloads/aviator/downloads/${TIMESTAMP} ccb-compute2.cs.uni-saarland.de:/data/aviator_telemetry/
python3 trigger_db_update.py --token cRbtB795TX-Xa6Q8_yXf
ssh ccb-compute2.cs.uni-saarland.de rm -rf /data/aviator_telemetry/${TIMESTAMP}
rsync -a /local/downloads/aviator/downloads/${TIMESTAMP} main:/share/data/Aviator/crawler_data/ && rm -rf /local/downloads/aviator/downloads/${TIMESTAMP}

