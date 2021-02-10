import gitlab
from gitlab.exceptions import GitlabGetError

import sys
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument("--token", help="Token for gitlab account that will create users and repositories.", required=True)

args = parser.parse_args()

gl = gitlab.Gitlab("https://ccb-gitlab.cs.uni-saarland.de", private_token=args.token)

project = gl.projects.get("pascal/aviator")
jobs = project.jobs.list(page=1, per_page=100)

for j in jobs:
    if j.name == "update_websitecalls_and_stats-production":
        break
else:
    print("No update job found...")
    sys.exit(1)

try:
    j.play()
except GitlabGetError as e:
    print("Could not play job... retrying...")
    j.retry()

# update in case we got a new id
time.sleep(10)
jobs = project.jobs.list(page=1, per_page=100)

for j in jobs:
    if j.name == "update_websitecalls_and_stats-production":
        break

# wait until job is done
j.refresh()
while j.status == "running":
    time.sleep(10)
    j.refresh()

print("Updating ended with status: {}".format(j.status))
