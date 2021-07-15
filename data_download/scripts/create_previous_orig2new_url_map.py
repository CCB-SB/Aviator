import pandas as pd
import json

df = pd.concat([pd.read_csv(f, sep='\t') for f in snakemake.input])

urlmap = df.groupby("new").agg({"orig": lambda els: list(set(els))}).to_dict()["orig"]

with open(snakemake.output[0], 'w', encoding='ascii') as f:
    json.dump(urlmap, f)

