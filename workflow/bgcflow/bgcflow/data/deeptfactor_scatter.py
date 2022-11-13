import pandas as pd
from pathlib import Path
import json
import sys

def deeptfactor_to_json(prediction, outfile, keep_false=True):
    genome_id = Path(prediction).parent.stem
    df = pd.read_csv(prediction, sep="\t")
    df.loc[:, "genome_id"] = genome_id
    # rename columns
    df = df.rename(columns={"prediction" : "deeptfactor_prediction", "score" : "deeptfactor_score"})
    df = df.set_index("sequence_ID", drop=True)
    if not keep_false:
        df = df[df.deeptfactor_prediction == True]
    df.T.to_json(outfile)
    return

if __name__ == "__main__":
    deeptfactor_to_json(sys.argv[1], sys.argv[2])
