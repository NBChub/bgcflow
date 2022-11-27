import json
import logging
import sys

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def summarize_gtdb_json(accession_list, df_gtdb_output):
    """
    Given a string of json filepaths generated by gtdb_prep.py,
    merge them together into one csv table
    """

    # Reading input
    logging.info("Reading GTDB metadata .json files...")
    accession = accession_list.split()
    out = []
    for a in accession:
        with open(a, "r") as f:
            out.append(json.load(f))
    df = pd.DataFrame(out).set_index("genome_id", drop=False)

    # Getting taxonomic information
    logging.info("Getting taxonomic information...")
    df_taxonomy = pd.DataFrame.from_dict(
        {i: df.loc[i, "gtdb_taxonomy"] for i in df.index}
    ).T
    # Split between species (specific epithet) and organism
    df_taxonomy["organism"] = df_taxonomy["species"]
    for idx in df_taxonomy.index:
        try:
            df_taxonomy.loc[idx, "species"] = df_taxonomy.loc[idx, "species"].split(
                " "
            )[1]
        except IndexError:  # leave blank for empty taxonomy
            pass
    # Tidying up
    df_taxonomy.columns = [c.title() for c in df_taxonomy.columns]

    # Getting other metadata
    try:
        logging.info("Getting metadata into table...")
        metadata = pd.DataFrame.from_dict(
            {i: df.loc[i, "metadata"] for i in df.index}
        ).T

        metadata_container = []
        for c in metadata.columns:
            value = {i: metadata.loc[i, c] for i in metadata.index}
            try:
                metadata_container.append(pd.DataFrame.from_dict(value).T)
            except ValueError:
                metadata_container.append(
                    pd.DataFrame([value]).T.rename(columns={0: c})
                )
        df_metadata = pd.concat(metadata_container, axis=1)

        logging.info("Finalizing table...")
        df_final = pd.concat(
            [df.loc[:, ["genome_id", "gtdb_release"]], df_taxonomy, df_metadata], axis=1
        )

    except KeyError:
        logging.info("No additional metadata found.")
        logging.info("Finalizing table...")
        df_final = pd.concat(
            [df.loc[:, ["genome_id", "gtdb_release"]], df_taxonomy], axis=1
        )

    # save to file
    logging.info(f"Writing to file: {df_gtdb_output}")
    df_final.to_csv(df_gtdb_output, index=False)
    return None


if __name__ == "__main__":
    summarize_gtdb_json(sys.argv[1], sys.argv[2])
