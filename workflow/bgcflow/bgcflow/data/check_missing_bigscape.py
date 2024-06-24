import ast
import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def correct_bigscape_results(
    input_dir,
    antismash_region_table,
    bigscape_log,
    time_stamp,
    output_dir,
    cutoffs=["0.3", "0.4", "0.5"],
):
    """
    Corrects BiG-SCAPE results by adding dropped BGCs back into the analysis.

    This function reads the BiG-SCAPE log to identify BGCs that were removed from the analysis
    due to having no domains found. It then creates new family assignments for these BGCs and
    updates the relevant BiG-SCAPE output tables accordingly.

    Parameters:
    - input_dir (Path): Directory containing the BiG-SCAPE output files.
    - antismash_region_table (str): Path to the CSV file containing antiSMASH region data.
    - bigscape_log (str): Path to the BiG-SCAPE log file.
    - time_stamp (str): Timestamp used in the naming of BiG-SCAPE output files.
    - output_dir (Path): Directory where corrected output files will be saved.
    - cutoffs (list of str): List of similarity cutoffs used in BiG-SCAPE analysis.

    Returns:
    None
    """
    logging.info("Checking BiG-SCAPE log for dropped BGCs...")
    with open(bigscape_log, "r") as f:
        bigscape_log_content = f.readlines()

    antismash_regions = pd.read_csv(antismash_region_table)
    removed_from_bigscape = [
        i.strip("'  No domains where found in ").strip(
            ".domtable. Removing it from further analysis\n"
        )
        for i in bigscape_log_content
        if "Removing it from further analysis" in i
    ]
    dropped_bgcs = antismash_regions[
        antismash_regions.bgc_id.isin(removed_from_bigscape)
    ]

    for cutoff in cutoffs:
        bigscape_families_table = input_dir / f"{time_stamp}_df_families_{cutoff}.csv"
        bigscape_family_presence_table = (
            input_dir / f"{time_stamp}_df_family_presence_{cutoff}.csv"
        )
        bigscape_cluster_table = input_dir / f"{time_stamp}_df_clusters_{cutoff}.csv"
        bigscape_network_table = input_dir / f"{time_stamp}_df_network_{cutoff}.csv"
        bigscape_known_table = input_dir / f"{time_stamp}_df_known_{cutoff}.csv"

        df_bigscape_families = pd.read_csv(bigscape_families_table)
        df_bigscape_family_presence = pd.read_csv(bigscape_family_presence_table)
        df_bigscape_cluster = pd.read_csv(bigscape_cluster_table)
        df_bigscape_network = pd.read_csv(bigscape_network_table, index_col=0)

        if len(dropped_bgcs) > 0:
            logging.debug(
                f"WARNING: Found {len(dropped_bgcs)} BGCs that are dropped by BiG-SCAPE"
            )
            logging.info("Correcting BiG-SCAPE results...")
            corrected_bgc_cluster_assignment = {}
            for num, i in enumerate(dropped_bgcs.index):
                dropped_bgc = dropped_bgcs.loc[i, "bgc_id"]
                genome_id = dropped_bgcs.loc[i, "genome_id"]
                logging.debug(
                    f"Adding {dropped_bgc} from {genome_id} to the result table with cutoff {cutoff}..."
                )
                last_fam_id = df_bigscape_families.iloc[-1, 0]
                new_fam_id = last_fam_id + 1
                df_bigscape_families.loc[len(df_bigscape_families.index)] = [
                    new_fam_id,
                    "unassigned",
                    f"u_unassigned_{new_fam_id}",
                    1,
                    np.nan,
                ]
                corrected_bgc_cluster_assignment[num] = {
                    "bgc_id": dropped_bgc,
                    "product": ".".join(
                        ast.literal_eval(dropped_bgcs.loc[i, "product"])
                    ),
                    "genome_id": genome_id,
                    "accn_id": dropped_bgcs.loc[i, "accession"],
                    f"fam_id_{cutoff}": new_fam_id,
                    f"fam_type_{cutoff}": "unassigned",
                    f"fam_known_compounds_{cutoff}": f"u_unassigned_{new_fam_id}",
                    "bigscape_class": "unassigned",
                }
                new_edge_id = df_bigscape_network.index[-1] + 1
                df_bigscape_network.loc[new_edge_id, "Clustername 1"] = dropped_bgc
                df_bigscape_network.loc[new_edge_id, "Clustername 2"] = dropped_bgc
                family_presence_index = df_bigscape_family_presence[
                    df_bigscape_family_presence.genome_id == genome_id
                ].index
                df_bigscape_family_presence.loc[family_presence_index, new_fam_id] = 1
                df_bigscape_family_presence[new_fam_id] = (
                    df_bigscape_family_presence[new_fam_id].fillna(0).astype(int)
                )

            df_bigscape_cluster_missing = pd.DataFrame.from_dict(
                corrected_bgc_cluster_assignment
            ).T
            df_bigscape_cluster = pd.concat(
                [df_bigscape_cluster, df_bigscape_cluster_missing]
            )
        else:
            logging.info(f"No BGCs were dropped by BiG-SCAPE in cutoff {cutoff}")
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        logging.debug(f"Writing updated tables for cutoff {cutoff} to {output_dir}...")
        df_bigscape_families.to_csv(
            output_dir / bigscape_families_table.name, index=False
        )
        df_bigscape_family_presence.to_csv(
            output_dir / bigscape_family_presence_table.name, index=False
        )
        df_bigscape_cluster.to_csv(
            output_dir / bigscape_cluster_table.name, index=False
        )
        df_bigscape_network.to_csv(output_dir / bigscape_network_table.name)

        logging.debug("Copying known compounds table for cutoff {cutoff}...")
        shutil.copy(bigscape_known_table, output_dir / bigscape_known_table.name)

    logging.debug("Copying genome summary table...")
    summary_table = Path(input_dir) / "df_genome_antismash_summary.csv"
    shutil.copy(summary_table, output_dir / summary_table.name)
    logging.info("Job done!")


def get_timestamp(
    input_dir,
    antismash_region_table,
    bigscape_log,
    output_dir,
    cutoffs=["0.30", "0.40", "0.50"],
):
    """
    Identifies the most recent BiG-SCAPE run and corrects its results.

    This function scans the input directory for BiG-SCAPE output files, identifies the most
    recent run based on the timestamp in the filenames, and then calls `correct_bigscape_results`
    to correct the results for that run.

    Parameters:
    - input_dir (Path): Directory containing the BiG-SCAPE output files.
    - antismash_region_table (str): Path to the CSV file containing antiSMASH region data.
    - bigscape_log (str): Path to the BiG-SCAPE log file.
    - output_dir (Path): Directory where corrected output files will be saved.
    - cutoffs (list of str): List of similarity cutoffs used in BiG-SCAPE analysis.

    Returns:
    None
    """
    input_dir = Path(input_dir)
    bigscape_runs = {}
    for table_path in input_dir.glob("*df_*_0.*.csv"):
        logging.debug(f"Found extracted table: {table_path}")
        selected_run_folder = table_path.name
        selected_run_time = selected_run_folder.split("_df")[0]
        selected_run_time = datetime.strptime(selected_run_time, "%Y-%m-%d %H_%M_%S")
        bigscape_runs[selected_run_time] = table_path
    run_times = list(bigscape_runs.keys())
    run_times.sort(reverse=True)

    # make sure run times has values
    assert len(run_times) > 0
    selected_run = run_times[0]

    logging.debug(f"Time stamp: {selected_run}")
    selected_run = str(selected_run).replace(":", "_")
    correct_bigscape_results(
        input_dir,
        antismash_region_table,
        bigscape_log,
        selected_run,
        output_dir,
        cutoffs=cutoffs,
    )


if __name__ == "__main__":
    get_timestamp(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
