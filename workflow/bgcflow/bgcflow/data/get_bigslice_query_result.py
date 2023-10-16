import logging
import shutil
import sqlite3
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def get_bigslice_query(
    project_name, output_folder, result_folder="resources/bigslice/full_run_result/"
):
    bigslice_full_run = Path(result_folder)
    query_id = find_query_result(project_name, bigslice_full_run)
    fetch_query_tables(query_id, output_folder, bigslice_full_run)
    logging.info("Copying SQL database of the run...")
    try:
        shutil.copy(
            bigslice_full_run / f"reports/{query_id}/data.db",
            Path(output_folder) / f"{query_id}.db",
        )
    except BlockingIOError:
        subprocess.run(
            ["cp", str(bigslice_full_run / f"reports/{query_id}/data.db"), str(Path(output_folder) / f"{query_id}.db")]
        )
    logging.info("Job done!")
    return


def find_query_result(project_name, bigslice_full_run):
    logging.info("Connecting to BiG-SLICE SQL Reports...")
    conn = sqlite3.connect(bigslice_full_run / "reports/reports.db")
    df = pd.read_sql_query("select * from reports;", conn)

    # select only run result of the project
    match = []
    for i in df.index:
        if project_name in df.loc[i, "name"]:
            match.append(i)
    df = df.loc[match, :]
    logging.debug(f"Found resulting match to project {project_name}:\n{df}")

    # select the latest run for the project
    df["creation_date"] = pd.to_datetime(df["creation_date"])
    filter_time = df.creation_date == df.creation_date.max()
    df = df.where(filter_time).dropna()
    logging.debug(f"Extracting information from the latest run:{df.name.values}")
    return int(df.id.values)


def fetch_query_tables(query_id, output_folder, bigslice_full_run):
    # Connect to database
    conn = sqlite3.connect(bigslice_full_run / f"reports/{query_id}/data.db")

    # get all table name from sql
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    sql_table_names = [i[0] for i in cursor.fetchall()]

    output_folder = Path(output_folder)

    try:
        output_folder.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        logging.debug(f"Directory {output_folder} already exist")
    else:
        logging.debug(f"Creating directory: {output_folder}...")

    logging.debug(f"Extracting tables to {output_folder}...")

    # convert table to pandas df
    df_tables = {}
    for t in sql_table_names:
        df = pd.read_sql_query(f"select * from {t};", conn)
        outfile = output_folder / f"{t}.csv"
        df.to_csv(outfile, index=False)
        df_tables.update({t: df})
    return


if __name__ == "__main__":
    get_bigslice_query(sys.argv[1], sys.argv[2], sys.argv[3])
