import sqlite3
import pandas as pd
from pathlib import Path
import shutil
import sys

def get_bigslice_query(project_name, output_folder, result_folder="resources/bigslice/full_run_result/"):
    bigslice_full_run = Path(result_folder)
    query_id = find_query_result(project_name, bigslice_full_run)
    fetch_query_tables(query_id, output_folder, bigslice_full_run)
    print("copying sql database...")
    shutil.copy(bigslice_full_run / f"reports/{query_id}/data.db", Path(output_folder) / f"{query_id}.db")
    return

def find_query_result(project_name, bigslice_full_run):
    conn = sqlite3.connect(bigslice_full_run / "reports/reports.db")
    df = pd.read_sql_query(f"select * from reports;", conn)

    # select only run result of the project
    filter_name = df.name == project_name
    df = df.where(filter_name).dropna()
    
    # select the latest run for the project
    df['creation_date'] = pd.to_datetime(df['creation_date'])
    filter_time = df.creation_date == df.creation_date.max()
    df = df.where(filter_time).dropna()
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
        print(f"Directory {output_folder} already exist")
    else:
        print(f"Creating directory: {output_folder}...")
    
    print(f"Writing results to {output_folder}...")
    
    # convert table to pandas df
    df_tables = {}
    for t in sql_table_names:
        df = pd.read_sql_query(f"select * from {t};", conn)     
        outfile = output_folder / f"{t}.csv"
        df.to_csv(outfile, index=False)
        df_tables.update({t : df})
    return

    
if __name__ == "__main__":
    get_bigslice_query(sys.argv[1], sys.argv[2], sys.argv[3])