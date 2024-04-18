import argparse
import duckdb
import logging

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def upload_to_motherduck(database_path, database_name):
    """
    Connects to a local DuckDB database, loads a remote DuckDB database, 
    and creates or replaces a remote database from the current local database.

    Args:
        database_path (str): The path to the local DuckDB database.
        database_name (str): The name of the remote database to create or replace.
    """
    logging.info('Connecting to the local DuckDB database...')
    local_conn = duckdb.connect(database_path)

    logging.info('Loading the remote DuckDB database...')
    local_conn.execute("LOAD 'motherduck'")

    logging.info('Creating or replacing the remote database from the current local database...')
    local_conn.execute(f"CREATE OR REPLACE DATABASE {database_name} FROM CURRENT_DATABASE()")

    logging.info('Closing the connection...')
    local_conn.close()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description='Upload a local DuckDB database to a remote DuckDB database.')
    parser.add_argument('database_path', type=str, help='The path to the local DuckDB database.')
    parser.add_argument('--database_name', type=str, default='bgcflow', help='The name of the remote database to create or replace. Default is "bgcflow".')
    args = parser.parse_args()

    upload_to_motherduck(args.database_path, args.database_name)
