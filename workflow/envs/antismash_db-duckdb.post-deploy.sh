set -e

# Create a directory for resources if it doesn't already exist
mkdir -p resources
rm -rf resources/antismash_db-schema_duckdb

# create logs directory
mkdir -p logs/antismash_db-duckdb
LOG="logs/antismash_db-duckdb/antismash_db-duckdb_template.log"

# Clone the antismash_db-schema_duckdb repository into the resources directory
git clone https://github.com/NBChub/antismash_db-schema_duckdb.git resources/antismash_db-schema_duckdb 2>> $LOG

# Clone the db-schema repository into the antismash_db-schema_duckdb directory
git clone https://github.com/antismash/db-schema.git resources/antismash_db-schema_duckdb/db-schema 2>> $LOG

# Download the antiSMASH databases using a custom script or command
download-antismash-databases 2>> $LOG

# Initialize the DuckDB database with the schema from the db-schema directory
(cd resources/antismash_db-schema_duckdb && python init_duckdb.py db-schema duckdb-schema) 2>> $LOG

# Download the NCBI taxonomy dump files into the ncbi-taxdump directory
(cd resources/antismash_db-schema_duckdb && wget -P ncbi-taxdump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -nc) 2>> $LOG

# Extract the downloaded NCBI taxonomy dump files
(cd resources/antismash_db-schema_duckdb/ncbi-taxdump && tar -xvf new_taxdump.tar.gz) 2>> $LOG

# Install the asdb-taxa tool using cargo (Rust's package manager)
cargo install asdb-taxa 2>> $LOG

# Add cargo's bin directory to the PATH to ensure asdb-taxa can be executed
export PATH="$HOME/.cargo/bin:$PATH"

# Clone the db-import repository into the antismash_db-schema_duckdb directory
git clone git@github.com:matinnuhamunada/db-import.git resources/antismash_db-schema_duckdb/db-import 2>> $LOG

# Checkout a specific branch of the db-import repository
(cd resources/antismash_db-schema_duckdb/db-import && git checkout -b v4.0.0-duckdb v4.0.0-duckdb) 2>> $LOG
