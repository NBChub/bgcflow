set -e

# Create a directory for resources if it doesn't already exist
mkdir -p resources
rm -rf resources/antismash_db-schema_duckdb

# Clone the antismash_db-schema_duckdb repository into the resources directory
git clone https://github.com/NBChub/antismash_db-schema_duckdb.git resources/antismash_db-schema_duckdb

# Clone the db-schema repository into the antismash_db-schema_duckdb directory
git clone https://github.com/antismash/db-schema.git resources/antismash_db-schema_duckdb/db-schema

# Download the antiSMASH databases using a custom script or command
download-antismash-databases

# Initialize the DuckDB database with the schema from the db-schema directory
(cd resources/antismash_db-schema_duckdb && python init_duckdb.py db-schema duckdb-schema)

# Download the NCBI taxonomy dump files into the ncbi-taxdump directory
(cd resources/antismash_db-schema_duckdb && wget -P ncbi-taxdump https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -nc)

# Extract the downloaded NCBI taxonomy dump files
(cd resources/antismash_db-schema_duckdb/ncbi-taxdump && tar -xvf new_taxdump.tar.gz)

# Clone the db-import repository into the antismash_db-schema_duckdb directory
git clone git@github.com:matinnuhamunada/db-import.git resources/antismash_db-schema_duckdb/db-import

# Checkout a specific branch of the db-import repository
(cd resources/antismash_db-schema_duckdb/db-import && git checkout -b v4.0.0-duckdb v4.0.0-duckdb)
