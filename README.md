# snakemake-BGC
Snakemake workflow to combine internal &amp; public dataset for downstream analysis

## General steps for the workflow
1. Download public genome fasta files from NCBI using input table public.csv 
2. Download internal genome fasta files from Azure Data Lake using input table internal.csv
3. Run prokka annotation and save output files
4. Run GTDB-tk on internal genomes and download GTDB taxonomy for public genomes
5. Run antiSMASH on prokka generated gbk files
6. (In batch) Run autoMLST-wrapper to get phylogenetic tree
7. (In batch) Run antiSMASH-db import to setup sql server
8. (In batch) Run BiG-SCAPE to get similarity network
9. (In batch) Run BiG-SLICE to get GCFs
10. (In batch) Run BiG-FAM db to setup sql server
11. (In batch) Run Roary to get pangenome presence absence matrix
12. (In batch) Run BGC analytics to get dataframes df_bgcs, df_genomes, df_genes populated with relevant information
13. (In batch) Run BGC analytics to get visualization of pangenome, bgc data
14. (In batch) Run iTOL and save the figures in directories
