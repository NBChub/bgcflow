# Start your workflow at - *raw* data

To start `bgcflow` analysis, drop the fasta files of all genomes from your project in the directory called `fasta`. 

## Notes on formatting of fasta input files

Make sure that each fasta file has unique genome ID including version of assembly as the filename with `.fna` as extension. For example, NCBI refseq genomes shall be in format - `GCF_000000000.1.fna`, internal genomes shall be `NBC00000.1.fna` format, and PATRIC genomes can be in `100226.137.fna` format

Make sure that accesion IDs of each contig are uniquely represented along with full name of organism in the fasta header. For example:

1. `>NZ_CP042324.1 Streptomyces coelicolor A3(2) strain CFB_NBC_0001 chromosome, complete genome`
2. `>NBC00000.1.1 Streptomyces sp. strain NBC00000 chromosome, complete genome`
3. `>CP042324   CP042324.1   [Streptomyces coelicolor A3(2) strain CFB_NBC_0001 | 100226.137]`

## Direct NCBI download

In case you want to analyze public genomes at NCBI RefSeq database or internally sequenced genomes, please provide a tab separated text file in the configuration folder `config/ncbi.tsv` or `config/internal.tsv`. Both of these files must have a column named `genome_id` with with indexes as NCBI RefSeq genome accession IDs or internal genome IDs with versions, for exmaple `GCF_000000000.1` or `NBC00000.1`. The remaining columns on metadata of the strains can also be left empty. We recommend to add `organism`, `genus`, `species`, `strain` as the four additional columns.

If you want to analyze further customized options for starting dataset apart from above please contact Omkar Mohite (omkmoh@biosustain.dtu.dk) or Matin Nuhamunada (matinnu@biosustain.dtu.dk).