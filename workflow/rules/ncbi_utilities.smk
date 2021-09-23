#############################
# Use genomes from NCBI

## set up strains
df = pd.read_csv(config["ncbi"], sep="\t")

# filter duplicates
df['refseq'] = df[['strain', 'assembly','refseq']].groupby(['strain', 'assembly'])['refseq'].transform(lambda x: ','.join(x))
ncbi_genomes = df[['strain', 'assembly','refseq']].drop_duplicates().reset_index()

## generate taxonomy
for num, tax in enumerate(ncbi_genomes.strain.to_list()):
    taxonomy = tax.split(" ",2)
    ncbi_genomes.loc[num, "genus"] = taxonomy[0]
    ncbi_genomes.loc[num, "species"] = taxonomy[1]
    ncbi_genomes.loc[num, "strain_info"] = taxonomy[-1]
    
## fix indexing
ncbi_genomes.set_index("assembly", drop=False)
ncbi_genomes.index.names = ["assembly_id"]

## watermark snakemake version to samples
for i in ncbi_genomes.index:
    ncbi_genomes.loc[i, "assembly_snakemake"] = str(ncbi_genomes.loc[i, "assembly"]) + __version__

#validate
validate(ncbi_genomes, schema="../schemas/ncbi.schema.yaml")

wildcard_constraints:
    assembly="|".join(ncbi_genomes.assembly),

# wildcards
NCBI_GENOMES = ncbi_genomes.assembly.to_list()

#############################
rule ncbi_genome_download:
    output:
        "data/interim/fasta/{assembly}.fna",
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s refseq -F fasta -A {wildcards.assembly} -o data/raw/ncbi/download bacteria --verbose
        gunzip -c data/raw/ncbi/download/refseq/bacteria/{wildcards.assembly}/*.fna.gz > {output}
        #rm -rf data/raw/ncbi/download/refseq/bacteria/{wildcards.assembly}
        """

rule prokka_ncbi:
    input: 
        fna = "data/interim/fasta/{assembly}.fna",
    output:
        gff = "data/interim/prokka/{assembly}/{assembly}.gff",
        gbk = "data/interim/prokka/{assembly}/{assembly}.gbk",
        genus = temp("data/interim/{assembly}/genus"),
        species = temp("data/interim/{assembly}/species")
    conda:
        "../envs/prokka.yaml"
    threads: 12
    shell:
        """
        head -1 {input} | cut -d' ' -f2 > {output.genus}
        head -1 {input} | cut -d' ' -f3 > {output.species}
        prokka --outdir data/interim/prokka/{wildcards.assembly} --force --prefix {wildcards.assembly} --genus `cat genus` --species `cat species` --strain {wildcards.assembly} --cdsrnaolap --cpus {threads} --rnammer --increment 10 --evalue 1e-05 {input.fna}
        """