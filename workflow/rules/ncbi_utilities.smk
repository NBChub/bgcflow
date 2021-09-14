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
        "results/genomes/{assembly}{version}/prokka_input.fna",
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        ncbi-genome-download -s refseq -F fasta -A {wildcards.assembly} -o results/NCBI/ bacteria --verbose
        gunzip -c results/NCBI/refseq/bacteria/{wildcards.assembly}/*.fna.gz > {output}
        rm -rf results/NCBI/refseq/bacteria/{wildcards.assembly}/
        """

rule prokka_ncbi:
    input: 
        "results/genomes/{assembly}{version}/prokka_input.fna",
    output:
        directory("results/genomes/{assembly}{version}/{assembly}{version}_prokka_actinoannotPFAM")
    conda:
        "../envs/prokka.yaml"
    threads: 12
    shell:
        """
        head -1 {input} | cut -d' ' -f2 > genus
        head -1 {input} | cut -d' ' -f3 > species
        prokka --outdir {output} --prefix {wildcards.assembly}{wildcards.version} --genus `cat genus` --species `cat species` --strain {wildcards.assembly} --cdsrnaolap --cpus {threads} --rnammer --increment 10 --evalue 1e-05 {input[0]}
        """