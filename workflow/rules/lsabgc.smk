rule install_lsabgc_db:
    output:
        directory("resources/lsaBGC-database")
    log:
        "logs/lsabgc/install-db.log"
    conda:
        "../envs/lsabgc.yaml"
    params:
        kofam="--no_ko"
    shell:
        """
        mkdir -p {output}
        setup_annotation_dbs.py -p {output} {params.kofam} &>> {log}
        """

rule lsabgc_prepare:
    input:
        fna=lambda wildcards: get_fasta_inputs(wildcards.name, DF_SAMPLES),
        bgc_mapping="data/interim/bgcs/{name}/{name}_antismash_{version}.csv",
    output:
        Primary_Genomes_Listing = "data/interim/lsabgc/{name}/as_{version}/Primary_Genomes_Listing.txt",
        Primary_Genome_BGC_Genbanks_Listing = "data/interim/lsabgc/{name}/as_{version}/Primary_Genome_BGC_Genbanks_Listing.txt",
    log:
        "logs/lsabgc/lsabgc_prepare-{name}-as_{version}.log",
    run:
        mapping_path = Path(input.bgc_mapping)
        fasta_dir = Path("data/interim/fasta")
        df = pd.read_csv(mapping_path)
        for i in df.index:
            genome_id = df.loc[i, "genome_id"]
            bgc_id = df.loc[i, "bgc_id"]
            bgc_path = mapping_path.parent / wildcards.version / genome_id / (str(bgc_id) + ".gbk")
            fna_path = fasta_dir / f"{genome_id}.fna"
            df.loc[i, "bgc_path"] = str(bgc_path)
            df.loc[i, "fna_path"] = str(fna_path)
            assert bgc_path.is_file(), f"File {bgc_path} does not exist"
            assert fna_path.is_file(), f"File {fna_path} does not exist"
        df.loc[:, ["genome_id", "bgc_path"]].to_csv(output.Primary_Genome_BGC_Genbanks_Listing, sep="\t", index=False, header=False)
        df.loc[:, ["genome_id", "fna_path"]].drop_duplicates(subset=["genome_id"]).to_csv(output.Primary_Genomes_Listing, sep="\t", header=False, index=False)

rule lsabgc_ready:
    input:
        database="resources/lsaBGC-database",
        Primary_Genomes_Listing = "data/interim/lsabgc/{name}/as_{version}/Primary_Genomes_Listing.txt",
        Primary_Genome_BGC_Genbanks_Listing = "data/interim/lsabgc/{name}/as_{version}/Primary_Genome_BGC_Genbanks_Listing.txt",
        bigscape_index="data/interim/bigscape/{name}_antismash_{version}/index.html",
    output:
        lsabgc_output=directory("data/processed/{name}/lsabgc/as_{version}/lsaBGC_Ready_Results")
    conda:
        "../envs/lsabgc.yaml"
    threads: 16
    log:
        "logs/lsabgc/lsabgc_ready-{name}-as_{version}.log",
    params:
        annotate="" #"--annotate",
    shell:
        """
        bigscape_dir=$(dirname {input.bigscape_index})
        lsaBGC-Ready.py --genome_listing {input.Primary_Genomes_Listing} \
            --bgc_genbank_listing {input.Primary_Genome_BGC_Genbanks_Listing} \
            --bgc_prediction_software antiSMASH \
            --bigscape_results $bigscape_dir \
            --cpus {threads} \
            {params.annotate} --run_gtotree --lsabgc_cluster --lsabgc_expansion \
            --output_directory {output.lsabgc_output} 2>> {log}
        """

rule lsabgc_prepare_tax:
    input:
        tax="data/processed/{name}/tables/df_gtdb_meta.csv"
    output:
        tax="data/interim/lsabgc/{name}/as_{version}/Genome_to_Species_Mapping.txt",
    log:
        "logs/lsabgc/lsabgc_prepare_tax-{name}-as_{version}.log",
    run:
        df = pd.read_csv(input.tax)
        df.loc[:, ["genome_id", "Organism"]].to_csv(output.tax, sep="\t", index=False, header=False)

rule lsabgc_autoanalyze:
    input:
        tax="data/interim/lsabgc/{name}/as_{version}/Genome_to_Species_Mapping.txt",
        lsabgc_ready="data/processed/{name}/lsabgc/as_{version}/lsaBGC_Ready_Results"
    output:
        lsabgc_output=directory("data/processed/{name}/lsabgc/as_{version}/lsaBGC_AutoAnalyze_Results")
    conda:
        "../envs/lsabgc.yaml"
    threads: 16
    log:
        "logs/lsabgc/lsabgc_autoanalyze-{name}-as_{version}.log",
    shell:
        """
        lsaBGC-AutoAnalyze.py \
            -i {input.lsabgc_ready}/Final_Results/Primary_Sample_Annotation_Files.txt \
            -g {input.lsabgc_ready}/Final_Results/GCF_Listings/ \
            -m {input.lsabgc_ready}/Final_Results/Orthogroups.tsv \
            -s {input.lsabgc_ready}/Final_Results/GToTree_output.tre \
            -w {input.lsabgc_ready}/Final_Results/GToTree_Expected_Similarities.txt \
            -u {input.tax} -c {threads} -o {output} 2>> {log}
        """
