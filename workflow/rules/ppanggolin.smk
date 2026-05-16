rule ppanggolin_bgc_prep:
    input:
        antismash_bgc = "data/interim/bgcs/{name}/{version}"
    output:
        table = "data/interim/ppanggolin/BGC/dataset/{name}_{version}.tsv",
    run:
        bgc = Path(input.antismash_bgc)
        df = pd.DataFrame([[i.name, i / f"{i.name}.gbk"] for i in bgc.glob("*") if not i.name.startswith(".snakemake")])
        df.to_csv(output.table, sep="\t", index=False, header=False)

rule ppanggolin_BGC:
    input:
        table = "data/interim/ppanggolin/BGC/dataset/{name}_{version}.tsv",
    output:
        ppanggolin = directory("data/ppanggolin/BGC/{name}_{version}")
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/ppanggolin/ppanggolin_{name}_{version}.log"
    threads: 16
    shell:
        """
        ppanggolin panrgp --anno {input} --cpu {threads} --output {output} &>> {log}
        """

rule ppanggolin_genome:
    input:
        gff=lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES),
    output:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        rarefaction = directory("data/processed/{name}/ppanggolin/genome/rarefaction"),
        msa_core = directory("data/processed/{name}/ppanggolin/genome/msa_core"),
        msa_phylo = directory("data/processed/{name}/ppanggolin/genome/msa_phylo"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/ppanggolin_{name}.log"
    threads: 16
    params:
        cluster_identity = 0.8,
        cluster_coverage = 0.8,
        rarefaction_depth = 10,
        rarefaction_min = 5,
        rarefaction_max = 30,
        msa_partition = "core",
    shell:
        """
        TMPDIR="data/interim/tmp/ppanggolin/{wildcards.name}"
        mkdir -p $TMPDIR
        OUTDIR="data/processed/{wildcards.name}/ppanggolin/genome"
        mkdir -p $OUTDIR
        INPUT_GFF="$TMPDIR/prokka_gff.txt"
        ORGANISM_ANNOTATION_LIST="$TMPDIR/ORGANISM_ANNOTATION_LIST.txt"
        echo '{input.gff}' > $INPUT_GFF

        echo "Preparing input data for ppanggolin..." >> {log}
        python workflow/bgcflow/bgcflow/data/prep_ppanggolin_from_roary.py "$INPUT_GFF" "$ORGANISM_ANNOTATION_LIST" 2>> {log}

        echo "\n##### 1. Running ppanggolin annotate... ##### " >> {log}
        ppanggolin annotate \
            --anno $ORGANISM_ANNOTATION_LIST \
            --output $OUTDIR \
            --cpu {threads} \
            --force \
            --verbose 1 &>> {log}

        echo "\n##### 2. Building gene families with ppanggolin cluster... #####" >> {log}
        ppanggolin cluster \
            -p {output.ppanggolin} \
            --identity {params.cluster_identity} \
            --coverage {params.cluster_coverage} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n##### 3. Building pangenome graph with ppanggolin graph... #####" >> {log}
        ppanggolin graph \
            -p {output.ppanggolin} \
            --verbose 1 &>> {log}

        echo "\n##### 4. Partitioning graph with ppanggolin partition... #####" >> {log}
        ppanggolin partition \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 5. Predicting region of genome plasticity with ppanggolin rgp... ####" >> {log}
        ppanggolin rgp \
            -p {output.ppanggolin} \
            --verbose 1 &>> {log}

        echo "\n#### 6. Finding spots of insertion with ppanggolin spot... ####" >> {log}
        ppanggolin spot \
            -p {output.ppanggolin} \
            --verbose 1 &>> {log}

        echo "\n#### 7. Finding conserved modules with ppanggolin module... ####" >> {log}
        ppanggolin module \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 8. Calculating rarefaction curves with ppanggolin rarefaction... ####" >> {log}
        ppanggolin rarefaction -f \
            --min {params.rarefaction_min} \
            --max {params.rarefaction_max} \
            --depth {params.rarefaction_depth} \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --output {output.rarefaction} &>> {log}

        echo "\n#### 9. Building trees using ppanggolin msa... ####" >> {log}
        ppanggolin msa -f -p {output.ppanggolin} --partition {params.msa_partition} --output {output.msa_core} &>> {log}
        ppanggolin msa -f -p {output.ppanggolin} --phylo --output {output.msa_phylo} &>> {log}
        """

rule ppanggolin_write_pangenome:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/msa_phylo"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/write_pangenome"),
    conda:
        "../envs/ppanggolin.yaml"
    threads: 16
    log:
        "logs/ppanggolin/genome/write_pangenome_{name}.log"
    shell:
        """
        ppanggolin write_pangenome -f -p {input.ppanggolin} \
            --gexf \
            --light_gexf \
            --json \
            --csv \
            --stats \
            --partitions \
            --families_tsv \
            --regions \
            --regions_families \
            --spots \
            --borders \
            --modules \
            --spot_modules \
            --cpu {threads} \
            --output {output.folder} &>> {log}
        """

rule ppanggolin_write_genomes:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/write_pangenome"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/write_genomes"),
    conda:
        "../envs/ppanggolin.yaml"
    threads: 16
    log:
        "logs/ppanggolin/genome/write_genomes_{name}.log"
    shell:
        """
        ppanggolin write_genomes -f -p {input.ppanggolin} \
            --cpu {threads} \
            --proksee \
            --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_draw:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/write_genomes"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/draw"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/draw_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} \
            --tile_plot \
            --ucurve \
            --draw_spots \
            --nocloud \
            --add_dendrogram \
            --add_metadata \
            --output {output.folder} &>> {log}
        """
