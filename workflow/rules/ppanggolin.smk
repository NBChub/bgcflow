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
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n##### 3. Partitioning graph with ppanggolin partition... #####" >> {log}
        ppanggolin partition \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 4. Predicting region of genome plasticity with ppanggolin rgp... ####" >> {log}
        ppanggolin rgp \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 5. Finding spots of insertion with ppanggolin spot... ####" >> {log}
        ppanggolin spot \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 6. Finding conserved modules with ppanggolin module... ####" >> {log}
        ppanggolin module \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --verbose 1 &>> {log}

        echo "\n#### 7. Calculating rarefaction curves with ppanggolin rarefaction... ####" >> {log}
        ppanggolin rarefaction -f \
            --min {params.rarefaction_min} \
            --max {params.rarefaction_max} \
            --depth {params.rarefaction_depth} \
            -p {output.ppanggolin} \
            --cpu {threads} \
            --output {output.rarefaction} &>> {log}

        echo "\n#### 8. Building trees using ppanggolin msa... ####" >> {log}
        ppanggolin msa -f -p {output.ppanggolin} --partition {params.msa_partition} --output {output.msa_core} &>> {log}
        ppanggolin msa -f -p {output.ppanggolin} --phylo --output {output.msa_phylo} &>> {log}
        """

rule ppanggolin_genome_write_regions:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/rarefaction"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/regions"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/write_regions_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --regions --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_write_stats:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/regions"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/stats"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/write_stats_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --stats --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_write_spots:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/stats"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/spots"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/write_spots_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --spots --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_draw_ucurve:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/spots"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/ucurve"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/draw_ucurve_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --ucurve --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_draw_tile_plot:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/ucurve"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/tile_plot"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/tile_plot_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --tile_plot --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_draw_tile_plot_nocloud:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/tile_plot"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/tile_plot_nocloud"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/tile_plot_nocloud_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --tile_plot --nocloud --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_draw_spots:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/tile_plot_nocloud"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/spots_draw"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/draw_spots_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --spots all --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_gexf:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/spots_draw"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/gexf"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/gexf_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --light_gexf --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --gexf --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --json --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_gene_pres_abs:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/gexf"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/gene_pres_abs"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/gene_pres_abs_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --Rtab --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --csv --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_partitions:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/gene_pres_abs"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/partitions"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/partitions_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --partitions --output $(dirname {input.ppanggolin}) &>> {log}
        """

rule ppanggolin_genome_projection:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/partitions"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/projection"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/projection_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --projection --output $(dirname {input.ppanggolin}) &>> {log}
        """

rule ppanggolin_genome_families:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/projection"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/families"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/families_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --families_tsv --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_borders:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/families"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/borders"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/borders_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --borders --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_modules:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/borders"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/modules"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/modules_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --modules --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_spot_modules:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome/modules"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome/spot_modules"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome/spot_modules_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --spot_modules --output {output.folder} &>> {log}
        """
