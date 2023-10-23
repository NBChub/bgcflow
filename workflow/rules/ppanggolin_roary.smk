rule ppanggolin_genome_roary:
    input:
        gff=lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES),
        roary = "data/interim/roary/{name}"
    output:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        rarefaction = directory("data/processed/{name}/ppanggolin/genome_roary/rarefaction"),
        msa_core = directory("data/processed/{name}/ppanggolin/genome_roary/msa_core"),
        msa_phylo = directory("data/processed/{name}/ppanggolin/genome_roary/msa_phylo"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/ppanggolin_roary_{name}.log"
    threads: 16
    params:
        cluster_identity = 0.5,
        cluster_coverage = 0.8,
        rarefaction_depth = 30,
        rarefaction_min = 5,
        rarefaction_max = 50,
        msa_partition = "core",
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        OUTDIR="data/processed/{wildcards.name}/ppanggolin/genome"
        mkdir -p $OUTDIR
        INPUT_GFF="$TMPDIR/prokka_gff.txt"
        ORGANISM_ANNOTATION_LIST="$TMPDIR/ORGANISM_ANNOTATION_LIST.txt"
        CLUSTERS_FILE="$TMPDIR/CLUSTERS_FILE.txt"
        echo '{input.gff}' > $INPUT_GFF

        echo "Preparing input data for ppanggolin..." >> {log}
        python workflow/bgcflow/bgcflow/data/prep_ppanggolin_from_roary.py "$INPUT_GFF" "$ORGANISM_ANNOTATION_LIST" 2>> {log}

        echo "\nConverting Roary output to mmseqs2 format" >> {log}
        python workflow/bgcflow/bgcflow/data/prep_roary_cluster_to_mmseqs2_format.py {input.roary}/clustered_proteins $CLUSTERS_FILE 2>> {log}

        echo "\n##### 1 Running ppanggolin annotate... #####" >> {log}s
        ppanggolin annotate \
            --cpu {threads} \
            --anno $ORGANISM_ANNOTATION_LIST \
            --output $OUTDIR \
            --force \
            --verbose 1 &>> {log}

        echo "\n##### 2 Running ppanggolin cluster... #####" >> {log}
        ppanggolin cluster \
            --cpu {threads} \
            --identity {params.cluster_identity} \
            --coverage {params.cluster_coverage} \
            -p {output.ppanggolin}/pangenome.h5 \
            --clusters $CLUSTERS_FILE \
            --infer_singletons \
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

rule ppanggolin_genome_roary_write_regions:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/rarefaction"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/regions"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/write_regions_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --regions --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_write_stats:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/regions"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/stats"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/write_stats_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --stats --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_write_spots:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/stats"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/spots"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/write_spots_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --spots --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_draw_ucurve:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/spots"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/ucurve"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/draw_ucurve_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --ucurve --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_draw_tile_plot:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/ucurve"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/tile_plot"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/tile_plot_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --tile_plot --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_draw_tile_plot_nocloud:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/tile_plot"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/tile_plot_nocloud"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/tile_plot_nocloud_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --tile_plot --nocloud --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_draw_spots:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/tile_plot_nocloud"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/spots_draw"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/draw_spots_{name}.log"
    shell:
        """
        ppanggolin draw -f -p {input.ppanggolin} --spots all --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_gexf:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/spots_draw"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/gexf"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/gexf_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --light_gexf --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --gexf --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --json --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_gene_pres_abs:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/gexf"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/gene_pres_abs"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/gene_pres_abs_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --Rtab --output {output.folder} &>> {log}
        ppanggolin write -f -p {input.ppanggolin} --csv --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_partitions:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/gene_pres_abs"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/partitions"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/partitions_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --partitions --output $(dirname {input.ppanggolin}) &>> {log}
        """

rule ppanggolin_genome_roary_projection:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/partitions"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/projection"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/projection_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --projection --output $(dirname {input.ppanggolin}) &>> {log}
        """

rule ppanggolin_genome_roary_families:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/projection"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/families"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/families_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --families_tsv --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_borders:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/families"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/borders"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/borders_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --borders --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_modules:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/borders"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/modules"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/modules_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --modules --output {output.folder} &>> {log}
        """

rule ppanggolin_genome_roary_spot_modules:
    input:
        ppanggolin = "data/processed/{name}/ppanggolin/genome_roary/pangenome.h5",
        previous = "data/processed/{name}/ppanggolin/genome_roary/modules"
    output:
        folder = directory("data/processed/{name}/ppanggolin/genome_roary/spot_modules"),
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "logs/ppanggolin/genome_roary/spot_modules_{name}.log"
    shell:
        """
        ppanggolin write -f -p {input.ppanggolin} --spot_modules --output {output.folder} &>> {log}
        """
