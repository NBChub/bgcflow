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
        "workflow/report/logs/ppanggolin/ppanggolin/ppanggolin_{name}_{version}.log"
    threads: 16
    shell:
        """
        ppanggolin panrgp --anno {input} --cpu {threads} --output {output} &>> {log}
        """

rule ppanggolin_genome_roary:
    input:
        gff=lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES),
        roary = "data/interim/roary/{name}"
    output:
        ppanggolin = directory("data/processed/{name}/ppanggolin/genome_roary"),
        panrgp = directory("data/processed/{name}/ppanggolin/panrgp_roary")
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "workflow/report/logs/ppanggolin/genome/ppanggolin_{name}.log"
    threads: 16
    params:
        rarefaction_depth = 30,
        rarefaction_min = 5,
        rarefaction_max = 50
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_GFF="$TMPDIR/prokka_gff.txt"
        ORGANISM_ANNOTATION_LIST="$TMPDIR/ORGANISM_ANNOTATION_LIST.txt"
        CLUSTERS_FILE="$TMPDIR/CLUSTERS_FILE.txt"
        echo '{input.gff}' > $INPUT_GFF

        echo "Preparing input data for ppanggolin..." >> {log}
        python workflow/bgcflow/bgcflow/data/prep_ppanggolin_from_roary.py "$INPUT_GFF" "$ORGANISM_ANNOTATION_LIST" 2>> {log}

        echo "\nConverting Roary output to mmseqs2 format" >> {log}
        python workflow/bgcflow/bgcflow/data/prep_roary_cluster_to_mmseqs2_format.py {input.roary}/clustered_proteins $CLUSTERS_FILE 2>> {log}

        echo "\nRunning ppanggolin annotate..." >> {log}
        ppanggolin annotate --cpu {threads} --anno $ORGANISM_ANNOTATION_LIST --output {output.ppanggolin} &>> {log}

        echo "\nRunning ppanggolin cluster..." >> {log}
        ppanggolin cluster --cpu {threads} -p {output.ppanggolin}/pangenome.h5 --clusters $CLUSTERS_FILE &>> {log}

        echo "\nRunning ppanggolin graph..." >> {log}
        ppanggolin graph -p {output.ppanggolin}/pangenome.h5 --cpu {threads} &>> {log}

        echo "\nRunning ppanggolin partition..." >> {log}
        ppanggolin partition -p {output.ppanggolin}/pangenome.h5 --cpu {threads} &>> {log}

        echo "\nRunning ppanggolin rgp..." >> {log}
        ppanggolin rgp -p {output.ppanggolin}/pangenome.h5 --cpu {threads} &>> {log}

        echo "\nRunning ppanggolin spot..." >> {log}
        ppanggolin spot -p {output.ppanggolin}/pangenome.h5 --cpu {threads} &>> {log}

        echo "\nDrawing rarefaction cureve..." >> {log}
        ppanggolin rarefaction -f \
            --min {params.rarefaction_min} \
            --max {params.rarefaction_max} \
            --depth {params.rarefaction_depth} \
            -p {output.ppanggolin}/pangenome.h5 \
            --cpu {threads} \
            --output {output.ppanggolin} &>> {log}

        echo "\nWriting ppanggolin stats..." >> {log}
        ppanggolin write -f -p {output.ppanggolin}/pangenome.h5 --stats --output {output.panrgp} &>> {log}

        echo "\nWriting ppanggolin regions output..." >> {log}
        ppanggolin write -f -p {output.ppanggolin}/pangenome.h5 --regions --output {output.panrgp} &>> {log}

        echo "\nWriting ppanggolin spots output..." >> {log}
        ppanggolin write -f -p {output.ppanggolin}/pangenome.h5 --spots --output {output.panrgp} &>> {log}

        echo "\nDrawing figures..." >> {log}
        ppanggolin draw -f -p {output.ppanggolin}/pangenome.h5 --ucurve --output {output.ppanggolin} &>> {log}
        ppanggolin draw -f -p {output.ppanggolin}/pangenome.h5 --tile_plot --output {output.ppanggolin} &>> {log}
        ppanggolin draw -f -p {output.ppanggolin}/pangenome.h5 --tile_plot --nocloud --output {output.ppanggolin} &>> {log}
        ppanggolin draw -f -p {output.ppanggolin}/pangenome.h5 --spots all --output {output.panrgp} &>> {log}
        """
