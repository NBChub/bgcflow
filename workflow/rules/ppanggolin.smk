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

rule ppanggolin_genome:
    input:
        gff=lambda wildcards: get_prokka_outputs(wildcards.name, DF_SAMPLES),
        roary = "data/interim/roary/{name}"
    output:
        ppanggolin = directory("data/processed/{name}/ppanggolin/genome")
    conda:
        "../envs/ppanggolin.yaml"
    log:
        "workflow/report/logs/ppanggolin/genome/ppanggolin_{name}.log"
    threads: 16
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_GFF="$TMPDIR/prokka_gff.txt"
        ORGANISM_ANNOTATION_LIST="$TMPDIR/ORGANISM_ANNOTATION_LIST.txt"
        CLUSTERS_FILE="$TMPDIR/CLUSTERS_FILE.txt"
        echo '{input.gff}' > $INPUT_GFF
        python workflow/bgcflow/bgcflow/data/prep_ppanggolin_from_roary.py "$INPUT_GFF" "$ORGANISM_ANNOTATION_LIST" 2>> {log}
        python workflow/bgcflow/bgcflow/data/prep_roary_cluster_to_mmseqs2_format.py {input.roary}/clustered_proteins $CLUSTERS_FILE 2>> {log}
        ppanggolin annotate --cpu {threads} --anno $ORGANISM_ANNOTATION_LIST --output {output} &>> {log}
        ppanggolin cluster --cpu {threads} -p {output}/pangenome.h5 --clusters $CLUSTERS_FILE &>> {log}
        """
