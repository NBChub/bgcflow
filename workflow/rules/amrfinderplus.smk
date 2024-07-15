rule amrfinderplus:
    input:
        fna="data/interim/fasta/{strains}.fna",
        faa="data/interim/prokka/{strains}/{strains}.faa",
        gff="data/interim/prokka/{strains}/{strains}.gff",
    output:
        table="data/interim/amrfinderplus/{strains}.tsv",
    log:
        "logs/amrfinderplus/amrfinderplus/{strains}.log"
    conda:
        "../envs/amrfinderplus.yaml"
    params:
        annotation="prokka",
        ident_min="-1",
        coverage_min="0.5",
        translation_table="11",
    threads: 4
    shell:
        """
        amrfinder -p {input.faa} \
            -n {input.fna} \
            -g {input.gff} \
            -a {params.annotation} \
            --plus \
            --ident_min {params.ident_min} \
            --coverage_min {params.coverage_min} \
            --translation_table {params.translation_table} \
            --threads {threads} \
            --log {log} > {output.table}
        """

rule amrfinder_gather:
    input:
        tables=lambda wildcards: expand(
            "data/interim/amrfinderplus/{strains}.tsv",
            strains=[s for s in PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique()],
        ),
    output:
        table="data/processed/{name}/tables/df_amrfinderplus.csv"
    conda:
        "../envs/bgc_analytics.yaml"
    log:
        "logs/amrfinderplus/gather/amrfinderplus_gather_{name}.log"
    shell:
        """
        TMPDIR="data/interim/tmp/{wildcards.name}"
        mkdir -p $TMPDIR
        INPUT_TSV="$TMPDIR/df_amrfinderplus.txt"
        echo '{input.tables}' > $INPUT_TSV
        python workflow/bgcflow/bgcflow/data/gather_amrfinder.py $INPUT_TSV {output.table} 2>> {log}
        rm $INPUT_TSV
        """
