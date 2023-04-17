rule prepare_aa_mmseqs2:
    input:
        gbk = lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
    output:
        fasta = "data/interim/mmseqs2/{name}/{name}_{version}.faa",
        gbk = "data/interim/minimap2/{name}/{name}_{version}.gbk",
    log:
        "workflow/report/logs/mmseqs2/create_aa_{name}_{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/misc/create_aa.py '{input.gbk}' {output.fasta} 2>> {log}
        python workflow/bgcflow/bgcflow/features/prep_mmseqs2.py '{input.gbk}' {output.gbk} 2>> {log}
        """

rule minimap2:
    input:
        gbk = "data/interim/minimap2/{name}/{name}_{version}.gbk",
    output:
        fna = "data/interim/minimap2/{name}/{name}_{version}.fna",
        paf = "data/interim/minimap2/{name}/{name}_{version}.paf",
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "workflow/report/logs/minimap2/{name}_{version}.log",
    shell:
        """
        any2fasta {input.gbk} > {output.fna} 2>> {log}
        minimap2 -X -N 50 -p 0.1 -c {output.fna} {output.fna} > {output.paf} 2>> {log}
        """

rule mmseqs2_easy_cluster:
    input:
        fasta = "data/interim/mmseqs2/{name}/{name}_{version}.faa",
    output:
        tsv = "data/interim/mmseqs2/{name}/{name}_{version}_cluster.tsv",
    log:
        "workflow/report/logs/mmseqs2/build_{name}_{version}.log",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 4
    shell:
        """
        tmpdir="data/interim/mmseqs2/tmp/{wildcards.name}_{wildcards.version}"
        mkdir -p $tmpdir &>> {log}
        mmseqs easy-cluster {input.fasta} data/interim/mmseqs2/{wildcards.name}/{wildcards.name}_{wildcards.version} $tmpdir -e 1e-5 -c 0.7 &>> {log}
        """

rule mmseqs2:
    input:
        fasta = "data/interim/mmseqs2/{name}/{name}_{version}.faa",
    output:
        db = "data/interim/mmseqs2/{name}/{name}_{version}.db",
        clusterdb = "data/interim/mmseqs2/{name}/cluster_{name}_{version}.db.index",
        dbtype  = "data/interim/mmseqs2/{name}/cluster_{name}_{version}.db.dbtype"
    log:
        "workflow/report/logs/mmseqs2/build_{name}_{version}.log",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 4
    shell:
        """
        mmseqs createdb {input.fasta} {output.db} &>> {log}
        tmpdir="data/interim/mmseqs2/{wildcards.name}/tmp/{wildcards.name}_{wildcards.version}"
        mkdir -p $tmpdir &>> {log}
        mmseqs cluster --threads {threads} {output.db} data/interim/mmseqs2/{wildcards.name}/cluster_{wildcards.name}_{wildcards.version}.db $tmpdir &>> {log}
        """

rule mmseqs2_extract:
    input:
        db = "data/interim/mmseqs2/{name}/{name}_{version}.db",
        clusterdb = "data/interim/mmseqs2/{name}/cluster_{name}_{version}.db.index",
    output:
        edge_table = "data/processed/{name}/mmseqs2/edge_{version}.tsv",
        seqdb = "data/interim/mmseqs2/{name}/seq_{name}_{version}.db.index",
        msa = "data/interim/mmseqs2/{name}/msa_{name}_{version}.db",
    log:
        "workflow/report/logs/mmseqs2/extract_{name}_{version}.log",
    conda:
        "../envs/mmseqs2.yaml"
    threads: 4
    shell:
        """
        mmseqs createtsv --threads {threads} {input.db} {input.db} data/interim/mmseqs2/{wildcards.name}/cluster_{wildcards.name}_{wildcards.version}.db {output.edge_table} &>> {log}
        mmseqs result2msa --threads {threads} {input.db} {input.db} data/interim/mmseqs2/{wildcards.name}/cluster_{wildcards.name}_{wildcards.version}.db {output.msa} --msa-format-mode 3 &>> {log}
        mmseqs createseqfiledb {input.db} data/interim/mmseqs2/{wildcards.name}/cluster_{wildcards.name}_{wildcards.version}.db data/interim/mmseqs2/{wildcards.name}/seq_{wildcards.name}_{wildcards.version}.db &>> {log}
        """

rule mmseq_all:
    input:
        msa = "data/interim/mmseqs2/{name}/msa_{name}_{version}.db",
        edge_table = "data/processed/{name}/mmseqs2/edge_{version}.tsv",
        tsv = "data/interim/mmseqs2/{name}/{name}_{version}_cluster.tsv",
        paf = "data/interim/minimap2/{name}/{name}_{version}.paf",
    output:
        tag = "data/processed/{name}/mmseqs2/{version}.tag",
    shell:
        """
        touch {output.tag}
        """
