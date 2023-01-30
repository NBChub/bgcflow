def get_single_bgc_input(pep_object, bgc_id, antismash_version):
    """
    Given a PEP Object, get all genbank files based on the bgc_id column
    """
    antismash_path = Path(f"data/interim/antismash/{antismash_version}")
    df = pep_object.sample_tables
    df = df.set_index("bgc_id", drop=False)

    genome_id = df.loc[bgc_id, "genome_id"]
    # override with custom path
    assert 'gbk_path' in df.columns
    custom_path = df.loc[bgc_id, "gbk_path"]
    #print(custom_path, type(custom_path), custom_path != None, file=sys.stderr)
    if custom_path != None:
        gbk_path = custom_path
    else:
        gbk_path = antismash_path / genome_id / f"{bgc_id}.gbk"
    #print(bgc_id, gbk_path, file=sys.stderr)
    return gbk_path

interproscan_version = "5.60-92.0"

rule install_interproscan:
    output:
        interproscan = f"resources/interproscan/interproscan-{interproscan_version}"
    log:
        "workflow/report/logs/interproscan/install_interproscan.log",
    conda:
        "../envs/utilities.yaml"
    params:
        version = interproscan_version
    shell:
        """
        java --version >> {log}
        #sudo apt-get install libpcre3-dev
        wget -P resources https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{params.version}/interproscan-{params.version}-64-bit.tar.gz -nc 2>> {log}
        wget -P resources https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{params.version}/interproscan-{params.version}-64-bit.tar.gz.md5 -nc 2>> {log}
        (cd resources && md5sum -c interproscan-{params.version}-64-bit.tar.gz.md5) 2>> {log}
        (cd resources && tar -pxvzf interproscan-{params.version}-64-bit.tar.gz && mv interproscan-{params.version} interproscan) 2>> {log}
        (cd resources/interproscan/interproscan-{params.version} && python3 setup.py interproscan.properties) 2>> {log}
        (cd resources/interproscan/interproscan-{params.version} && ./interproscan.sh -i test_all_appl.fasta -f tsv -dp) &>> {log}
        (cd resources/interproscan/interproscan-{params.version} && ./interproscan.sh -i test_all_appl.fasta -f tsv) &>> {log}
        """

rule prepare_aa_interproscan:
    input:
        gbk = lambda wildcards: get_single_bgc_input(PEP_PROJECTS[wildcards.name], wildcards.bgc, wildcards.version),
    output:
        fasta = temp("data/interim/interproscan/{bgc}_{name}_{version}.faa"),
    log:
        "workflow/report/logs/interproscan/create_aa_{bgc}_{name}_{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/misc/create_aa.py {input.gbk} {output.fasta} 2>> {log}
        """

rule interproscan:
    input:
        fasta = "data/interim/interproscan/{bgc}_{name}_{version}.faa",
        interproscan = f"resources/interproscan/interproscan-{interproscan_version}"
    output:
        tsv = "data/interim/interproscan/{bgc}_{name}_{version}.faa.tsv",
        json = "data/interim/interproscan/{bgc}_{name}_{version}.faa.json",
    log:
        "workflow/report/logs/interproscan/{bgc}_{name}_{version}.log",
    conda:
        "../envs/utilities.yaml"
    threads: 4
    params:
        appl = "TIGRFAM,PFAM"
    shell:
        """
        {input.interproscan}/interproscan.sh -appl {params.appl} -cpu {threads} -d data/interim/interproscan -f TSV,JSON -i {input.fasta} --verbose &>> {log}
        """
