interproscan_version = "5.60-92.0"

rule install_interproscan:
    output:
        interproscan = directory(f"resources/interproscan/interproscan-{interproscan_version}")
    log:
        "logs/interproscan/install_interproscan.log",
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
        if [ ! -d "{output.interproscan}" ]; then
            mkdir -p resources/interproscan
            (cd resources && md5sum -c interproscan-{params.version}-64-bit.tar.gz.md5) &>> {log}
            (cd resources && tar -pxvzf interproscan-{params.version}-64-bit.tar.gz && mv interproscan-{params.version} interproscan/) &>> {log}
        else
            echo "InterProScan is extracted in {output.interproscan}." >> {log}
        fi
        (cd resources/interproscan/interproscan-{params.version} && python3 setup.py interproscan.properties) 2>> {log}
        (cd resources/interproscan/interproscan-{params.version} && ./interproscan.sh -i test_all_appl.fasta -f tsv -dp) &>> {log}
        (cd resources/interproscan/interproscan-{params.version} && ./interproscan.sh -i test_all_appl.fasta -f tsv) &>> {log}
        """

rule prepare_aa_interproscan:
    input:
        gbk = lambda wildcards: get_bgc_inputs(PEP_PROJECTS[wildcards.name], wildcards.version),
    output:
        fasta = temp("data/interim/interproscan/{name}_{version}.faa"),
    log:
        "logs/interproscan/create_aa_{name}_{version}.log",
    conda:
        "../envs/bgc_analytics.yaml"
    shell:
        """
        python workflow/bgcflow/bgcflow/misc/create_aa.py '{input.gbk}' {output.fasta} 2>> {log}
        """

rule interproscan:
    input:
        fasta = "data/interim/interproscan/{name}_{version}.faa",
        interproscan = f"resources/interproscan/interproscan-{interproscan_version}"
    output:
        tsv = "data/processed/{name}/tables/interproscan_as_{version}.tsv",
        json = "data/processed/{name}/tables/interproscan_as_{version}.json",
    log:
        "logs/interproscan/{name}_{version}.log",
    conda:
        "../envs/utilities.yaml"
    threads: 4
    params:
        appl = "TIGRFAM,PFAM"
    shell:
        """
        {input.interproscan}/interproscan.sh -appl {params.appl} -cpu {threads} -d data/interim/interproscan -f TSV,JSON -i {input.fasta} --verbose &>> {log}
        mv data/interim/interproscan/{wildcards.name}_{wildcards.version}.faa.tsv {output.tsv}
        mv data/interim/interproscan/{wildcards.name}_{wildcards.version}.faa.json {output.json}
        """
