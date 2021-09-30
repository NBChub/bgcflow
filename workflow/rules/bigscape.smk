rule install_bigscape:
    output:
        directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip)
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/)
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz)
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm)
        """

rule bigscape:
    input: 
        "resources/BiG-SCAPE/",
        expand("data/interim/antismash/{strains}", strains=STRAINS)
    output:
        directory("data/interim/bigscape/all/")
    conda:
        "../envs/bigscape.yaml"
    params:
        label = "all",
    log:
        "workflow/report/bigscape.log"
    threads: 12
    shell:
        """
        python bigscape.py -i interim/antismash/ -o interim/antismash/ -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}
        """