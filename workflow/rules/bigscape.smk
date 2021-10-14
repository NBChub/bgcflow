rule install_bigscape:
    output:
        directory("resources/BiG-SCAPE/")
    conda:
        "../envs/bigscape.yaml"
    shell:
        """
        (cd resources && wget https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/archive/master/BiG-SCAPE-master.zip)
        (cd resources && unzip BiG-SCAPE-master.zip && mv BiG-SCAPE-master/ BiG-SCAPE/ && rm BiG-SCAPE-master.zip)
        (cd resources/BiG-SCAPE && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz)
        (cd resources/BiG-SCAPE && hmmpress Pfam-A.hmm)
        """

rule bigscape:
    input: 
        "resources/BiG-SCAPE",
        expand("data/interim/antismash/{strains}", strains=STRAINS)
    output:
        "data/interim/bigscape/index.html",
    conda:
        "../envs/bigscape.yaml"
    params:
        label = "all",
    log:
        "workflow/report/bigscape.log"
    threads: 16
    shell:
        """
        python resources/BiG-SCAPE/bigscape.py -i data/interim/antismash/ -o data/interim/bigscape/ -c {threads} --cutoff 0.3 0.4 0.5 --include_singletons --label {params.label} --hybrids-off --mibig --verbose > {log}
        """