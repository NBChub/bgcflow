rule databases_setup:
    output:
        "resources/antismash_db.txt",
        "resources/Actinos_6species.gbff"
    conda:
        "../envs/antismash.yaml"
    shell:  
        """
        download-antismash-databases
        antismash --check-prereqs > {output[0]}
        (cd resources && unzip Actinos_6species.zip)
        """

rule automlst_setup:
    output:
        directory("resources/automlst/")
    conda:
        "../envs/automlst.yaml"
    shell:  
        """
        (cd resources && git clone https://bitbucket.org/ziemertlab/automlst.git)
        #wget -c http://automlst.ziemertlab.com/static/refseqreduced.db -P resources/automlst_db
        """