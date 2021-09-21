rule autoMLST:
    input:
        "results/genomes/{strains}/final_polish.fasta",
        "resources/automlst/"   
    output:
        directory("results/genomes/{strains}/automlst"),        
    conda: "../envs/automlst.yaml"
    threads: 12
    shell:
        """
        mkdir {output[0]}
        cp {input[0]} {output[0]}/.
        python resources/automlst/automlst.py -ref resources/automlst_db/refseqreduced.db -rd resources/automlst_db --cpu {threads} {output[0]} {output[0]}
        """