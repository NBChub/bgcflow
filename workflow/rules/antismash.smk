"""
author: matinnu
date: 2021-06-12
Run:

"""
from pathlib import Path
import pandas as pd

TARGET = Path('/data/a/matinnu/data')
df_target = pd.read_csv(TARGET / 'metadata.txt', sep='\t')
STRAINS = df_target.strain.to_list()

rule all:
    input:
        expand(TARGET / "genomes/antiSMASH/{strains}_antiSMASH", strains = STRAINS),
        Path('../databases')

rule antismash_db:
    input: 
        TARGET / 'metadata.txt'
    output:
        directory(Path('../databases'))
    conda:
        "envs/antismash.yaml"
    threads: 12
    shell:
        """
        download-antismash-databases --database-dir {output}
        """

rule antismash:
    input: 
        TARGET / "genomes/{strains}_prokka_actinoannotPFAM",
        Path('../databases')
    output:
        directory(TARGET / "genomes/antiSMASH/{strains}_antiSMASH")
    conda:
        "envs/antismash.yaml"
    threads: 12
    shell:
        """
        antismash --genefinding-tool prodigal --output-dir {output} --cb-general --cb-subclusters --cb-knownclusters -c {threads} {input}/{wildcards.strains}.gbk
        """