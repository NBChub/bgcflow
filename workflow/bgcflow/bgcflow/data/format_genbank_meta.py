import subprocess
import sys
import pandas as pd
from Bio import SeqIO
from datetime import datetime

from pathlib import Path

def get_git_version():
    """
    Get the sha1 of the current git version
    """
    
    git_version = ""
    try:
        version_cmd = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],universal_newlines=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        status_cmd = subprocess.run(['git', 'status', '--porcelain'],universal_newlines=True,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        git_version = str(version_cmd.stdout.strip())
        changes = str(status_cmd.stdout).strip()
        if changes:
            git_version += "(changed)"
    except OSError:
        pass
    
    return git_version


def get_version(version):
    """
    Get the current version string
    """
    
    git_version = get_git_version()
    if git_version:
        version += "-%s" % git_version
    return version


def add_bgcflow_comments(gbk_in_path, version, json_path, genome_id, gbk_out_path):
    """ 
    Add bgcflow meta-annotation to genbank output
    """

    version = get_version(version)
    records = SeqIO.parse(gbk_in_path, 'genbank')

    df_gtdb_meta = pd.read_json(json_path, orient="index").T
    df_gtdb_meta.fillna('Unclassified', inplace=True)
    df_tax = pd.DataFrame([i for i in df_gtdb_meta.gtdb_taxonomy])
    for col in df_tax.columns:
        df_tax[col] = df_tax[col].apply(lambda x: x.split('__')[1])
    df_tax.columns = [c.title() for c in df_tax.columns]
    tax_levels = ['Domain','Phylum','Class','Order','Family','Genus','Species']
    taxonomy_str = df_Tax.loc[0, tax_levels].tolist()
    print('Taxonomy found:', taxonomy_str)

    bgcflow_comment = (
        "##BGCflow-Data-START##\n"
        "Version      :: {version}\n"
        "Run date     :: {date}\n"
        "##BGCflow-Data-END##"
        ).format(
            version=version,
            date=str(datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        )

    new_records = []
    for record in records:
        comment = bgcflow_comment
        if 'comment' in record.annotations:
            record.annotations['comment'] += '\n' + comment
        else:
            record.annotations['comment'] = comment

        if 'organism' in record.annotations:
            organism = record.annotations['organism']
            if 'Unclassified' in organism:
                record.annotations['organism'] = organism.split(' Unclassified')[0].strip()
        
        record.annotations['taxonomy'] = taxonomy_str

        new_records.append(record)
    
    with open(gbk_out_path, "w") as output_handle:
        SeqIO.write(new_records, output_handle, "genbank")


if __name__ == "__main__":
    add_bgcflow_comments(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    