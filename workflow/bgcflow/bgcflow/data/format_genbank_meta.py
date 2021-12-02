import subprocess
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


def add_bgcflow_comments(gbk_in_path, gbk_out_path, version):
    """ 
    Add bgcflow meta-annotation to genbank output
    """
    gbk_in_path = Path(gbk_in_path)
    gbk_out_path = Path(gbk_out_path)

    version = get_version(version)
    records = SeqIO.parse(gbk_in_path, 'genbank')

    antismash_comment = (
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
        comment = antismash_comment
        if 'comment' in record.annotations:
            record.annotations['comment'] += '\n' + comment
        else:
            record.annotations['comment'] = comment
        new_records.append(record)
    
    with open(gbk_out_path, "w") as output_handle:
        SeqIO.write(new_records, output_handle, "genbank")

add_bgcflow_comments(snakemake.input.gbk_prokka, snakemake.output.gbk_processed, snakemake.params.version)