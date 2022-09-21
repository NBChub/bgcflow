from pathlib import Path
from Bio import SeqIO
import sys
from alive_progress import alive_bar
import json
import logging

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def count_bgcs(gbk_file_path, genome_id=False, outfile=False):
    gbk_file_path = Path(gbk_file_path)
    
    if not type(genome_id) == str:
        genome_id = gbk_file_path.stem
    
    logging.debug(f"Summarize BGCs info for {gbk_file_path}")
    logging.debug(f"Genome id is: {genome_id}")
    records_list = SeqIO.parse(gbk_file_path, 'genbank')

    bgc_stats = {}

    # statistics counter
    bgc_cntr = 0
    protoclusters_cntr = 0
    cand_clusters_cntr = 0
    contig_edge_cntr = 0

    # product capture
    bgc_type_dict = {}

    for rec in records_list:
        # Information on the number of BGCs, protoclusters and candidate clusters
        for feat in rec.features:
            if feat.type == 'region':
                # get region counts
                bgc_cntr = bgc_cntr + 1
                if feat.qualifiers['contig_edge'][0] == 'True':
                    contig_edge_cntr = contig_edge_cntr + 1

                # get product counts
                bgc_type = '.'.join(sorted(feat.qualifiers['product']))
                if bgc_type in bgc_type_dict.keys():
                    bgc_type_dict[bgc_type] = bgc_type_dict[bgc_type] + 1
                else:
                    bgc_type_dict[bgc_type] = 1

                if feat.type == 'protocluster':
                    protoclusters_cntr = protoclusters_cntr + 1
                if feat.type == 'cand_cluster':
                    cand_clusters_cntr = cand_clusters_cntr + 1

                bgc_stats['bgcs_count'] = bgc_cntr
                bgc_stats['bgcs_on_contig_edge'] = contig_edge_cntr
                bgc_stats['protoclusters_count'] = protoclusters_cntr
                bgc_stats['cand_clusters_count'] = cand_clusters_cntr
    result = {genome_id : bgc_stats | {"products" : bgc_type_dict}}
    
    if not type(outfile) == str:
        return result
    else:
        logging.debug(f"Writing output to: {outfile}")
        with open(outfile, "w") as f:
            json.dump(result, f, indent=2)
        return

if __name__ == "__main__":
    count_bgcs(sys.argv[1], sys.argv[2], sys.argv[3])