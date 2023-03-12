import json
import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def region_table_builder(f, accession):
    """
    Given a feature of a record, return value to build database table
    """
    # grab region values
    region_number = f["qualifiers"]["region_number"][0]
    region_id = f"{accession}.region{str(region_number).zfill(3)}"
    location = f["location"].strip("[").strip("]")
    start_pos, end_pos = location.split(":")
    contig_edge = f["qualifiers"]["contig_edge"][0]
    # fill values
    value = {
        region_id: {
            "accession": accession,
            "region_number": region_number,
            "location": location,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "contig_edge": contig_edge,
            "product": f["qualifiers"]["product"],
            "rules": f["qualifiers"]["rules"],
        }
    }
    return value


def get_antismash_overview(json_path, outfile, genome_id=False, n_hits=1):
    """
    Read antismash json output and get the summary page into a json format
    """
    path = Path(json_path)
    with open(path, "r") as f:
        data = json.load(f)

    if not genome_id:
        genome_id = data["input_file"].strip(".gbk")
    else:
        pass

    # iterating over record
    output = {}
    for r, record in enumerate(data["records"]):
        logging.info(f"Getting antismash regions from record: {record['id']}")
        region_db = {}

        region_feat = [i for i in record["features"] if i["type"] == "region"]
        for f in region_feat:
            region_db.update(region_table_builder(f, record["id"]))

        for c, area in enumerate(record["areas"]):
            cluster_id = f"{r+1}.{c+1}"
            output_cluster = {}
            logging.info(f"Grabbing information from region {cluster_id}")
            # _from = area['start']
            # _to = area['end']
            # products = area['products']

            knownclusterblast = record["modules"]["antismash.modules.clusterblast"][
                "knowncluster"
            ]["results"][c]

            assert n_hits > 0

            output_hits = []

            for n, hits in enumerate(knownclusterblast["ranking"]):
                if n + 1 <= (n_hits):
                    most_similar_mibig_id = hits[0]["accession"]
                    most_similar_mibig_description = hits[0]["description"]
                    most_similar_mibig_clustertype = hits[0]["cluster_type"]
                    n_genes_in_target = len(hits[0]["tags"])
                    n_genes_hits = hits[1]["hits"]
                    hit_similarity = n_genes_hits / n_genes_in_target
                    output_hits.append(
                        {
                            "most_similar_known_cluster_id": most_similar_mibig_id,
                            "most_similar_known_cluster_description": most_similar_mibig_description,
                            "most_similar_known_cluster_type": most_similar_mibig_clustertype,
                            "similarity": hit_similarity,
                        }
                    )
                else:
                    pass

            bgc_id = f"{record['id']}.region{str(c+1).zfill(3)}"
            output_cluster = {
                "genome_id": data["input_file"].strip(".gbk"),
                "region": cluster_id,
            }

            for column in [
                "accession",
                "start_pos",
                "end_pos",
                "contig_edge",
                "product",
            ]:
                output_cluster[column] = region_db[bgc_id][column]
            try:
                output_cluster["region_length"] = int(output_cluster["end_pos"]) - int(
                    output_cluster["start_pos"]
                )
            except ValueError:
                logging.warning(
                    f'Error calculating region length. Region might be incomplete: {output_cluster["start_pos"]}:{output_cluster["end_pos"]}'
                )
                start_pos = "".join(
                    [s for s in output_cluster["start_pos"] if s.isdigit()]
                )
                logging.warning(
                    f'Correcting start position from {output_cluster["start_pos"]} to {start_pos}'
                )
                output_cluster["start_pos"] = start_pos
                output_cluster["region_length"] = int(output_cluster["end_pos"]) - int(
                    output_cluster["start_pos"]
                )

            if len(output_hits) == 1:
                for k in output_hits[0].keys():
                    output_cluster[k] = output_hits[0][k]

            output[bgc_id] = output_cluster

    with open(outfile, "w") as f:
        logging.info(f"Writing results to {outfile}")
        json.dump(output, f, indent=2)

    return output


if __name__ == "__main__":
    get_antismash_overview(sys.argv[1], sys.argv[2], sys.argv[3])
