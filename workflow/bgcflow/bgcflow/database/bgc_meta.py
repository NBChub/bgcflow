import json
import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.INFO)


def handle_join_locus_location(raw_location):
    """
    Function to handle locus location with two or more possibilities (PGAP?).
    Will go through all possible start and stop location, then select the range that covers all.
    Example:
        handle_join_locus_location("join{[282863:283224](-), [282406:282864](-)}")
    Returns a string:
        "[282406:283224](-)"
    """
    container = {"+": {"start": [], "stop": []}, "-": {"start": [], "stop": []}}

    if raw_location.startswith("join"):
        raw_location = raw_location.strip("join")

    q = [
        ("".join(c for c in s if c not in ")[]{}<>")).split("(")
        for s in raw_location.split(", ")
    ]
    for item in q:
        start, stop = item[0].split(":")
        strand = item[1]
        if strand == "+":
            container["+"]["start"].append(int(start))
            container["+"]["stop"].append(int(stop))
        elif strand == "-":
            container["-"]["start"].append(int(start))
            container["-"]["stop"].append(int(stop))

    for k, v in container.items():
        if (len(v["start"]) == 0) or (len(v["stop"]) == 0):
            pass
        else:
            location = f"[{min(v['start'])}:{max(v['stop'])}]({k})"
    return location


# regions
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


def cdss_table_builder(f, cds_id):
    # location, strand = f['location'].strip("[").strip(")").split("](")
    try:
        gene_function = f["qualifiers"]["gene_functions"]
    except KeyError:
        logging.debug(
            f"{cds_id} does not have gene_function. Available values:\n{f['qualifiers'].keys()}"
        )
        gene_function = None

    try:
        name = f["qualifiers"]["gene"][0]
    except KeyError:
        logging.debug(
            f"{cds_id} does not have gene name. Available values:\n{f['qualifiers'].keys()}"
        )
        name = None

    try:
        gene_kind = f["qualifiers"]["gene_kind"][0]
    except KeyError:
        logging.debug(
            f"{cds_id} does not have gene kind. Available values:\n{f['qualifiers'].keys()}"
        )
        gene_kind = None

    try:
        EC_number = f["qualifiers"]["EC_number"][0]
    except KeyError:
        logging.debug(
            f"{cds_id} does not have EC_number. Available values:\n{f['qualifiers'].keys()}"
        )
        EC_number = None

    value = {
        cds_id: {
            "gene_function": gene_function,
            "locus_tag": f["qualifiers"]["locus_tag"][0],
            "name": name,
            "product": f["qualifiers"]["product"][0],
            "translation": f["qualifiers"]["translation"][0],
            "location": f["location"],
            "gene_kind": gene_kind,
            "codon_start": f["qualifiers"]["codon_start"][0],
            "EC_number": EC_number,
        }
    }
    return value


def region_finder(cdss_id, location_raw, regions_container):
    try:
        location, strand = location_raw.strip("[").strip(")").split("](")
    except ValueError:
        logging.warning(f"Unusual location format: {location_raw}")
        location_raw = handle_join_locus_location(location_raw)
        location, strand = location_raw.strip("[").strip(")").split("](")

    q_start, q_stop = [i for i in location.split(":")]
    try:
        q_start, q_stop = int(q_start), int(q_stop)
    except ValueError:
        logging.warning(
            f"Unusual location format: {cdss_id} {q_start}:{q_stop} ({strand})"
        )
        q_start = int("".join([s for s in q_start if s.isdigit()]))
        q_stop = int("".join([s for s in q_stop if s.isdigit()]))
        logging.info(f"Keeping only digits: {cdss_id} {q_start}:{q_stop} ({strand})")

    query = range(q_start, q_stop)

    hits = []

    for region_id in regions_container.keys():
        start_pos = int(
            "".join(
                [s for s in regions_container[region_id]["start_pos"] if s.isdigit()]
            )
        )
        end_pos = int(
            "".join([s for s in regions_container[region_id]["end_pos"] if s.isdigit()])
        )
        target = range(start_pos, end_pos)
        try:
            value = range(max(query[0], target[0]), min(query[-1], target[-1]) + 1)
            if len(value) > 1:
                logging.debug(f"{cdss_id} overlaps with {region_id} at {value}")
                hits.append(region_id)
            # else:
            #    hits.append(np.nan)
        except IndexError:
            pass
    return hits


def get_dna_sequences(record, genome_id):
    """
    Given a sequence record, return DNA sequences and information
    """
    sequence_id = record["id"]
    logging.info(f"Getting dna_sequences information for {sequence_id}")
    record_container = {}
    record_container["seq"] = record["seq"]["data"]
    record_container["description"] = record["description"]
    record_container["molecule_type"] = record["annotations"]["molecule_type"]
    record_container["topology"] = record["annotations"]["topology"]
    if len(record["annotations"]["accessions"]) != 1:
        logging.warning(
            f'More than one accession in record: {record["annotations"]["accessions"]}'
        )
        logging.debug(
            f'Grabbing only the first accession: {record["annotations"]["accessions"][0]}'
        )
    record_container["accessions"] = record["annotations"]["accessions"][0]
    record_container["genome_id"] = genome_id
    return sequence_id, record_container


def get_region_information(record, genome_id, r, table_regions, n_hits=1):
    region_db = {}
    region_feat = [i for i in record["features"] if i["type"] == "region"]
    for f in region_feat:
        region_db.update(region_table_builder(f, record["id"]))
    for c, area in enumerate(record["areas"]):
        cluster_id = f"{r+1}.{c+1}"
        output_cluster = {}
        logging.info(f"Grabbing information from region {cluster_id}")
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
            "genome_id": genome_id,
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

        q_start, q_stop = output_cluster["start_pos"], output_cluster["end_pos"]
        try:
            q_start, q_stop = int(q_start), int(q_stop)
        except ValueError:
            logging.warning(f"Unusual location format: {q_start}:{q_stop}")
            q_start = int("".join([s for s in q_start if s.isdigit()]))
            q_stop = int("".join([s for s in q_stop if s.isdigit()]))
            logging.info(f"Keeping only digits: {q_start}:{q_stop}")
        output_cluster["region_length"] = q_stop - q_start

        if len(output_hits) == 1:
            for k in output_hits[0].keys():
                output_cluster[k] = output_hits[0][k]

        table_regions[bgc_id] = output_cluster

    return


def get_cdss_information(record, genome_id, table_regions, table_cdss, accession):
    # cds_db = {}
    cds_ctr = 1
    cds_feat = [i for i in record["features"] if i["type"] == "CDS"]
    logging.info(f"Grabbing feature information from {len(cds_feat)} CDS")
    for feature in cds_feat:
        cds_id = f"{accession}-cds_{cds_ctr}"
        cdss_data = cdss_table_builder(feature, cds_id)
        cdss_data[cds_id]["accessions"] = accession
        cdss_data[cds_id]["genome_id"] = genome_id
        region_hits = region_finder(
            cds_id, cdss_data[cds_id]["location"], table_regions
        )
        if len(region_hits) > 0:
            cdss_data[cds_id]["region_id"] = region_hits[0]
        else:
            cdss_data[cds_id]["region_id"] = None
        cdss_data[cds_id]["genome_id"] = genome_id
        table_cdss.update(cdss_data)
        cds_ctr = cds_ctr + 1
    return


def antismash_json_exporter(json_path, output_dir, genome_id=False, n_hits=1):
    """
    Read antismash json output and get region, dna_sequences, and cdss information
    """
    # Create output containers
    table_dna_sequences = {}
    table_regions = {}
    table_cdss = {}

    # load json file
    path = Path(json_path)
    with open(path, "r") as f:
        data = json.load(f)

    if not genome_id:
        genome_id = data["input_file"].strip(".gbk")
    else:
        pass

    # Iterating over record
    logging.info(f"Extracting information from {data['input_file']}")
    for r, record in enumerate(data["records"]):
        # Grab DNA information
        sequence_id, record_container = get_dna_sequences(record, genome_id)
        table_dna_sequences[sequence_id] = record_container
        # ---------------------------------------------------

        # Grab region information
        logging.info(f"Extracting regions and cds information from {sequence_id}")
        get_region_information(record, genome_id, r, table_regions)
        # ---------------------------------------------------

        # Grab cdss information
        accession = record_container["accessions"]
        get_cdss_information(record, genome_id, table_regions, table_cdss, accession)
        # ---------------------------------------------------

    # write json files
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    target_jsons = {
        "dna_sequences": table_dna_sequences,
        "regions": table_regions,
        "cdss": table_cdss,
    }
    for k in target_jsons.keys():
        outfile_name = output_dir / f"{genome_id}_{k}.json"
        with open(outfile_name, "w") as output_file:
            logging.info(f"Writing {k} output to {outfile_name}")
            json.dump(target_jsons[k], output_file, indent=2)
    return


if __name__ == "__main__":
    antismash_json_exporter(sys.argv[1], sys.argv[2], sys.argv[3])
