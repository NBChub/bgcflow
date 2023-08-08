import json
import logging
import sys
from pathlib import Path

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.INFO)


def handle_join_locus_location(raw_location):
    """
    Process a locus location string with multiple possible ranges and strands, selecting a consolidated range (for example in PGAP).

    This function takes a locus location string that may contain multiple possible start-stop pairs, each associated
    with a strand. It calculates and returns a single, consolidated locus range that covers all possibilities.

    Args:
        raw_location (str): The locus location string to be processed, possibly starting with 'join'.

    Returns:
        str: A consolidated locus location string representing the encompassing range.

    Example:
        >>> handle_join_locus_location("join{[282863:283224](-), [282406:282864](-)}")
        '[282406:283224](-)'
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


def region_table_builder(f, accession):
    """
    Build a database table entry from a feature of a record representing a genomic region.

    This function takes a feature of a genomic record and extracts relevant information to create an entry for a
    database table representing a genomic region. The function returns a dictionary with the extracted information.

    Args:
        f (dict): A dictionary representing a genomic feature containing information about a region.
        accession (str): The accession number associated with the genomic record.

    Returns:
        dict: A dictionary containing the extracted information structured for database entry.

    Example:
        >>> feature = {...}  # A dictionary representing a genomic feature
        >>> accession = "NC_12345"
        >>> result = region_table_builder(feature, accession)
        >>> print(result)
        {
            'NC_12345.region001': {
                'accession': 'NC_12345',
                'region_number': '001',
                'location': '[100:200]',
                'start_pos': '100',
                'end_pos': '200',
                'contig_edge': 'yes',
                'product': ['Some Product'],
                'rules': ['Some Rule']
            }
        }
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
    """
    Build a database table entry from a feature of a record representing a coding sequence (CDS).

    This function takes a feature of a genomic record representing a coding sequence (CDS) and extracts relevant
    information to create an entry for a database table. The function returns a dictionary with the extracted information.

    Args:
        f (dict): A dictionary representing a genomic feature containing information about a CDS.
        cds_id (str): The unique identifier associated with the CDS.

    Returns:
        dict: A dictionary containing the extracted information structured for database entry.

    Example:
        >>> feature = {...}  # A dictionary representing a genomic feature
        >>> cds_id = "CDS12345"
        >>> result = cdss_table_builder(feature, cds_id)
        >>> print(result)
        {
            'CDS12345': {
                'gene_function': ['Some Gene Function'],
                'locus_tag': 'ABC123',
                'name': 'Gene Name',
                'product': 'Some Product',
                'translation': 'ATGGCCTTAAG...',
                'location': '[100:200](+)',
                'gene_kind': 'Some Gene Kind',
                'codon_start': '1',
                'EC_number': 'EC-12345'
            }
        }
    """
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
    """
    Find regions that overlap with a given coding sequence (CDS) location.

    This function takes a coding sequence (CDS) identifier, a raw location string, and a dictionary of regions,
    and identifies regions from the provided dictionary that overlap with the specified CDS location. It returns
    lists of overlapping region IDs and the corresponding overlapping ranges.

    Args:
        cdss_id (str): The identifier of the coding sequence (CDS) being analyzed.
        location_raw (str): The raw location string of the CDS, including start and stop positions and strand.
        regions_container (dict): A dictionary containing region information with start and stop positions.

    Returns:
        list: A list of region IDs that overlap with the CDS location.
        list: A list of overlapping ranges between the CDS and each overlapping region.

    Example:
        >>> cdss_id = "CDS12345"
        >>> location_raw = "[100:200](+)"
        >>> regions = {
        ...     'Region001': {'start_pos': 150, 'end_pos': 250},
        ...     'Region002': {'start_pos': 50, 'end_pos': 120},
        ...     ...
        ... }
        >>> overlaps, overlapping_ranges = region_finder(cdss_id, location_raw, regions)
        >>> print(overlaps)
        ['Region001', 'Region002']
        >>> print(overlapping_ranges)
        [range(150, 201), range(100, 121)]
    """
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
    hits_value = []

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
                hits_value.append(value)
            # else:
            #    hits.append(np.nan)
        except IndexError:
            pass
    return hits, hits_value


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
    if "accessions" in record["annotations"].keys():
        if len(record["annotations"]["accessions"]) != 1:
            logging.warning(
                f'More than one accession in record: {record["annotations"]["accessions"]}'
            )
            logging.debug(
                f'Grabbing only the first accession: {record["annotations"]["accessions"][0]}'
            )
        record_container["accessions"] = record["annotations"]["accessions"][0]
    else:
        logging.warning(
            f'record["annotations"] does not have "accessions" information.'
        )
        logging.debug(f'Using sequence_id: {sequence_id} as "accessions"')
        record_container["accessions"] = sequence_id
    record_container["genome_id"] = genome_id
    return sequence_id, record_container


def get_region_information(record, genome_id, r, table_regions, n_hits=1):
    """
    Retrieve DNA sequences and associated information from a sequence record.

    This function takes a sequence record and a genome identifier, and extracts relevant DNA sequence data and metadata
    from the record. It returns a tuple containing the sequence identifier and a dictionary with extracted information.

    Args:
        record (dict): A dictionary representing a sequence record, typically obtained from a sequence database.
        genome_id (str): The identifier of the genome to which the sequence belongs.

    Returns:
        tuple: A tuple containing the sequence identifier and a dictionary with extracted sequence and metadata.

    Example:
        >>> sequence_record = {...}  # A dictionary representing a sequence record
        >>> genome_id = "Genome123"
        >>> sequence_id, sequence_info = get_dna_sequences(sequence_record, genome_id)
        >>> print(sequence_id)
        'Sequence789'
        >>> print(sequence_info)
        {
            'seq': 'ATCGATCG...',
            'description': 'Some sequence description',
            'molecule_type': 'DNA',
            'topology': 'linear',
            'accessions': 'Accession123',
            'genome_id': 'Genome123'
        }
    """
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
    """
    Extract CDS (coding sequence) information from a sequence record and update a database table.

    This function takes a sequence record, a genome identifier, and existing database tables for regions and CDSs.
    It processes CDS features in the sequence record, extracts relevant information, and updates the CDS database table.
    It also assigns overlapping regions to each CDS based on their genomic locations.

    Args:
        record (dict): A dictionary representing a sequence record.
        genome_id (str): The identifier of the genome to which the sequence belongs.
        table_regions (dict): A dictionary containing region information with start and stop positions.
        table_cdss (dict): A dictionary representing the current CDS database table.
        accession (str): The accession number associated with the genomic record.

    Returns:
        None

    Example:
        >>> sequence_record = {...}  # A dictionary representing a sequence record
        >>> genome_id = "Genome123"
        >>> regions_table = {...}  # A dictionary representing a database table for regions
        >>> cdss_table = {...}  # A dictionary representing a database table for CDSs
        >>> accession = "NC_12345"
        >>> get_cdss_information(sequence_record, genome_id, regions_table, cdss_table, accession)
        # The cdss_table dictionary is updated with CDS information
    """
    starting_size = len(table_cdss)
    cds_ctr = len(table_cdss)
    logging.debug(f"Starting CDS counter from {cds_ctr}")
    cds_feat = [i for i in record["features"] if i["type"] == "CDS"]
    logging.info(f"Grabbing feature information from {len(cds_feat)} CDS")
    for feature in cds_feat:
        cds_id = f"{accession}-cds_{cds_ctr}"
        assert (
            cds_id not in table_cdss.keys()
        ), f"CDS ID already exists! {cds_id}: {table_cdss[cds_id]}"
        cdss_data = cdss_table_builder(feature, cds_id)
        cdss_data[cds_id]["accessions"] = accession
        cdss_data[cds_id]["genome_id"] = genome_id
        region_hits, hits_value = region_finder(
            cds_id, cdss_data[cds_id]["location"], table_regions
        )
        if len(region_hits) > 1:
            overlaps_dict = {i: hits_value[num] for num, i, in enumerate(region_hits)}
            logging.warning(
                f"Unable to decide regions for {cdss_data[cds_id]['locus_tag']}:{cdss_data[cds_id]['location']} | {overlaps_dict}"
            )
            region_hits = max(overlaps_dict, key=lambda k: len(overlaps_dict[k]))
            logging.info(
                f"Picking {region_hits} as it has bigger overlap size: {len(overlaps_dict[region_hits])} bp"
            )
            cdss_data[cds_id]["region_id"] = region_hits
        elif len(region_hits) == 1:
            cdss_data[cds_id]["region_id"] = region_hits[0]
        else:
            cdss_data[cds_id]["region_id"] = None
        cdss_data[cds_id]["genome_id"] = genome_id
        cdss_data = {v["locus_tag"]: v for k, v in cdss_data.items()}
        table_cdss.update(cdss_data)
        cds_ctr = cds_ctr + 1
    logging.info(
        f"Original table size: {starting_size}, cdss_data: {len(cds_feat)}, size after addition: {len(table_cdss)}"
    )
    return


def antismash_json_exporter(json_path, output_dir, genome_id=False, n_hits=1):
    """
    Read AntiSMASH JSON output and extract region, DNA sequence, and CDS information.

    This function reads an AntiSMASH JSON output file, processes its contents, and extracts relevant information
    including DNA sequences, regions, and coding sequences (CDSs). It then generates separate JSON files for each
    extracted information type and saves them to the specified output directory.

    Args:
        json_path (str): Path to the AntiSMASH JSON output file.
        output_dir (str): Directory where the extracted information JSON files will be saved.
        genome_id (str, optional): Identifier for the genome associated with the processed data.
            If not provided, the function uses the input filename without the '.gbk' extension.
        n_hits (int, optional): The number of hits to consider when processing region information.
            Default is 1.

    Returns:
        None

    Example:
        >>> json_path = 'antismash_output.json'
        >>> output_dir = 'output_data'
        >>> genome_id = 'Genome123'
        >>> n_hits = 2
        >>> antismash_json_exporter(json_path, output_dir, genome_id, n_hits)
        # JSON files for DNA sequences, regions, and CDSs are generated and saved in the output directory.
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
        region_map = {
            k: v for k, v in table_regions.items() if v["accession"] == record["id"]
        }
        get_cdss_information(record, genome_id, region_map, table_cdss, accession)
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
