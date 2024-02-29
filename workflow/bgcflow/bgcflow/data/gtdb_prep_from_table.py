import json
import logging
import sys
from pathlib import Path

import pandas as pd

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

metadata_keywords = {
    "metadata_nucleotide": [
        "trna_aa_count",
        "contig_count",
        "n50_contigs",
        "longest_contig",
        "scaffold_count",
        "n50_scaffolds",
        "longest_scaffold",
        "genome_size",
        "gc_percentage",
        "ambiguous_bases",
        "gc_count",
        "l50_contigs",
        "l50_scaffolds",
        "mean_contig_length",
        "mean_scaffold_length",
    ],
    "metadata_gene": [
        "checkm_completeness",
        "checkm_contamination",
        "checkm_strain_heterogeneity",
        "lsu_5s_count",
        "ssu_count",
        "lsu_23s_count",
        "protein_count",
        "coding_density",
        "lsu_23s_contig_len",
        "lsu_23s_length",
        "lsu_23s_query_id",
        "lsu_5s_contig_len",
        "lsu_5s_length",
        "lsu_5s_query_id",
        "lsu_silva_23s_blast_align_len",
        "lsu_silva_23s_blast_bitscore",
        "lsu_silva_23s_blast_evalue",
        "lsu_silva_23s_blast_perc_identity",
        "lsu_silva_23s_blast_subject_id",
        "lsu_silva_23s_taxonomy",
        "checkm_marker_count",
        "checkm_marker_lineage",
        "checkm_marker_set_count",
        "coding_bases",
        "mimag_high_quality",
        "mimag_low_quality",
        "mimag_medium_quality",
        "ssu_contig_len",
        "ssu_gg_blast_align_len",
        "ssu_gg_blast_bitscore",
        "ssu_gg_blast_evalue",
        "ssu_gg_blast_perc_identity",
        "ssu_gg_blast_subject_id",
        "ssu_gg_taxonomy",
        "ssu_length",
        "ssu_query_id",
        "ssu_silva_blast_align_len",
        "ssu_silva_blast_bitscore",
        "ssu_silva_blast_evalue",
        "ssu_silva_blast_perc_identity",
        "ssu_silva_blast_subject_id",
        "ssu_silva_taxonomy",
        "total_gap_length",
        "trna_count",
        "trna_selenocysteine_count",
    ],
    "metadata_ncbi": [
        "ncbi_genbank_assembly_accession",
        "ncbi_strain_identifiers",
        "ncbi_assembly_level",
        "ncbi_assembly_name",
        "ncbi_assembly_type",
        "ncbi_bioproject",
        "ncbi_biosample",
        "ncbi_country",
        "ncbi_date",
        "ncbi_genome_category",
        "ncbi_genome_representation",
        "ncbi_isolate",
        "ncbi_isolation_source",
        "ncbi_lat_lon",
        "ncbi_molecule_count",
        "ncbi_protein_count",
        "ncbi_refseq_category",
        "ncbi_seq_rel_date",
        "ncbi_spanned_gaps",
        "ncbi_species_taxid",
        "ncbi_ssu_count",
        "ncbi_submitter",
        "ncbi_taxid",
        "ncbi_total_gap_length",
        "ncbi_translation_table",
        "ncbi_trna_count",
        "ncbi_unspanned_gaps",
        "ncbi_contig_count",
        "ncbi_contig_n50",
        "ncbi_ncrna_count",
        "ncbi_organism_name",
        "ncbi_rrna_count",
        "ncbi_scaffold_count",
        "ncbi_scaffold_l50",
        "ncbi_scaffold_n50",
        "ncbi_scaffold_n75",
        "ncbi_scaffold_n90",
        "ncbi_total_length",
        "ncbi_ungapped_length",
        "ncbi_wgs_master",
    ],
    "metadata_type_material": [
        "gtdb_type_designation_ncbi_taxa",
        "gtdb_type_designation_ncbi_taxa_sources",
        "gtdb_type_species_of_genus",
    ],
    "metadataTaxonomy": [
        "ncbi_taxonomy",
        "ncbi_taxonomy_unfiltered",
        "gtdb_representative",
        "gtdb_genome_representative",
        "ncbi_type_material_designation",
    ],
}


def empty_template(genome_id):
    empty_template = {
        "genome_id": genome_id,
        "gtdb_url": None,
        "gtdb_release": None,
        "gtdb_taxonomy": {
            "domain": "d__",
            "phylum": "p__",
            "class": "c__",
            "order": "o__",
            "family": "f__",
            "genus": "g__",
            "species": "s__",
        },
    }
    return empty_template


def create_gtdb_metadata_from_table(
    genome_id,
    gtdb_metadata_path,
    gtdb_release,
    outfile,
    metadata_keywords=metadata_keywords,
):
    """
    This function creates GTDB metadata from a given table and writes it to a JSON file.

    Parameters:
    genome_id (str): The genome ID to be processed.
    gtdb_metadata_path (str): The path to the GTDB metadata file.
    gtdb_release (str): The GTDB release version.
    outfile (str): The path to the output file.
    metadata_keywords (dict): The metadata keywords to be processed.

    Returns:
    None
    """
    gtdb_mapping = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    logging.info("Reading GTDB metadata from %s", gtdb_metadata_path)
    df_gtdb_metadata = pd.read_csv(gtdb_metadata_path, sep="\t", low_memory=False)
    df_gtdb_metadata = df_gtdb_metadata.set_index("accession")

    if genome_id.startswith("GCF_"):
        logging.debug("Detected Refseq entry")
        accession = f"RS_{genome_id}"
    elif genome_id.startswith("GCA_"):
        logging.debug("Detected Refseq entry")
        accession = f"GB_{genome_id}"

    if accession in df_gtdb_metadata.index:
        logging.info("Accession %s found in metadata", accession)
        metadata = df_gtdb_metadata.loc[accession, :].to_dict()
        output = populate_from_gtdb_table(
            genome_id, metadata, gtdb_release, metadata_keywords, gtdb_mapping
        )
    else:
        logging.warning(
            "Accession %s not found in metadata. Returning empty json.", accession
        )
        output = empty_template(genome_id)

    logging.info("Metadata creation completed for genome_id %s", genome_id)

    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    logging.debug("Writing metadata to %s", outfile)
    with open(outfile, "w") as file:
        json.dump(output, file, indent=2)
    return


def populate_from_gtdb_table(
    genome_id, metadata, gtdb_release, metadata_keywords, gtdb_mapping
):
    output = {}
    output["genome_id"] = genome_id
    output[
        "gtdb_url"
    ] = f"https://api.gtdb.ecogenomic.org/genome/{genome_id}/taxon-history"
    output["gtdb_release"] = gtdb_release
    output["gtdb_taxonomy"] = {
        gtdb_mapping[v.split("__")[0]]: v for v in metadata["gtdb_taxonomy"].split(";")
    }
    output["metadata_url"] = f"https://api.gtdb.ecogenomic.org/genome/{genome_id}/card"
    output["metadata"] = {}
    output["metadata"]["genome"] = {"accession": genome_id, "name": genome_id}
    for keywords, values in metadata_keywords.items():
        output["metadata"][keywords] = {k: metadata[k] for k in values}
    output["metadata"]["detail"] = "genome_found"
    return output


if __name__ == "__main__":
    create_gtdb_metadata_from_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
