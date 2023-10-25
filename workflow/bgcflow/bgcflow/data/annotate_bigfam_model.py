import logging
import math
import sqlite3
import sys
from pathlib import Path

import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def shannon_div(data):
    """
    Computes the Shannon Diversity Index for the given data.

    Args:
        data (dict): A dictionary containing the counts of different categories.

    Returns:
        float: The computed Shannon Diversity Index.
    """
    total = sum(data.values())
    container = []
    for k, v in data.items():
        p = v / total
        val = p * (math.log(p))
        container.append(val)
    shannon = abs(sum(container))
    return shannon


def get_bigfam_class(gcf_id, conn, threshold="<= 900", rank=0):
    """
    Retrieves the chemical class and chemical subclass for a given gcf_id.

    Args:
        gcf_id (int): The ID of the gene cluster family to retrieve information for.
        conn (sqlite3.Connection): The SQLite connection object to the database.
        threshold (str, optional): The threshold for membership_value, defaults to "<= 900".
        rank (int, optional): The rank of the gcf_membership, defaults to 0.

    Returns:
        dict: A dictionary containing the chemical class, chemical subclass, and related information.
    """
    sql_query = f"""
    SELECT
        "gcf_membership"."gcf_id" AS "gcf_id",
        "gcf_membership"."bgc_id" AS "bgc_id",
        "gcf_membership"."membership_value" AS "membership_value",
        "gcf_membership"."rank" AS "rank",
        "Bgc Class - Bgc"."bgc_id" AS "Bgc Class - Bgc__bgc_id",
        "Bgc Class - Bgc"."chem_subclass_id" AS "Bgc Class - Bgc__chem_subclass_id",
        "Chem Subclass"."id" AS "Chem Subclass__id",
        "Chem Subclass"."class_id" AS "Chem Subclass__class_id",
        "Chem Subclass"."name" AS "Chem Subclass__name",
        "Chem Class - Class"."id" AS "Chem Class - Class__id",
        "Chem Class - Class"."name" AS "Chem Class - Class__name"
    FROM "gcf_membership"
    LEFT JOIN "bgc_class" "Bgc Class - Bgc" ON "gcf_membership"."bgc_id" = "Bgc Class - Bgc"."bgc_id" LEFT JOIN "chem_subclass" "Chem Subclass" ON "Bgc Class - Bgc"."chem_subclass_id" = "Chem Subclass"."id" LEFT JOIN "chem_class" "Chem Class - Class" ON "Chem Subclass"."class_id" = "Chem Class - Class"."id"
    WHERE ("gcf_membership"."gcf_id" = {gcf_id}
       AND "gcf_membership"."membership_value" {threshold} AND "gcf_membership"."rank" = {rank})
    """
    df = pd.read_sql_query(sql_query, conn)
    non_singleton_map = df.bgc_id.value_counts().to_dict()
    bgc = [k for k, v in non_singleton_map.items() if v > 1]
    df.loc[df[~df.bgc_id.isin(bgc)].index, "score"] = 1
    for i in df[df.bgc_id.isin(bgc)].index:
        score = 1 / non_singleton_map[df.loc[i, "bgc_id"]]
        df.loc[i, "score"] = score

    chemical_class = df.groupby("Chem Class - Class__name").sum().score / len(
        df.bgc_id.unique()
    )
    chemical_class = chemical_class.to_dict()
    chemical_subclass = df.groupby("Chem Subclass__name").sum().score / len(
        df.bgc_id.unique()
    )
    chemical_subclass = chemical_subclass.to_dict()
    top_chemical_class = max(chemical_class, key=chemical_class.get)
    top_chemical_subclass = max(chemical_subclass, key=chemical_subclass.get)
    bgc_member = len(df.bgc_id.unique())

    result = {
        "bgc_member": bgc_member,
        "chemical_class_hits": len(df),
        "top_chemical_class": top_chemical_class,
        "top_chemical_class_proportion": chemical_class[top_chemical_class],
        "top_chemical_subclass": top_chemical_subclass,
        "top_chemical_subclass_proportion": chemical_subclass[top_chemical_subclass],
        "chemical_class": chemical_class,
        "chemical_subclass": chemical_subclass,
    }
    return result


def get_bigfam_taxa(gcf_id, conn, threshold="<= 900", rank=0, level=5):
    """
    Retrieves a summary of chemical and taxonomic information for a given gcf_id.

    Args:
        gcf_id (int): The ID of the gene cluster family to retrieve information for.
        conn (sqlite3.Connection): The SQLite connection object to the database.
        threshold (str, optional): The threshold for membership_value, defaults to "<= 900".
        rank (int, optional): The rank of the gcf_membership, defaults to 0.
        level (int, optional): The taxonomic level to retrieve, defaults to 5.

    Returns:
        dict: A dictionary containing a summary of chemical and taxonomic information.
    """
    sql_all_hits = f"""
    SELECT
        "gcf_membership"."gcf_id" AS "gcf_id",
        "gcf_membership"."bgc_id" AS "bgc_id",
        "gcf_membership"."membership_value" AS "membership_value",
        "gcf_membership"."rank" AS "rank"
    FROM "gcf_membership"
    WHERE ("gcf_membership"."gcf_id" = {gcf_id}
        AND "gcf_membership"."membership_value" {threshold} AND "gcf_membership"."rank" = {rank})
    """

    sql_query = f"""
    SELECT
        "gcf_membership"."gcf_id" AS "gcf_id",
        "gcf_membership"."bgc_id" AS "bgc_id",
        "gcf_membership"."membership_value" AS "membership_value",
        "gcf_membership"."rank" AS "rank",
        "Bgc Taxonomy - Bgc"."bgc_id" AS "Bgc Taxonomy - Bgc__bgc_id",
        "Bgc Taxonomy - Bgc"."taxon_id" AS "Bgc Taxonomy - Bgc__taxon_id",
        "Taxon"."id" AS "Taxon__id", "Taxon"."level" AS "Taxon__level",
        "Taxon"."name" AS "Taxon__name", "Taxon Class - Level"."id" AS "Taxon Class - Level__id",
        "Taxon Class - Level"."level" AS "Taxon Class - Level__level",
        "Taxon Class - Level"."name" AS "Taxon Class - Level__name"
    FROM "gcf_membership"
    LEFT JOIN "bgc_taxonomy" "Bgc Taxonomy - Bgc" ON "gcf_membership"."bgc_id" = "Bgc Taxonomy - Bgc"."bgc_id" LEFT JOIN "taxon" "Taxon" ON "Bgc Taxonomy - Bgc"."taxon_id" = "Taxon"."id" LEFT JOIN "taxon_class" "Taxon Class - Level" ON "Taxon"."level" = "Taxon Class - Level"."level"
    WHERE ("gcf_membership"."gcf_id" = {gcf_id}
        AND "gcf_membership"."membership_value" {threshold}
        AND "gcf_membership"."rank" = {rank}
        AND "Taxon"."level" = {level})
    """
    df_all_hits = pd.read_sql_query(sql_all_hits, conn)
    df = pd.read_sql_query(sql_query, conn)
    taxonomic_level = df["Taxon Class - Level__name"].unique()
    if len(taxonomic_level) > 0:
        taxonomic_level = taxonomic_level[0]
    else:
        logging.info(f"No taxonomic level found for {gcf_id}")
        taxonomic_level = "Unassigned"
    df = df_all_hits.merge(
        df, on=["bgc_id", "gcf_id", "membership_value", "rank"], how="outer"
    ).fillna("Unassigned")
    taxa_distribution = df.Taxon__name.value_counts().to_dict()
    top_genus = max(taxa_distribution, key=taxa_distribution.get)
    taxonomic_hits = len(df.bgc_id.unique())
    result = {
        "taxonomic_hits": taxonomic_hits,
        "taxonomic_level": taxonomic_level,
        "H-index": shannon_div(taxa_distribution),
        "richness": len(taxa_distribution.keys()),
        "top_taxa": top_genus,
        "top_taxa_proportion": taxa_distribution[top_genus] / taxonomic_hits,
        "taxa_distribution": taxa_distribution,
    }
    return result


def get_bigfam_summary(gcf_id, conn, threshold="<= 900", rank=0, level=5):
    """
    Retrieves a summary of chemical and taxonomic information for a given gcf_id.

    Args:
        gcf_id (int): The ID of the gene cluster family to retrieve information for.
        conn (sqlite3.Connection): The SQLite connection object to the database.
        threshold (str, optional): The threshold for membership_value, defaults to "<= 900".
        rank (int, optional): The rank of the gcf_membership, defaults to 0.
        level (int, optional): The taxonomic level to retrieve, defaults to 5.

    Returns:
        dict: A dictionary containing a summary of chemical and taxonomic information.
    """
    result = get_bigfam_class(gcf_id, conn)
    result.update(get_bigfam_taxa(gcf_id, conn))
    return {gcf_id: result}


def annotate_bigfam_models(bigfam_model_hits, bigfam_db, output):
    """
    Annotates BigFam models and writes the results to an output file.

    Args:
        bigfam_model_hits (str): The path to the BigFam model hits CSV file.
        bigfam_db (str): The path to the BigFam SQLite database.
        output (str): The path to the output file where the results will be written.
    """
    df_bigfam_model = pd.read_csv(bigfam_model_hits)
    logging.info(f"Reading BiG-FAM database: {bigfam_db}")
    conn = sqlite3.connect(bigfam_db)

    result = {}
    logging.info("Annotating models...")
    with alive_bar(
        len(df_bigfam_model.gcf_id), title="Annotating BiG-FAM models:"
    ) as bar:
        for i in df_bigfam_model.gcf_id:
            value = get_bigfam_summary(i, conn)
            result.update(value)
            bar()

    logging.info(f"Writing results to: {output}")
    df = pd.DataFrame.from_dict(result).T
    df.index.name = "gcf_id"
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output)


if __name__ == "__main__":
    annotate_bigfam_models(sys.argv[1], sys.argv[2], sys.argv[3])
