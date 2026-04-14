import logging
import os
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import networkx as nx
import pandas as pd
from alive_progress import alive_bar

log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)


def get_cluster_dataframes(
    df_genomes, df_nodes, mibig_bgc_table, as_dir="../data/antismash/"
):
    """
    Returns two dataframes of clusters with information from genomes and MIBIG
    """

    df_clusters = pd.DataFrame(
        columns=["product", "bigscape_class", "genome_id", "accn_id"]
    )
    df_known = pd.DataFrame(
        columns=[
            "product",
            "bigscape_class",
            "biosyn_class",
            "compounds",
            "chem_acts",
            "accession",
            "completeness",
            "organism_name",
            "ncbi_tax_id",
            "publications",
            "evidence",
        ]
    )
    df_mibig_bgcs = pd.read_csv(mibig_bgc_table, index_col="mibig_id")

    # Generate bgcs dataframe with metadata from df_nodes and df_genomes
    logging.info(
        "Generating bgcs dataframe with metadata from df_nodes and df_genomes..."
    )
    with alive_bar(len(df_genomes.index)) as bar:
        for genome_id in df_genomes.index:
            logging.info(f"Processing BGCs in the genome: {genome_id}")
            if genome_id in os.listdir(as_dir):
                genome_dir = os.path.join(as_dir, genome_id)
                bgc_id_list = []
                for region in os.listdir(genome_dir):
                    if "region" in region and region.endswith(".gbk"):
                        region_number = int(region.split("region")[1].split(".")[0])
                        bgc_id = region + "_region_" + str(region_number)
                        bgc_id_list.append(bgc_id)
                for bgc_id in bgc_id_list:
                    logging.debug(f"Processing {bgc_id}")
                    if bgc_id in df_nodes.index:
                        df_clusters.loc[bgc_id, "genome_id"] = genome_id
                        df_clusters.loc[bgc_id, "product"] = df_nodes.loc[
                            bgc_id, "Category"
                        ]
                        df_clusters.loc[bgc_id, "bigscape_class"] = df_nodes.loc[
                            bgc_id, "Class"
                        ]
                        df_clusters.loc[bgc_id, "accn_id"] = df_nodes.loc[
                            bgc_id, "GBK"
                        ]
                    else:
                        logging.debug(f"{bgc_id} not in df_nodes")
            else:
                logging.warning(f"{genome_id} not in directory!")
            bar()

    # Generate separate table for known BGCs from MIBIG
    logging.info("Generating separate table for known BGCs from MIBIG")
    with alive_bar(len(df_nodes.index)) as bar:
        for bgc_id in df_nodes.index:
            if "BGC" in bgc_id:
                df_known.loc[bgc_id, "product"] = df_nodes.loc[
                    bgc_id, "Category"
                ]
                df_known.loc[bgc_id, "bigscape_class"] = df_nodes.loc[
                    bgc_id, "Class"
                ]
                if "." in bgc_id:
                    mibig_id = bgc_id.split(".")[0]
                else:
                    mibig_id = bgc_id

                if mibig_id in df_mibig_bgcs.index:
                    for col in df_mibig_bgcs.columns:
                        df_known.loc[bgc_id, col] = df_mibig_bgcs.loc[mibig_id, col]
                else:
                    logging.debug(f"{bgc_id} not in df_mibig_bgcs dataframe index")
            bar()

    return df_known, df_clusters


def add_bigscape_families(df_clusters, df_known, net_data_path, cutoff):
    """
    Adds GCC and GCF numbers detected by BiG-SCAPE clustering for different cut-offs
    """

    cluster_class_set = [
        cluster_class
        for cluster_class in os.listdir(net_data_path)
        if cluster_class != "record_annotations.tsv" and (net_data_path / cluster_class).is_dir()
    ]

    for cluster_class in cluster_class_set:
        logging.info(f"Processing all BGCs from {cluster_class}")
        class_dir = net_data_path / cluster_class
        clust_path = class_dir / (cluster_class + "_clustering_c" + cutoff + ".tsv")
        if clust_path.exists():
            df_clusters = read_clustering_df(df_clusters, clust_path, cutoff)
            df_known = read_clustering_df(df_known, clust_path, cutoff)

    return df_clusters, df_known


def read_clustering_df(df_clusters, clust_path, cutoff):
    """
    Adds GCF (Gene Cluster Family) and GCC (Clan) numbers for each BGC
    """

    df_clust = pd.read_csv(clust_path, sep="\t", index_col="Record", dtype=str)
    col_name = "gcf_" + cutoff

    for bgc in df_clusters.index:
        if bgc in df_clust.index:
            df_clusters.loc[bgc, col_name] = df_clust.loc[bgc, "Family"]
            df_clusters.loc[bgc, "Clan Number"] = df_clust.loc[bgc, "CC"]

    return df_clusters


def get_bigscape_network(net_data_path, cutoff="0.3"):
    """
    Reads similarity network for a particular cutoff
    """

    cluster_class_set = [
        cluster_class
        for cluster_class in os.listdir(net_data_path)
        if cluster_class != "record_annotations.tsv" and (net_data_path / cluster_class).is_dir()
    ]

    df_network = pd.DataFrame()
    for cluster_class in cluster_class_set:
        class_dir = net_data_path / cluster_class
        net_file = cluster_class + "_c" + cutoff + ".network"
        net_path = class_dir / net_file
        if net_path.exists():
            df_class_net = pd.read_csv(net_path, sep="\t")
            df_network = pd.concat([df_network, df_class_net], ignore_index=True)

    return df_network


def get_network_graph(df_network, weight="distance"):
    """
    Returns networkX graph for a given network
    """

    G_clusters = nx.from_pandas_edgelist(
        df_network, "Record_a", "Record_b", weight
    )
    G_families = nx.connected_components(G_clusters)
    family_nodes = [c for c in sorted(G_families, key=len, reverse=True)]

    return G_clusters, family_nodes


def remove_single_mibig(df_network, df_known, family_nodes):
    """
    Removes singleton MIBIG BGCs from the network
    """
    logging.info("Removing singleton MIBIG BGCs from the network")
    nodes_to_remove = []

    for fam in family_nodes:
        single_mibig = True

        for node in fam:
            if "BGC" not in node:
                single_mibig = False

        if single_mibig:
            for node in fam:
                if node not in nodes_to_remove:
                    nodes_to_remove.append(node)
    logging.debug(
        f"{len(nodes_to_remove)} number of MIBIG nodes will be removed from analysis due to no similarity"
    )
    df_network = df_network[~df_network["Record_a"].isin(nodes_to_remove)]
    df_network = df_network[~df_network["Record_b"].isin(nodes_to_remove)]
    for node in nodes_to_remove:
        if node in df_known.index:
            df_known = df_known.drop(node)

    return df_network, df_known


def get_family_graph(G_clusters):
    """
    Returns families list as networkx graph
    """

    # Find connected components or cluster families
    Families_list = list(nx.connected_components(G_clusters))
    # Sort families in decreasing order of occurrence
    family_size = [len(family) for family in Families_list]
    orig_index = list(range(len(family_size)))
    orig_index_fam_size_dict = dict(zip(orig_index, family_size))

    sorted_index_fam_size_dict = OrderedDict(
        sorted(orig_index_fam_size_dict.items(), key=lambda x: x[1])
    )
    new_index = list(range(len(family_size) - 1, -1, -1))
    orig_index_sorted = list(sorted_index_fam_size_dict.keys())
    new_orig_index_dict = dict(zip(new_index, orig_index_sorted))

    # Ordered family graphs
    family_graphs = [
        Families_list[new_orig_index_dict[fam_id]]
        for fam_id in range(len(Families_list))
    ]

    return family_graphs


def update_cluster_family(
    df_clusters, df_known, family_nodes, mibig_bgc_table, cutoff="0.30"
):
    """
    Updates df_clusters with family ids (connected components)
    """
    df_mibig_bgcs = pd.read_csv(mibig_bgc_table, index_col="mibig_id")

    df_families = pd.DataFrame(
        columns=["fam_type", "fam_name", "clusters_in_fam", "mibig_ids"]
    )
    for cntr in range(len(family_nodes)):
        fam_id = cntr + 1
        family = family_nodes[cntr]
        known_bgcs = [bgc for bgc in family if bgc in df_mibig_bgcs.index]  #

        if len(known_bgcs) > 0:
            df_families.loc[fam_id, "fam_type"] = "known_family"
            logging.debug(
                f'Family {fam_id} has {len(known_bgcs)} known BGCs: {", ".join(known_bgcs)}'
            )

            # handle missing entry from MIBIG release
            try:
                known_compounds = ";".join(
                    df_known.loc[known_bgcs, "compounds"].tolist()
                )
            except TypeError:
                logging.warning(
                    f"MIBIG entries {known_bgcs} returned no compounds: {df_known.loc[known_bgcs, 'compounds'].values}"
                )
                logging.warning(
                    "This is likely because of missing JSON entry in the downloaded MIBIG release."
                )
                logging.warning(
                    "Check if this is the case in the resources/mibig/df_mibig_bgcs.csv folder"
                )
                logging.warning("Replacing compound value with link to MIBIG...")
                for h in known_bgcs:
                    link_mibig = f"https://mibig.secondarymetabolites.org/repository/{h.split('.')[0]}/index.html"
                    logging.warning(
                        f"{h}: Replacing compound {df_known.loc[h, 'compounds']} with {link_mibig}"
                    )
                    df_known.loc[h, "compounds"] = link_mibig

            df_families.loc[fam_id, "fam_name"] = known_compounds
            df_families.loc[fam_id, "clusters_in_fam"] = len(family)
            df_families.loc[fam_id, "mibig_ids"] = ";".join(known_bgcs)
        else:
            df_families.loc[fam_id, "fam_type"] = "unknown_family"
            user_bgcs_in_family = [bgc for bgc in family if bgc in df_clusters.index]
            bgc_class = ",".join(
                df_clusters.loc[user_bgcs_in_family, "bigscape_class"].unique().tolist()
            )
            df_families.loc[fam_id, "fam_name"] = "u_" + bgc_class + "_" + str(fam_id)
            df_families.loc[fam_id, "clusters_in_fam"] = len(family)

        for bgc in family:
            if bgc in df_clusters.index:
                df_clusters.loc[bgc, "fam_id_" + cutoff] = str(fam_id)
                if len(known_bgcs) > 0:
                    df_clusters.loc[bgc, "fam_type_" + cutoff] = "known_family"
                    known_compounds = ";".join(
                        df_known.loc[known_bgcs, "compounds"].tolist()
                    )
                    df_clusters.loc[
                        bgc, "fam_known_compounds_" + cutoff
                    ] = known_compounds
                else:
                    df_clusters.loc[bgc, "fam_type_" + cutoff] = "unknown_family"
                    df_clusters.loc[bgc, "fam_known_compounds_" + cutoff] = (
                        "u_" + bgc_class + "_" + str(fam_id)
                    )

            elif bgc in df_known.index:
                df_known.loc[bgc, "fam_id_" + cutoff] = str(fam_id)
                if len(known_bgcs) > 0:
                    df_known.loc[bgc, "fam_type_" + cutoff] = "known_family"
                    known_compounds = ";".join(
                        df_known.loc[known_bgcs, "compounds"].tolist()
                    )
                    df_known.loc[bgc, "fam_known_compounds_" + cutoff] = known_compounds
                else:
                    df_known.loc[bgc, "fam_type_" + cutoff] = "unknown_family"
                    df_known.loc[bgc, "fam_known_compounds_" + cutoff] = (
                        "u_" + bgc_class + "_" + str(fam_id)
                    )

    return df_clusters, df_known, df_families


def get_family_presence(df_clusters, df_genomes, df_families, cutoff):
    """
    Returns BGC family presence absence matrix across genomes
    """

    df_family_presence = pd.DataFrame(
        0, index=df_genomes.index, columns=df_families.index.astype(str)
    )
    for bgc_id in df_clusters.index:
        fam_id = str(df_clusters.loc[bgc_id, "fam_id_" + cutoff])
        genome_id = df_clusters.loc[bgc_id, "genome_id"]
        df_family_presence.loc[genome_id, fam_id] = 1

    return df_family_presence


def run_family_analysis(
    cutoff,
    net_data_path,
    df_clusters,
    df_genomes,
    df_known_all,
    output_dir,
    mibig_table,
    query_name,
):
    logging.info(f"Processing data from BiG-SCAPE with cutoff {cutoff}")
    df_network = get_bigscape_network(net_data_path, cutoff=cutoff)
    G_clusters, family_nodes = get_network_graph(df_network, weight="distance")
    df_network, df_known = remove_single_mibig(df_network, df_known_all, family_nodes)
    G_clusters, family_nodes = get_network_graph(df_network, weight="distance")
    singleton_bgc = [list(fam)[0] for fam in family_nodes if len(fam) == 1]
    family_graphs = get_family_graph(G_clusters)
    df_clusters, df_known, df_families = update_cluster_family(
        df_clusters, df_known, family_nodes, mibig_table, cutoff=cutoff
    )
    df_family_presence = get_family_presence(
        df_clusters, df_genomes, df_families, cutoff=cutoff
    )
    logging.debug(f"Number of genomes: {df_genomes.shape[0]}")
    logging.debug(f"Number of BGCs: {df_clusters.shape[0]}")
    logging.debug(f"Number of edges in the network: {df_network.shape[0]}")
    logging.debug(f"Number of total families: {len(family_nodes)}")
    logging.debug(
        f"Number of total non-single families: {len(family_nodes) - len(singleton_bgc)}"
    )
    logging.debug(f"Number of singleton families in the network: {len(singleton_bgc)}")
    logging.debug(
        f"Number of families with known BGCs: {df_families[df_families.fam_type == 'known_family'].shape[0]}"
    )
    logging.debug(f"Number of known BGCs in the network: {df_known.shape[0]}")
    logging.debug(
        f"Number of BGCs in top 10 families {[len(fam) for fam in family_graphs[:10]]}"
    )
    df_known_families = df_families[df_families.fam_type == "known_family"]
    logging.debug(
        f"Some of the common known BGCs{chr(10)}{chr(10)}{df_known_families.iloc[:20,1:3]}{chr(10)}"
    )
    df_unknown_families = df_families[df_families.fam_type == "unknown_family"]
    logging.debug(
        f"Some of the common unknown BGCs:{chr(10)}{chr(10)}{df_unknown_families.iloc[:20,1:3]}{chr(10)}"
    )

    # Assign index name
    df_network.index.name = "bigscape_edge_id"
    df_known.index.name = "bgc_id"
    df_families.index.name = f"fam_id_{cutoff}"
    df_clusters.index.name = "bgc_id"
    df_family_presence.index.name = "genome_id"

    # Save all the dataframes
    df_network.to_csv(f"{output_dir}/{query_name}_df_network_" + cutoff + ".csv")
    df_known.to_csv(f"{output_dir}/{query_name}_df_known_" + cutoff + ".csv")
    df_families.to_csv(f"{output_dir}/{query_name}_df_families_" + cutoff + ".csv")
    df_clusters.to_csv(f"{output_dir}/{query_name}_df_clusters_" + cutoff + ".csv")
    df_family_presence.to_csv(
        f"{output_dir}/{query_name}_df_family_presence_" + cutoff + ".csv"
    )

    return df_clusters, df_families, df_network


def process_bigscape_output(
    bigscape_directory, as_dir, df_genomes_path, mibig_bgc_table, output_dir
):
    bigscape_directory = Path(bigscape_directory)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # get the latest run
    bigscape_runs = {}
    logging.info("Getting BiG-SCAPE runs...")
    for net_data_path in bigscape_directory.glob("output_files/*"):
        logging.debug(f"Found network data: {net_data_path}")
        if net_data_path.is_dir() and "_c" in net_data_path.name:
            selected_run_folder = net_data_path.name
            selected_run_time_str = '_'.join(selected_run_folder.split("_c")[0].split('_')[-2:])
            selected_run_time = datetime.strptime(selected_run_time_str, "%Y-%m-%d_%H-%M-%S")
            bigscape_runs[selected_run_time] = selected_run_folder.split("_c")[0]  # the base
    run_times = list(bigscape_runs.keys())
    run_times.sort(reverse=True)

    # make sure run times has values
    selected_run = run_times[0]
    base = bigscape_runs[selected_run]

    logging.info(f"Processing {selected_run}")
    df_genomes = pd.read_csv(df_genomes_path)
    duplicated_ids = df_genomes[df_genomes.duplicated(subset=["genome_id"], keep=False)]["genome_id"].unique().tolist()
    if len(duplicated_ids) > 0:
        logging.warning(
            "Found duplicate genome_id values in genome summary input: %s",
            ", ".join(duplicated_ids),
        )
        df_genomes = df_genomes.drop_duplicates(subset=["genome_id"], keep="first")
    df_genomes = df_genomes.set_index("genome_id", drop=False)
    df_genomes.to_csv(f"{output_dir}/df_genome_antismash_summary.csv")

    # Process each cutoff
    for cutoff in ["0.3", "0.4", "0.5"]:
        net_data_path = bigscape_directory / "output_files" / (base + "_c" + cutoff)
        node_annot_path = net_data_path / "record_annotations.tsv"
        df_nodes = pd.read_csv(node_annot_path, index_col="Record", sep="\t")

        df_known_all, df_clusters = get_cluster_dataframes(
            df_genomes, df_nodes, mibig_bgc_table, as_dir
        )
        assert (
            len(df_clusters) > 0
        ), f"df_cluster is empty {df_clusters.shape}. Check folder data/interim/bgcs to have proper antismash region genbank files."
        # Enrich dataframes with BiGSCAPE information on GCCs and GCFs with cutoffs
        df_clusters, df_known_all = add_bigscape_families(
            df_clusters, df_known_all, net_data_path, cutoff
        )

        df_clusters, df_families, df_network = run_family_analysis(
            cutoff,
            net_data_path,
            df_clusters,
            df_genomes,
            df_known_all,
            output_dir,
            mibig_bgc_table,
            str(selected_run).replace(":", "_"),
        )
    return


if __name__ == "__main__":
    process_bigscape_output(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    )
