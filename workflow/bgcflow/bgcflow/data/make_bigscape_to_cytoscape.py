import pandas as pd
import os, sys
from pathlib import Path
import logging
import networkx as nx 
from collections import OrderedDict
from alive_progress import alive_bar

log_format = '%(levelname)-8s %(asctime)s   %(message)s'
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

def get_cluster_dataframes(df_genomes, df_nodes, as_dir = '../data/antismash/'):
    '''
    Returns two dataframes of clusters with information from genomes and MIBIG
    '''
    df_clusters = pd.DataFrame(columns=['product', 'bigscape_class', 'genome_id', 'genome_name', 'genus',  'accn_id'])
    df_known = pd.DataFrame(columns=['product', 'bigscape_class', 'genome_id', 'genome_name', 'genus',  'accn_id', 'compound'])
    genus_list = df_genomes.genus.unique().tolist()
    
    # Generate bgcs dataframe with metadata from df_nodes and df_genomes
    logging.info('Generating bgcs dataframe with metadata from df_nodes and df_genomes...')
    with alive_bar(len(df_genomes.index)) as bar:
        for genome_id in df_genomes.index:
            logging.info(f'Processing BGCs in the genome: {genome_id}')
            if genome_id in os.listdir(as_dir):
                genome_dir = os.path.join(as_dir, genome_id)
                bgc_id_list = [region[:-4] for region in os.listdir(genome_dir) if 'region0' in region]
                for bgc_id in bgc_id_list:
                    if bgc_id in df_nodes.index:
                        df_clusters.loc[bgc_id, 'genome_id'] = genome_id
                        df_clusters.loc[bgc_id, 'product'] = df_nodes.loc[bgc_id, 'Product Prediction']
                        df_clusters.loc[bgc_id, 'bigscape_class'] = df_nodes.loc[bgc_id, 'BiG-SCAPE class']
                        df_clusters.loc[bgc_id, 'accn_id'] = df_nodes.loc[bgc_id, 'Accession ID']
                        #df_clusters.loc[bgc_id, 'genome_name'] = df_genomes.loc[genome_id, 'genome_name']
                        df_clusters.loc[bgc_id, 'genus'] = df_genomes.loc[genome_id, 'genus']
                        df_clusters.loc[bgc_id, 'genome_len'] = df_genomes.loc[genome_id, 'genome_len']
                        df_clusters.loc[bgc_id, 'bgcs_count'] = df_genomes.loc[genome_id, 'bgcs_count']
                    else:
                        logging.debug(f'{bgc_id} not in df_nodes')
            else:
                logging.warning(f'{genome_id} not in directory!')
            bar()
    
    # Generate separate table for known BGCs from MIBIG
    logging.info('Generating separate table for known BGCs from MIBIG')
    with alive_bar(len(df_nodes.index)) as bar:
        for bgc_id in df_nodes.index:
            if 'BGC0' in bgc_id:
                df_known.loc[bgc_id, 'product'] = df_nodes.loc[bgc_id, 'Product Prediction']
                df_known.loc[bgc_id, 'bigscape_class'] = df_nodes.loc[bgc_id, 'BiG-SCAPE class']
                #df_known.loc[bgc_id, 'genome_id'] = df_nodes.loc[bgc_id, 'Accesion ID']
                genome_name = df_nodes.loc[bgc_id, 'Organism']
                df_known.loc[bgc_id, 'genome_name'] = genome_name
                genus = genome_name.split(' ')[0]
                df_known.loc[bgc_id, 'genus'] = genus
                #df_known.loc[bgc_id, 'accn_id'] = df_nodes.loc[bgc_id, 'Accesion ID']
                desc = df_nodes.loc[bgc_id, 'Description']
                compound_name = desc.split('biosynthetic gene cluster')[0].strip()
                df_known.loc[bgc_id, 'compound'] = compound_name
            bar()
        
    return df_known, df_clusters

def add_bigscape_families(df_clusters, df_known, net_data_path):
    '''
    Adds GCC and GCF numbers detected by BiG-SCAPE clustering for different cut-offs
    '''
    
    cluster_class_set = [cluster_class for cluster_class in os.listdir(net_data_path)
                         if '.tsv' not in cluster_class]
    
    for cluster_class in cluster_class_set:
        logging.info(f'Processing all BGCs from {cluster_class}')
        class_dir = os.path.join(net_data_path, cluster_class)
        gcf_cutoffs_files = [file for file in os.listdir(class_dir) if '_clustering_' in file]
        for gcf_file in gcf_cutoffs_files:
            cutoff = gcf_file[-8:-4]
            gcf_path = os.path.join(class_dir, gcf_file)
            df_clusters = read_gcf_df(df_clusters, gcf_path, cutoff)
            df_known = read_gcf_df(df_known, gcf_path, cutoff)
            
        clan_files = [file for file in os.listdir(class_dir) if '_clans_' in file]
        if len(clan_files) == 1:
            clan_select = clan_files[0]
            clan_path = os.path.join(class_dir, clan_select)
                   
            df_clusters = read_gcc_df(df_clusters, clan_path)
            df_known = read_gcc_df(df_known, clan_path)
            
    return df_clusters, df_known

def read_gcf_df(df_clusters, gcf_path, cutoff):
    '''
    Adds GCF (Gene Cluster Family) number for each BGC
    '''
        
    df_gcf = pd.read_csv(gcf_path, sep = '\t', index_col = '#BGC Name', dtype = str)
    col_name = 'gcf_' + cutoff
    
    for bgc in df_clusters.index:
        if bgc in df_gcf.index:
            df_clusters.loc[bgc, col_name] = df_gcf.loc[bgc, 'Family Number']
            
             
    return df_clusters

def read_gcc_df(df_clusters, clan_path):
    '''
    Adds GCC (Gene Cluster Clan) number for each BGC
    '''
    
    df_gcc = pd.read_csv(clan_path, sep = '\t', index_col = '#BGC Name', dtype=str)
    col_name = 'Clan Number'
    
    for bgc in df_clusters.index:
        if bgc in df_gcc.index:
            df_clusters.loc[bgc, col_name] = df_gcc.loc[bgc, 'Clan Number']
            
    return df_clusters

def run_family_analysis(cutoff, net_data_path, df_clusters, df_genomes, df_known_all, output_dir, query_name):
    logging.info(f'Processing data from BiG-SCAPE with cutoff {cutoff}')
    df_network = get_bigscape_network(net_data_path, cutoff = cutoff)
    G_clusters, family_nodes = get_network_graph(df_network, weight = 'Raw distance')
    df_network, df_known = remove_single_mibig(df_network, df_known_all, family_nodes)
    G_clusters, family_nodes = get_network_graph(df_network, weight = 'Raw distance')
    singleton_bgc = [list(fam)[0] for fam in family_nodes if len(fam) == 1]
    family_graphs = get_family_graph(G_clusters)
    df_clusters, df_known, df_families = update_cluster_family(df_clusters, df_known, family_nodes, cutoff = cutoff)
    logging.debug(f'Number of genomes: {df_genomes.shape[0]}')
    logging.debug(f'Number of BGCs: {df_clusters.shape[0]}')
    logging.debug(f'Number of edges in the network: {df_network.shape[0]}')
    logging.debug(f'Number of total families: {len(family_nodes)}')
    logging.debug(f'Number of total non-single families: {len(family_nodes) - len(singleton_bgc)}')
    logging.debug(f'Number of singleton families in the network: {len(singleton_bgc)}')
    logging.debug(f"Number of families with known BGCs: {df_families[df_families.fam_type == 'known_family'].shape[0]}")
    logging.debug(f'Number of known BGCs in the network: {df_known.shape[0]}')
    # print('BGCs in largest families:', family_size)
    logging.debug(f'Number of BGCs in top 10 families {[len(fam) for fam in family_graphs[:10]]}')
    df_known_families = df_families[df_families.fam_type == 'known_family']
    logging.debug(f'Some of the common known BGCs{chr(10)}{chr(10)}{df_known_families.iloc[:20,1:3]}{chr(10)}')
    df_unknown_families = df_families[df_families.fam_type == 'unknown_family']
    logging.debug(f'Some of the common unknown BGCs:{chr(10)}{chr(10)}{df_unknown_families.iloc[:20,1:3]}{chr(10)}')
    # Save all the dataframes
    df_network.to_csv(f'{output_dir}/{query_name}_df_network_' + cutoff + '.csv')
    df_known.to_csv(f'{output_dir}/{query_name}_df_known_all_' + cutoff + '.csv')
    df_families.to_csv(f'{output_dir}/{query_name}_df_families_' + cutoff + '.csv')
    df_clusters.to_csv(f'{output_dir}/{query_name}_df_clusters_' + cutoff + '.csv')
    
    return df_clusters, df_families, df_network

def get_bigscape_network(net_data_path, cutoff = '0.30'):
    '''
    Reads similarity network for a particular to a dataframe
    '''
   
    cluster_class_set = [cluster_class for cluster_class in os.listdir(net_data_path)
                         if '.tsv' not in cluster_class]
    
    df_network = pd.DataFrame()
    for cluster_class in cluster_class_set:
        class_dir = os.path.join(net_data_path, cluster_class)
        net_file = cluster_class + '_c' + cutoff + '.network'
        net_path = os.path.join(class_dir, net_file)
        df_class_net = pd.read_csv(net_path, sep='\t')
        df_network = pd.concat([df_network, df_class_net], ignore_index = True)
    
    return df_network

def get_network_graph(df_network, weight = 'Raw distance'):
    '''
    Returns networkX graph for a given network
    '''
    
    G_clusters = nx.from_pandas_edgelist(df_network, 'Clustername 1', 'Clustername 2', weight)
    G_families = nx.connected_components(G_clusters)
    family_nodes = [c for c in sorted(G_families, key=len, reverse=True)]
    
    return G_clusters, family_nodes

def remove_single_mibig(df_network, df_known, family_nodes):
    '''
    Removes singleton MIBIG BGCs from the network
    '''
    logging.info("Removing singleton MIBIG BGCs from the network")
    nodes_to_remove = []
    
    for fam in family_nodes:
        single_mibig = True
        
        for node in fam:
            if 'BGC' not in node:
                single_mibig = False
        
        if single_mibig:
            for node in fam:
                if node not in nodes_to_remove:
                    nodes_to_remove.append(node)
    logging.debug(f"{len(nodes_to_remove)} number of MIBIG nodes will be removed from analysis due to no similarity")
    df_network = df_network[~df_network['Clustername 1'].isin(nodes_to_remove)]
    df_network = df_network[~df_network['Clustername 2'].isin(nodes_to_remove)]
    for node in nodes_to_remove:
        if node in df_known.index:
            df_known = df_known.drop(node)

    return df_network, df_known

def get_family_graph(G_clusters):
    '''
    Returns families list as networkx graph 
    '''

    # Find connected components or cluster families
    Families_list = list(nx.connected_components(G_clusters))
    # Sort families in decreasing order of occurrence
    family_size = [len(family) for family in Families_list]
    orig_index = list(range(len(family_size)))
    orig_index_fam_size_dict = dict(zip(orig_index, family_size))

    sorted_index_fam_size_dict = OrderedDict(sorted(orig_index_fam_size_dict.items(), 
                                                    key=lambda x: x[1]))
    new_index = list(range(len(family_size)-1,-1,-1))
    orig_index_sorted = list(sorted_index_fam_size_dict.keys())
    new_orig_index_dict = dict(zip(new_index, orig_index_sorted))

    # Ordered family graphs
    family_graphs = [Families_list[new_orig_index_dict[fam_id]] for fam_id in range(len(Families_list))]

    return family_graphs

def update_cluster_family(df_clusters, df_known, family_nodes, cutoff = '0.30'):
    '''
    Updates df_clusters with family ids (connected components)
    '''
    
    df_families = pd.DataFrame(columns=['fam_type', 'fam_name', 'clusters_in_fam'])
    for cntr in range(len(family_nodes)):
        fam_id = cntr + 1
        family = family_nodes[cntr]
        known_bgcs = [bgc for bgc in family if 'BGC' in bgc]

        if len(known_bgcs) > 0:
            df_families.loc[fam_id, 'fam_type'] = 'known_family'
            known_compounds = ', '.join(df_known.loc[known_bgcs, 'compound'].tolist())
            df_families.loc[fam_id, 'fam_name'] = known_compounds
            df_families.loc[fam_id, 'clusters_in_fam'] = len(family)
        else:
            df_families.loc[fam_id, 'fam_type'] = 'unknown_family'
            bgc_class = ','.join(df_clusters.loc[list(family), 'bigscape_class'].unique().tolist())
            df_families.loc[fam_id, 'fam_name'] = 'u_' + bgc_class + '_' + str(fam_id)
            df_families.loc[fam_id, 'clusters_in_fam'] = len(family)

        for bgc in family:
            if bgc in df_clusters.index:
                df_clusters.loc[bgc, 'fam_id_' + cutoff] = str(fam_id)
                if len(known_bgcs) > 0:
                    df_clusters.loc[bgc, 'fam_type_' + cutoff] = 'known_family'
                    known_compounds = ', '.join(df_known.loc[known_bgcs, 'compound'].tolist())
                    df_clusters.loc[bgc, 'fam_known_compounds_' + cutoff] = known_compounds
                else:
                    df_clusters.loc[bgc, 'fam_type_' + cutoff] = 'unknown_family'
                    df_clusters.loc[bgc, 'fam_known_compounds_' + cutoff] = 'u_' + bgc_class + '_' + str(fam_id)
                    
            elif bgc in df_known.index:
                df_known.loc[bgc, 'fam_id'] = str(fam_id)
    
    return df_clusters, df_known, df_families

def process_bigscape_output(bigscape_directory, as_dir, df_genomes_path, output_dir):
    bigscape_directory = Path(bigscape_directory)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for net_data_path in bigscape_directory.glob("network_files/*glocal*"):
        selected_run = net_data_path.stem
        logging.info(f'Processing {selected_run}')
        node_annot_path = net_data_path / 'Network_Annotations_Full.tsv' # Read the BGC table

        df_nodes = pd.read_csv(node_annot_path, index_col='BGC', sep='\t')

        # Generate df_clusters and df_known dataframe
        df_genomes = pd.read_csv(df_genomes_path, index_col=0)

        df_known_all, df_clusters = get_cluster_dataframes(df_genomes, df_nodes, as_dir)
        # Enrich dataframes with BiGSCAPE information on GCCs and GCFs with cutoffs
        df_clusters, df_known = add_bigscape_families(df_clusters, df_known_all, net_data_path)

        # Get GCF data as per the cutoff
        for cutoff in  ['0.30', '0.40', '0.50']:
            df_clusters, df_families, df_network = run_family_analysis(cutoff, net_data_path, df_clusters, df_genomes, df_known_all, output_dir, selected_run)
        return

if __name__ == "__main__":
    process_bigscape_output(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
