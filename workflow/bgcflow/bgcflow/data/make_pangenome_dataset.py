import os
import sys
import pandas as pd
from shutil import copyfile
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt


def get_roary_data(roary_interim_folder, roary_processed_folder):
    '''
    Copy important files from Roary interim to proecessed directory 

    Parameters
    ---------- 
    1. roary_folder : str / path
        Location of the output from Roary in the interim directory
    
    Returns
    -------
    1. roary_processed_folder : str / path
        Location of the processed output directory for important Roary results
    '''
    
    gene_presence_path = os.path.join(roary_interim_folder, 'gene_presence_absence.csv')
    df_gene_presence_summary = pd.read_csv(gene_presence_path, index_col='Gene')
    
    # Extract gene annotation columns to separate dataframe
    gene_summary_columns = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences',
        'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment',
        'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
        'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']

    gene_summary_out = os.path.join(roary_processed_folder, 'pangene_summary.csv')
    df_gene_summary = df_gene_presence_summary[gene_summary_columns].fillna('')

    # Add locus_tag and pangenome_id from pan_genome_reference.fa file
    pan_fasta_path = os.path.join(roary_interim_folder, 'pan_genome_reference.fa')
    records = SeqIO.parse(pan_fasta_path, format='fasta')
    for rec in records:
        locus_tag = rec.id
        pan_gene_desc = rec.description.split(' ')
        pan_gene_id = ' '.join([pan_gene_id for pan_gene_id in pan_gene_desc[1:]])
        if pan_gene_id in df_gene_summary.index:
            df_gene_summary.loc[pan_gene_id, 'locus_tag'] = locus_tag
    
    df_gene_summary.to_csv(gene_summary_out)

    # Extract locus tags 
    df_gene_presence = df_gene_presence_summary.drop(columns=gene_summary_columns).fillna('')
    gene_presence_out = os.path.join(roary_processed_folder, 'gene_presence_absence.csv')
    df_gene_presence.to_csv(gene_presence_out)

    # Save the gene presence absence binary matrix
    gene_presence_binary_path = os.path.join(roary_interim_folder, 'gene_presence_absence.Rtab')
    gene_presence_binary_out = os.path.join(roary_processed_folder, 'gene_presence_absence_binary.csv')
    df_gene_presence_binary = pd.read_csv(gene_presence_binary_path, index_col='Gene', sep='\t')
    df_gene_presence_binary.to_csv(gene_presence_binary_out)

    # Copy other output files to processed directory
    copyfile(os.path.join(roary_interim_folder, 'conserved_vs_total_genes.png'),
        os.path.join(roary_processed_folder, 'conserved_vs_total_genes.png'))

    copyfile(os.path.join(roary_interim_folder, 'Rplots.pdf'),
        os.path.join(roary_processed_folder, 'Rplots.pdf'))

    copyfile(os.path.join(roary_interim_folder, 'pan_genome_reference.fa'),
        os.path.join(roary_processed_folder, 'pan_genome_reference.fa'))
    
    copyfile(os.path.join(roary_interim_folder, 'summary_statistics.txt'),
        os.path.join(roary_processed_folder, 'summary_statistics.txt'))

    copyfile(os.path.join(roary_interim_folder, 'number_of_conserved_genes.Rtab'),
        os.path.join(roary_processed_folder, 'number_of_conserved_genes.Rtab'))
    
    copyfile(os.path.join(roary_interim_folder, 'number_of_genes_in_pan_genome.Rtab'),
        os.path.join(roary_processed_folder, 'number_of_genes_in_pan_genome.Rtab'))
    
    return df_gene_presence, df_gene_presence_binary, df_gene_summary


def plot_core_pan_curve(roary_interim_folder, roary_processed_folder):
    """
    Plots pangenome curves for core, pan, new and unique gene against number of genomes

    Parameters
    ---------- 
    1. roary_folder : str / path
        Location of the output from Roary in the interim directory
    
    Returns
    -------
    1. roary_processed_folder : str / path
        Location of the processed output directory where two plots will be saved
    """

    core_path = os.path.join(roary_interim_folder,'number_of_conserved_genes.Rtab')
    pan_path = os.path.join(roary_interim_folder,'number_of_genes_in_pan_genome.Rtab')
    new_path = os.path.join(roary_interim_folder,'number_of_new_genes.Rtab')
    uniq_path = os.path.join(roary_interim_folder,'number_of_unique_genes.Rtab')
    
    df_core_cnt = pd.read_table(core_path, sep='\t')
    df_pan_cnt = pd.read_table(pan_path, sep='\t')
    df_new_cnt = pd.read_csv(new_path, sep='\t')
    df_uniq_cnt = pd.read_csv(uniq_path, sep='\t')
    
    avg_core = df_core_cnt.sum(0)/df_core_cnt.shape[0]
    max_core = df_core_cnt.max(0)
    min_core = df_core_cnt.min(0)
    
    avg_pan = df_pan_cnt.sum(0)/df_core_cnt.shape[0]
    max_pan = df_pan_cnt.max(0)
    min_pan = df_pan_cnt.min(0)

    avg_new = df_new_cnt.sum(0)/df_uniq_cnt.shape[0]
    max_new = df_new_cnt.max(0)
    min_new = df_new_cnt.min(0)
        
    avg_uniq = df_uniq_cnt.sum(0)/df_uniq_cnt.shape[0]
    max_uniq = df_uniq_cnt.max(0)
    min_uniq = df_uniq_cnt.min(0)

    genome_numbers = list(range(1,df_core_cnt.shape[1]+1))
    
    # Plot pan and core genome curves in JPEG format
    fig_pan_core = go.Figure()
    fig_pan_core.add_trace(go.Scatter(x=genome_numbers, y=avg_core, fill='tonexty',
                             mode='lines+markers', name='Core', hoverinfo='name+x+y',
                             error_y=dict(type='data', symmetric=False, array=max_core - avg_core,
                              arrayminus=avg_core - min_core))) # fill down to xaxis
    
    fig_pan_core.add_trace(go.Scatter(x=genome_numbers, y=avg_pan, fill='tonexty',
                             mode='lines+markers', name='Pan', hoverinfo='name+x+y',
                            error_y=dict(type='data', symmetric=False, array=max_pan - avg_pan,
                              arrayminus=avg_pan - min_pan))) # fill to trace0 y
    
    fig_pan_core.update_layout(autosize=False, width=800, height=500, margin=dict(l=20, r=20, b=30, t=30, pad=4),
                      title_text='Pangenome and coregenome curve', xaxis_title_text='#Genomes', yaxis_title_text='#Genes', 
                      paper_bgcolor="White")

    fig_pan_core_path = os.path.join(roary_processed_folder, 'pan_core_curve.jpeg')
    if not os.path.isfile(fig_pan_core_path):
        fig_pan_core.write_image(fig_pan_core_path)
        
    # Plot new and unique genome curves in JPEG format
    fig_new_uniq = go.Figure()
    fig_new_uniq.add_trace(go.Scatter(x=genome_numbers, y=avg_new, fill='tonexty',
                             mode='lines+markers', name='New', hoverinfo='name+x+y',
                             error_y=dict(type='data', symmetric=False, array=max_new - avg_new,
                              arrayminus=avg_new - min_new))) # fill down to xaxis
    
    fig_new_uniq.add_trace(go.Scatter(x=genome_numbers, y=avg_uniq, fill='tonexty',
                             mode='lines+markers', name='Unique', hoverinfo='name+x+y',
                            error_y=dict(type='data', symmetric=False, array=max_uniq - avg_uniq,
                              arrayminus=avg_uniq - min_uniq))) # fill to trace0 y
    
    fig_new_uniq.update_layout(autosize=False, width=800, height=500, margin=dict(l=20, r=20, b=30, t=30, pad=4),
                      title_text='New and unique genes curve', xaxis_title_text='#Genomes', yaxis_title_text='#Genes', 
                      paper_bgcolor="White")

    fig_new_uniq_path = os.path.join(roary_processed_folder, 'new_unique_curve.jpeg')
    if not os.path.isfile(fig_new_uniq_path):
        fig_new_uniq.write_image(fig_new_uniq_path)

    return fig_pan_core.show(), fig_new_uniq.show()
    
def plot_pan_freq_plot(df_gene_presence_binary, roary_processed_folder):
    """
    Plots pangenome frequence plot for number of genes present in number of genomes

    Parameters
    ---------- 
    1. df_gene_presence_binary : pd.DataFrame
        Dataframe with presence absence matrix of genes in the pangenome
    
    Returns
    -------
    1. roary_processed_folder : str / path
        Location of the processed output directory where gene frequency plot will be saved
    """

    df_freq = pd.DataFrame(index=df_gene_presence_binary.index, columns=['#Genomes'])
    df_freq['#Genomes'] = df_gene_presence_binary.sum(axis=1)
    nbins = df_gene_presence_binary.shape[1]
    fig = px.histogram(df_freq, x="#Genomes", nbins=nbins)
    
    fig.update_layout(
        autosize=False, width=800, height=500, margin=dict(l=20, r=20, b=30, t=30, pad=4),
        title_text='Pangenome frequency plot', # title of plot
        xaxis_title_text='#Genomes', # xaxis label
        yaxis_title_text='#Genes', # yaxis label
        )
    fig_out_path = os.path.join(roary_processed_folder, 'gene_frequency.jpeg')
    if not os.path.isfile(fig_out_path):
        fig.write_image(fig_out_path)

    return fig.show()
    
    
def plot_pan_pie_chart(df_gene_presence_binary, roary_processed_folder, core=0.99, softcore=0.95, shell=0.15):
    """
    Plots pangenome frequence plot for number of genes present in number of genomes

    Parameters
    ---------- 
    1. df_gene_presence_binary : pd.DataFrame
        Dataframe with presence absence matrix of genes in the pangenome
    2. roary_processed_folder : str / path
        Location of the processed output directory where gene frequency plot will be saved
    3. core: int (0 to 1) 
        default = 0.99
        Fraction of genomes the gene must be present to be group in core (99% <= strains <= 100%)
    4. softcore: int (0 to 1) 
        default = 0.95
        Fraction of genomes the gene must be present to be group in softcore (95% <= strains < 99%)
    5. shell: int (0 to 1) 
        default = 0.15
        Fraction of genomes the gene must be present to be group in softcore (15% <= strains < 95%)
        Remaining genes will be in cloud
    
    Returns
    -------
    1. roary_processed_folder : str / path
        Location of the processed output directory where pangenome pie chart plot will be saved
    """

    roary = df_gene_presence_binary.copy()

    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    core = roary[roary.sum(axis=1) == roary.shape[1]].shape[0]
    softcore = roary[(roary.sum(axis=1) < roary.shape[1]) &
                     (roary.sum(axis=1) >= roary.shape[1]*0.85)].shape[0]
    shell = roary[(roary.sum(axis=1) < roary.shape[1]*0.85) &
                     (roary.sum(axis=1) >= roary.shape[1]*0.15)].shape[0]
    cloud = roary[roary.sum(axis=1) < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]
    
    def my_autopct(pct):
        val=int(pct*total/100.0)
        return '{v:d}'.format(v=val)

    fig = plt.pie([core, softcore, shell, cloud],
            labels=['core  (%d strains)'%roary.shape[1],
              'soft-core  (%d <= strains < %d)'%(roary.shape[1]*.85,
                                                 roary.shape[1]),
              'shell\n(%d <= strains < %d)'%(roary.shape[1]*.15,
                                             roary.shape[1]*.85),
              'cloud\n(strains < %d)'%(roary.shape[1]*.15)],
            explode=[0.1, 0.05, 0.02, 0], radius=0.9,
            colors=[(0, 0, 1, float(x)/total) for x in (core, softcore, shell, cloud)],
            autopct=my_autopct, textprops={'fontsize': 20, 'fontweight': 'bold'})
    
    fig_out_path = os.path.join(roary_processed_folder, 'pangenome_pie.jpeg')
    if not os.path.isfile(fig_out_path):
        plt.savefig(fig_out_path, facecolor='w', edgecolor='w', orientation='portrait', 
                format='jpeg', transparent=False, bbox_inches='tight', pad_inches=0.1)

    return plt.show()

if __name__ == "__main__":
    roary_interim_folder = sys.argv[1]
    roary_processed_folder = sys.argv[2]
    df_gene_presence, df_gene_presence_binary, df_gene_summary = get_roary_data(roary_interim_folder, roary_processed_folder)
    plot_core_pan_curve(roary_interim_folder, roary_processed_folder)
    plot_pan_freq_plot(df_gene_presence_binary, roary_processed_folder)
    plot_pan_pie_chart(df_gene_presence_binary, roary_processed_folder, core=0.99, softcore=0.95, shell=0.15)
