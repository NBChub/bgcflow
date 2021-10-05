import matplotlib.pyplot as plt # Visualization 
from matplotlib.ticker import StrMethodFormatter 

def plot_hist(df, column_name, bin_size, title_str, xlabel_str, ylabel_str, 
              to_path='../figures/bgc_dist.pdf'):
    
    plt.figure(figsize=(12,8))
    ax = df.hist(column= column_name, bins=bin_size, grid=False, color='#86bf91', zorder=2, rwidth=0.9)
    ax = ax[0]
    for x in ax:
        # Despine
        x.spines['right'].set_visible(False)
        x.spines['top'].set_visible(False)
        x.spines['left'].set_visible(False)

        # Switch off ticks
        x.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", 
                      left="off", right="off", labelleft="on")

        # Draw horizontal axis lines
        vals = x.get_yticks()
        for tick in vals:
            x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)

        # Remove title
        x.set_title(title_str, weight='bold', size=16)

        # Set x-axis label
        x.set_xlabel(xlabel_str, labelpad=20, weight='bold', size=12)

        # Set y-axis label
        x.set_ylabel(ylabel_str, labelpad=20, weight='bold', size=12)

        # Format y-axis label
        x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))
    plt.savefig(to_path, facecolor='w', edgecolor='w', orientation='portrait', 
                format='pdf', transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    
    return None


def plot_bgc_dist(df_genomes, col_select='bgcs_count', to_path='../figures/bgc_dist.pdf'):
    """
    Plot number of BGCs distribution over genomes
    """
    
    d = df_genomes[col_select].fillna(0).tolist()
    bin_size = int(max(d))
    title_str = 'Distribution of records per genome: '
    xlabel_str = col_select
    ylabel_str = 'No of genomes'
    
    plot_hist(df_genomes, col_select, bin_size, title_str, xlabel_str, ylabel_str, to_path)
    
    return None

def scatter_bgcs_len(df_genomes,to_path='../figures/scatter_bgcs_len.pdf'):
    
    scat_data = df_genomes[['genome_len', 'bgcs_count']].astype(float)
    scat_data['genome_len'] = round(scat_data['genome_len']/1000000,1)
    #     scat_data['GC Content'] = round(scat_data['GC Content'],1)
    scat_data['genus'] = df_genomes.genus
        
    df_counts = scat_data.groupby(['genome_len', 'bgcs_count',
                                   'genus']).size().reset_index(name='counts')
        
    # Create Fig and gridspec
    fig = plt.figure(figsize=(16, 10), dpi= 80)
    grid = plt.GridSpec(4, 4, hspace=0.6, wspace=0.2)

    # Define the axes
    ax_main = fig.add_subplot(grid[:-1, :-1])
    ax_right = fig.add_subplot(grid[:-1, -1], yticklabels=[])
    ax_bottom = fig.add_subplot(grid[-1, 0:-1], xticklabels=[])

    # Scatterplot on main ax
    ax_main.scatter('genome_len', 'bgcs_count', data=df_counts, s=df_counts.counts*40, 
                    alpha = 0.55, edgecolors='gray', linewidths=.5)
    # ax_main.scatter('GC Content', 'no_of_clusters', data=df_counts, s=df_counts.counts*33, c=cnt_colors, alpha = 0.55,
    #                edgecolors='gray', linewidths=.5)

    # histogram on the right
    # ax_bottom.hist(df['GC Content'], 40, histtype='stepfilled', orientation='vertical', color='deeppink')
    ax_bottom.hist(scat_data['genome_len'], 40, histtype='stepfilled', orientation='vertical',
                   color='deeppink')
    ax_bottom.invert_yaxis()

    # histogram in the bottom
    ax_right.hist(scat_data['bgcs_count'], int(max(scat_data.bgcs_count) + 1), 
                  histtype='stepfilled', orientation='horizontal', color='deeppink')

    # Decorations
    # ax_main.set(title='No. of clusters vs GC Content', 
    #             xlabel='GC Content(%)', ylabel='No of clusters per genome')
    ax_main.set(title='No. of clusters vs Genome Length', 
                xlabel='Genome Length (Mb)', ylabel='No of clusters per genome')
    ax_main.title.set_fontsize(24)

    ax_bottom.set(ylabel='No of genomes')
    ax_right.set(xlabel='No of genomes')

    # Annotate 
    from matplotlib.patches import Ellipse
    for item in ([ax_main.xaxis.label, ax_main.yaxis.label] + ax_main.get_xticklabels() + ax_main.get_yticklabels()):
        item.set_fontsize(18)
    for item in ([ax_right.xaxis.label, ax_right.yaxis.label] + ax_right.get_xticklabels() + ax_right.get_yticklabels()):
        item.set_fontsize(18)
    for item in ([ax_bottom.xaxis.label, ax_bottom.yaxis.label] + ax_bottom.get_xticklabels() + ax_bottom.get_yticklabels()):
        item.set_fontsize(18)

    plt.savefig(to_path, facecolor='w', edgecolor='w', orientation='portrait', 
                format='pdf', transparent=False, bbox_inches='tight', pad_inches=0.1)
    return plt.show()