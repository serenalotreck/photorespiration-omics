"""
Code to make an arrow table of opposite expression.

Author: Serena G. Lotreck
"""
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

def plot_opposite_expression(log2fc_df, cols_to_plot, gene_name_map=None,
                             scale_factor=50, title=None, dark=False):
    """
    Plot arrows for opposite expression. Assumes expression log2FC has been
    normalized between -1 and 1 and that the log2fc_df index is in sequential
    order starting with 0.

    parameters:
        log2fc_df, df: dataframe with log2FoldChange values to plot
        cols_to_plot, dict: values are the semantic column header name
            to be used, values are the column names in log2fc_df
        gene_name_map, dict: keys are TAIR locus ID's, values are some
            semantic name to be added to the table
        scale_factor, int: amount to multiply the normalized expression value
            to determine arrow thickness
        title, str: optional, string to use for the title
        dark, bool: whether or not to use matplotlib dark background

    returns: None
    """
    if dark:
        plt.style.use('dark_background')
    else:
        plt.rcdefaults()

    fig, ax = plt.subplots(figsize=(8,6))

    rows = len(log2fc_df)
    cols = len(cols_to_plot) + 1
    
    ax.set_ylim(-1, rows + 1)
    ax.set_xlim(0, cols + .5)
    
    for i, row in log2fc_df.iterrows():
    
        # Add gene name
        if gene_name_map == None:
            s = row.gene_id
        else:
            try:
                s = r"$\bf{" + gene_name_map[row.gene_id] + "}$ (" + row.gene_id + ')'
            except KeyError:
                s = row.gene_id
        ax.text(x=0.5, y=i, s=s)
    
        # Add arrows
        for j, to_plot in enumerate(cols_to_plot.keys()):
            
            arrow_dir = 'up' if row[cols_to_plot[to_plot]] > 0 else 'down'
            arrow_len = abs(row[cols_to_plot[to_plot]])
            if arrow_dir == 'up':
                start_coords = (j + 2, i - 0.5*arrow_len)
                end_coords = (j + 2, i + 0.5*arrow_len)
                color = 'blue'
            else:
                start_coords = (j + 2, i + 0.5*arrow_len)
                end_coords = (j + 2, i - 0.5*arrow_len)
                color = 'red'
            arrow = mpatches.FancyArrowPatch(start_coords, end_coords,
                                         mutation_scale=abs(row[cols_to_plot[to_plot]])*scale_factor,
                                         color=color)
            ax.add_patch(arrow)
    
    # Add column headers
    ax.text(0.7, rows - 0.25, 'Gene', weight='bold', ha='center')
    for k, col in enumerate(cols_to_plot.keys()):
        ax.text(k + 2, rows - 0.25, col, weight='bold', ha='center')
    
    # Add gridlines
    for row in range(rows):
        ax.plot(
        	[0, cols + 1],
        	[row -.5, row - .5],
        	ls=':',
        	lw='.5',
        	c='grey'
        )
    ax.plot([0, cols + 1], [rows - 0.5, rows - 0.5], lw='.5', c='black')
    
    # Turn of axes
    ax.axis('off')

    # Add title if requested
    if title is not None:
        ax.set_title(
        	title,
        	loc='center',
        	fontsize=16,
        	weight='bold'
        )