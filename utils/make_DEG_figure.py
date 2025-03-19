"""
Code to make omics expression figures.

Author: Serena G. Lotreck
"""
import pandas as pd
import matplotlib.pyplot as plt


def makeDEGfigure(deg_dfs, photosynth_sets, photosynth_colors, semantic_names, tair2gene,id_col='gene_id', title_name='genes'):
    """
    Maked stacked expresion figure.

    parameters:
        deg_dfs, dict: keys are comparison names, values are dfs with
            log2FC values for significant DEGs
        photosynth_sets, dict: keys are set names, values are lists of
            gene names
        photosynth_colors, dict: keys are photosynth gene categories,
            values are colors to use
        semantic_names, dict: keys are the keys of deg_dfs, values are
            semantic names to be used as subplot titles
        tair2gene, dict: keys are lwoercased TAIR ID's of the genes in the
            plot, values are the gene names
        id_col, str: name of the column with the gene IDs
    """
    # Capitalize the photosynthesis 
    # Get the overall list of all photosynth genes that appear in any DEG set
    # These will be the x axis for all plots
    all_photosynth = [str(g).upper() for s in photosynth_sets.values() for g in s]
    all_degs = [str(g).upper() for g_list in [df[id_col].tolist() for df in deg_dfs.values()] for g in g_list]
    x_genes = set(all_photosynth).intersection(all_degs)
    
    # Make a color dict for all the sets
    g_to_group = {}
    for ph, gs in photosynth_sets.items():
        for g in gs:
            g_to_group[g] = ph
    color_dict = {g: photosynth_colors[ph] for g, ph in g_to_group.items()}
    
    # Get the intersections of the photosynth sets with the DEG conditions
    to_plot = {}
    for deg_set_name, deg_set in deg_dfs.items():
        photosynth = deg_set[deg_set[id_col].isin(all_photosynth)][[id_col, 'log2FoldChange']]
        not_present = {id_col: [g for g in x_genes if g not in photosynth[id_col].tolist()]}
        not_present['log2FoldChange'] = [0]* len(not_present[id_col])
        to_plot_df = pd.concat([photosynth, pd.DataFrame(not_present)])
        to_plot_df[id_col] = to_plot_df[id_col].str.lower()
        to_plot[deg_set_name] = to_plot_df

    # Plot
    fig, axs = plt.subplots(len(to_plot), sharex=True, sharey=True, figsize=(20, len(to_plot)*5))

    for deg_set_name, ax in zip(to_plot.keys(), axs.flat):
        
        current_set = to_plot[deg_set_name]
        nonzero = len(current_set[current_set['log2FoldChange'] != 0])
        labels = [g_to_group[g] for g in current_set[id_col]]
        colors = [color_dict[g] for g in current_set[id_col]]
        ax.bar(current_set[id_col], current_set['log2FoldChange'], label=labels, color=colors)
        ax.set_title(semantic_names[deg_set_name] + f': {nonzero} DE {title_name}')
        ax.tick_params(axis='x', labelrotation=90, labelbottom=True)
        xtick_labels = ax.get_xticklabels()

    # Write over xtick labels with extra semantic stuff
    semantic_xticklabels = []
    for g in xtick_labels:
        try:
            semantic_xticklabels.append(f'{g.get_text().upper()}' + ' (' + r"$\bf{" + f'{tair2gene[g.get_text()]}' + "}$" + ')') 
        except KeyError:
            semantic_xticklabels.append(f'{g.get_text().upper()}')
    for ax in axs.flat:
        ax.set_xticklabels(semantic_xticklabels, rotation=90)

    fig.supylabel('log2FoldChange', x=ax.get_position().x0 - 0.05)
    
    # To remove duplicate labels
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    unique_labels = set(labels)
    legend_dict = dict(zip(labels, lines))
    unique_lines = [legend_dict[x] for x in unique_labels]
    legend_x = axs.flatten()[-1].get_position().xmax
    legend_y = axs.flatten()[0].get_position().ymax
    fig.legend(unique_lines, unique_labels, loc=(8.75*legend_x/10, legend_y + legend_y/20)) # This may not generalize well, did by trial and error

    plt.subplots_adjust(hspace=1)