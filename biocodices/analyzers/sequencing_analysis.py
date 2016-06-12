import matplotlib.pyplot as plt
import seaborn as sns


def plot_vcf_stats(vcf_stats_df, ax=None, palette_name='husl',
                   legend_title=None, figsize=(15, 5)):
    """Expects a DataFrame of VCF stats produced by the function vcf_stats.
    Will plot the values and draw lines to separate the chromosomes.
    Returns a matplotlib axes."""
    sns.set_palette(palette_name)
    sns.set_style('white')

    if not ax:
        _, ax = plt.subplots(1)

    chromosomes = vcf_stats_df.index.get_level_values(0)
    xticks_per_chrom = chromosomes.value_counts().sort_index()
    xlabels_offsets = xticks_per_chrom.cumsum() - (xticks_per_chrom / 2)
    # ^ Puts the chromsome labels in the middle of the chromosome section
    # of the plot.

    vcf_stats_df.plot(ax=ax, figsize=figsize, linewidth=2)
    sample_names = [col[0] for col in vcf_stats_df.columns]
    # ^ Columns are a multi-index, the first level is the sample name.
    ax.legend(sample_names, ncol=2, loc='best', title=legend_title)
    ax.set_xticks(xlabels_offsets)
    ax.set_xticklabels(xlabels_offsets.index)
    ax.set_xlabel('Chromosome')

    for limit in xticks_per_chrom.cumsum():
        ax.axvline(limit, color='#d9d9d9')

    sns.despine(left=True)
    return ax
