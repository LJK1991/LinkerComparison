import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import seaborn as sns
import pandas as pd

'''
ander kleur, de rest is ok.
'''

def main():
    mplstyle.use('fast')
    flist = ['/media/draco/lucask/linkercomparison/SPLIT/SPLiT_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/GtiT/GtoT_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/SHARE/SHARE_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/parse/parse_seq_linker_alignments.tsv']

    ABCD = ['A - SPLiT','B - G to T','C - SHARE','D - Parse']
    labs = ["Correctable\nlinker", "Uncorrectable\nlinker 1", "Uncorrectable\nlinker 2", "Uncorrectable\nlinker"]

    fig = plt.figure(figsize=(24, 16), constrained_layout=True)
    subfig = fig.subfigures(4, 2, hspace=0)

    for i in range(0, len(flist)):
        df = pd.read_csv(flist[i], sep='\t', usecols=[5, 6, 7, 8, 9, 10, 17])

        subfignest = subfig[i, 0].subfigures(1, 2, wspace=0, width_ratios=[1, 2])
        subfig[i, 0].suptitle(ABCD[i], x=0.05, ha='left', fontweight='bold')

        ax1 = subfignest[0].subplots(1, 1)
        ax2 = subfignest[1].subplots(1, 2, sharey=True, gridspec_kw={'wspace': 0})
        sns.violinplot(ax=ax1, x=df['linker_relation'], y=df['l2_len'], palette='tab10', scale='width')
        sns.violinplot(ax=ax2[0], x=df['linker_relation'], y=df['l2_start'], palette='tab10', scale='width')
        sns.violinplot(ax=ax2[1], x=df['linker_relation'], y=df['l2_stop'], palette='tab10', scale='width')
        ax1.set_title('linker 2 length')
        ax2[0].set_title('linker 2 start point')
        ax2[1].set_title('linker 2 stop point')
        ax1.set(ylabel='Length in bp', xlabel=None)
        ax2[0].set(ylabel='Base in read', xlabel=None)
        ax2[1].set(ylabel=None, xlabel=None)
        ax2[0].set_yticks([0, 18, 48, 56, 86, 100])
        ax2[1].get_yaxis().set_visible(False)
        ax1.set_xticklabels(labs, fontdict={'fontsize': 6})
        ax2[0].set_xticklabels(labs, fontdict={'fontsize': 6})
        ax2[1].set_xticklabels(labs, fontdict={'fontsize': 6})


        subfignest = subfig[i, 1].subfigures(1, 2, wspace=0, width_ratios=[1, 2])
        ax1 = subfignest[0].subplots(1, 1)
        ax2 = subfignest[1].subplots(1, 2, sharey=True, gridspec_kw={'wspace': 0})
        sns.violinplot(ax=ax1, x=df['linker_relation'], y=df['l1_len'], palette='tab10', scale='width')
        sns.violinplot(ax=ax2[0], x=df['linker_relation'], y=df['l1_start'], palette='tab10', scale='width')
        sns.violinplot(ax=ax2[1], x=df['linker_relation'], y=df['l1_stop'], palette='tab10', scale='width')
        ax1.set_title('linker 1 length')
        ax2[0].set_title('linker 1 start point')
        ax2[1].set_title('linker 1 stop point')
        ax1.set(ylabel='Length in bp', xlabel=None)
        ax2[0].set(ylabel='Base in read', xlabel=None)
        ax2[1].set(ylabel=None, xlabel=None)
        ax2[0].set_yticks([0, 18, 48, 56, 86, 100])
        ax2[1].get_yaxis().set_visible(False)
        ax1.set_xticklabels(labs, fontdict={'fontsize': 6})
        ax2[0].set_xticklabels(labs, fontdict={'fontsize': 6})
        ax2[1].set_xticklabels(labs, fontdict={'fontsize': 6})

    fig.savefig('/media/draco/lucask/linkercomparison/Supp_Fig_x.png')

if __name__ == '__main__': main()
