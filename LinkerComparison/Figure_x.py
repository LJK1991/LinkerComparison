import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
import numpy as np

'''
verander:
legende subfiguur B: 1 x en groter
kleur voor alle subfiguren anders
kleur subfiguur C specifiek anders dan B
fontsize van titels/andere text groter.
verwijder subfiguur A
'''

def main():

    mplstyle.use('fast')
    flist = ['/media/draco/lucask/linkercomparison/SPLIT/SPLiT_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/GtiT/GtoT_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/SHARE/SHARE_seq_linker_alignments.tsv',
             '/media/draco/lucask/linkercomparison/parse/parse_seq_linker_alignments.tsv']

    fig = plt.figure(figsize=(24, 16), constrained_layout=True)
    subfig = fig.subfigures(1, 2, width_ratios=[1, 1])


    plot_A = subfig[0].subplots(4, 2,gridspec_kw={'wspace': 0, 'hspace': 0})
    subfig[0].suptitle('A', x=0.05, ha='left', fontweight='bold', fontsize=18)

    for i in range(0, len(flist)):
        name = ['SPLiT', 'G to T','SHARE', 'Parse']
        df = pd.read_csv(flist[i], sep='\t', usecols=[3, 4, 17])
        labs = ["Correctable linker", "Uncorrectable linker 1", "Uncorrectable linker 2", "Uncorrectable linker"]
        #fig, axes = plt.subplots(1, 2, figsize=(10 ,5), sharey=True, layout='tight')
        sns.histplot(ax=plot_A[i, 1], data=df, x='l1_score', multiple='stack', hue='linker_relation', hue_order=[1, 2, 3, 4], legend=False, element='bars', binwidth=1, binrange=(0, 30), palette='tab10')
        sns.histplot(ax=plot_A[i, 0], data=df, x='l2_score', multiple='stack', hue='linker_relation', hue_order=[1, 2, 3, 4], legend=True, element='bars', binwidth=1, binrange=(0, 30), palette='tab10')
        legend = plot_A[i, 0].get_legend()
        handles = legend.legendHandles
        legend.remove()
        if i == 0:
            plot_A[i, 0].xaxis.set_label_position('top')
            plot_A[i, 0].set_ylabel('{}\nCount'.format(name[i]))
            plot_A[i, 0].set_xlabel('linker 2 score')
            plot_A[i, 0].legend(loc='upper left', title='Linker state', handles=handles, labels=labs)
            plot_A[i, 1].xaxis.set_label_position('top')
            plot_A[i, 1].set_xlabel('linker 1 score')
            plot_A[i, 1].axes.get_yaxis().set_visible(False)
            plot_A[i, 1].set(ylabel=None)
        if i == 3:
            plot_A[i, 0].set_ylabel('{}\nCount'.format(name[i]))
            plot_A[i, 0].legend(loc='upper left', title='Linker state', handles=handles, labels=labs)
            plot_A[i, 1].set(ylabel=None)
            plot_A[i, 1].axes.get_yaxis().set_visible(False)
        else:
            plot_A[i, 1].set(ylabel=None)
            plot_A[i, 1].axes.get_yaxis().set_visible(False)
            plot_A[i, 1].axes.get_xaxis().set_visible(False)

            plot_A[i, 0].set_ylabel('{}\nCount'.format(name[i]))
            plot_A[i, 0].set(xlabel=None)
            plot_A[i, 0].axes.get_xaxis().set_visible(False)
            plot_A[i, 0].legend(loc='upper left', title='Linker state', handles=handles, labels=labs)


    subfignest = subfig[1].subfigures(2, 1)
    plot_B = subfignest[0].subplots(2, 2, sharex=True)
    subfignest[0].suptitle('B', x=0.05, ha='left', fontweight='bold', fontsize=18)

    df = pd.read_csv('/media/draco/lucask/linkercomparison/Linker_comparison.csv', sep='\t')
    #df = pd.read_csv('/home/lucas/Documents/Work/2022/11-2022/Linker_comparison.csv', sep='\t')
    df['Count'] = df['Count'].astype(int)
    df['Percentage'] = (df['Count']/df.groupby('Dataset')['Count'].transform('sum')) * 100

    sns.barplot(ax=plot_B[0, 0], data=df.loc[df.Type == 'Correctable Linker'], x='Dataset', y='Percentage', palette='Set1', ci=None, edgecolor='black')
    plot_B[0, 0].set_ylabel('Percentage')
    plot_B[0, 0].set(xlabel=None)
    plot_B[0, 0].set_title('Correctable Linkers')
    sns.barplot(ax=plot_B[0, 1], data=df.loc[df.Type == 'Uncorrectable Linkers'], x='Dataset', y='Percentage', palette='Set1', ci=None, edgecolor='black')
    plot_B[0, 1].set_title('Uncorrectable Linkers')
    plot_B[0, 1].set(xlabel=None)
    sns.barplot(ax=plot_B[1, 0], data=df.loc[df.Type == 'Uncorrectable Linker 2'], x='Dataset', y='Percentage', palette='Set1', ci=None, edgecolor='black')
    plot_B[1, 0].set_ylabel('Percentage')
    plot_B[1, 0].set(xlabel=None)
    plot_B[1, 0].set_title('Uncorrectable Linker 2')
    sns.barplot(ax=plot_B[1, 1], data=df.loc[df.Type == 'Uncorrectable Linker 1'], x='Dataset', y='Percentage', palette='Set1', ci=None, edgecolor='black')
    plot_B[1, 1].set_title('Uncorrectable Linker 1')
    plot_B[1, 1].set(xlabel=None)


    df = pd.read_csv('/media/draco/lucask/linkercomparison/Linker_Dataset_variance.csv', sep='\t')
    #df = pd.read_csv('/home/lucas/Documents/Work/2022/11-2022/Linker_Dataset_variance.csv', sep='\t')
    df = df.dropna()
    df['Count'] = df['Count'].astype(int)
    df['Percentage'] = df['Count']/df.groupby('Dataset')['Count'].transform('sum')

    plot_C = subfignest[1].subplots(1, 1)
    subfignest[1].suptitle('C', x=0.05, ha='left', fontweight='bold', fontsize=18)
    sns.barplot(ax=plot_C, data=df, x='Type', y='Percentage', hue='Linker', ci="sd", palette='Set2', edgecolor='black')
    plot_C.set(xlabel=None)
    plot_C.legend(title=None)
    plot_C.set_title('Experiment Variance')

    plt.savefig("/media/draco/lucask/linkercomparison/figure_X.png")
    #plt.savefig("/home/lucas/Documents/Work/2022/11-2022/fig_x.png")

if __name__ == '__main__': main()