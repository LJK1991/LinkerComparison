import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype
import re
import numpy as np
from cycler import cycler
import matplotlib.colors as mcolors

from utils import *
from plotting import *

def main():
    flist = ['/media/draco/lucask/linkercomparison/SPLIT/SPLiT',
             '/media/draco/lucask/linkercomparison/GtiT/GtoT',
             '/media/draco/lucask/linkercomparison/SHARE/SHARE',
             '/media/draco/lucask/linkercomparison/parse/parse']

    tlist = ['SPliT','G to T', 'SHARE', 'Parse']

    pal = list(mcolors.TABLEAU_COLORS.values())


    fig = plt.figure(figsize=(24, 16), constrained_layout=True)
    subfig = fig.subfigures(1, 2, width_ratios=[1, 1])
    subfignest = subfig[1].subfigures(2, 1, height_ratios=[1, 2])

    ##########################################
    ### plot_A - per base per category    ####
    ##########################################
    plot_A = subfig[0].subfigures(4, 1, wspace=0)
    subfig[0].suptitle('A', x=0.05, ha='left', fontweight='bold', fontsize=18)
    a = 0
    for file in flist:
        file = file + '_seq_linker_halves_alignments.tsv'
        df = pd.read_csv(file, sep='\t', index_col=None, usecols=[0, 3, 4, 5, 6, 7, 8, 9, 10, 17])

        try:
            if re.search('parse', file, re.IGNORECASE):
                l1, l2 = get_linkers(linker_type='parse')
                l1_len = len(l1)
                l2_len = len(l2)
            elif re.search('share', file, re.IGNORECASE):
                l1, l2 = get_linkers(linker_type='share')
                l1_len = len(l1)
                l2_len = len(l2)
            elif re.search('split', file, re.IGNORECASE):
                l1, l2 = get_linkers(linker_type='split')
                l1_len = len(l1)
                l2_len = len(l2)
            elif re.search('gtot', file, re.IGNORECASE):
                l1, l2 = get_linkers(linker_type='GtoT')
                l1_len = len(l1)
                l2_len = len(l2)
            else:
                l1 = l2 = np.nan
                l1_len = len(l1)
                l2_len = len(l2)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))

        df['l1_pb'] = np.nan
        df['l2_pb'] = np.nan
        df['l1_pb'] = df['l1_pb'].astype(object)
        df['l2_pb'] = df['l2_pb'].astype(object)

        for i in df.index:
            l1_start = df.l1_start[i]
            l2_start = df.l2_start[i]

            l1_seq = df.seq[i][l1_start:l1_start + l1_len]
            l2_seq = df.seq[i][l2_start:l2_start + l2_len]

            score = []

            # for base in reference linker
            for j in range(0, l1_len):
                # if j bigger then length of extracted linker, base is automatically considerd bad
                if j >= len(l1_seq):
                    same_base = False
                # if j is shorter than extracted linker, compare it to refrence
                elif j < len(l1_seq):
                    same_base = l1_seq[j] == l1[j]
                # add to score
                score.append(same_base)

            # append to dataframe
            df.at[i, 'l1_pb'] = score

            score = []
            # for base in reference linker
            for j in range(0, l2_len):
                # if j bigger then length of extracted linker, base is automatically considerd bad
                if j >= len(l2_seq):
                    same_base = False
                # if j is shorter than extracted linker, compare it to refrence
                elif j < len(l2_seq):
                    same_base = l2_seq[j] == l2[j]
                # add to score
                score.append(same_base)

            # append to dataframe
            df.at[i, 'l2_pb'] = score

        #start plotting
        tmp_ax = plot_A[a].subplots(1, 2, sharey=True, gridspec_kw={'wspace':0})
        labs = ["Correctable\nlinker", "Uncorrectable\nlinker 1", "Uncorrectable\nlinker 2", "Uncorrectable\nlinker"]

        rel = [1, 2, 3, 4]
        for state in rel:
            tmp_df_l1 = pd.DataFrame.from_dict(dict(df.loc[df['linker_relation'] == state, 'l1_pb'])).T
            l1_prop = tmp_df_l1.sum(axis=0) / len(tmp_df_l1.index)
            sns.lineplot(ax=tmp_ax[1], data=l1_prop, legend=False, palette='tab10')

            tmp_df_l2 = pd.DataFrame.from_dict(dict(df.loc[df['linker_relation'] == state, 'l2_pb'])).T
            l2_prop = tmp_df_l2.sum(axis=0) / len(tmp_df_l2.index)
            sns.lineplot(ax=tmp_ax[0], data=l2_prop, legend=True, palette='tab10')

        if a == 0:
            tmp_ax[0].set_title('linker 2')
            tmp_ax[1].set_title('linker 1')
        else:
            tmp_ax[0].set_title(None)
            tmp_ax[1].set_title(None)
        if a == 3:
            tmp_ax[1].set_xlabel('base position')
            tmp_ax[0].set_xlabel('base position')
        else:
            tmp_ax[0].set_xlabel(None)
            tmp_ax[1].set_xlabel(None)

        tmp_ax[0].set_ylabel('proportion')
        tmp_ax[1].get_yaxis().set_visible(False)
        tmp_ax[1].legend(loc='center right', bbox_to_anchor=(0.9, 0.25, 0.5, 0.5), labels=labs, title='{}'.format(tlist[a]))
        tmp_ax[0].axvline(x=14.5, color='black', label='Site of ligation', ls=':', lw='2')
        tmp_ax[1].axvline(x=14.5, color='black', label='Site of ligation', ls=':', lw='2')
        a += 1


    ##########################################
    ### table_S2 - Lev_dist per category  ####
    ##########################################
    #Lev_dist per category fucking impossible to plot!

    for file in flist:
        file = file + '_corrected_bcs.tsv'
        df = pd.read_csv(file, sep='\t', index_col=None, usecols=[1, 17, 18, 20, 22])

        BC = [1, 2, 3]
        new_df = pd.DataFrame()
        for bc in BC:
            rname = df.read_name
            link_rel = df.linker_relation
            barcode = ['BC' + str(bc)] * len(df.index)
            lev_dist = df.iloc[:, bc]
            tmp_df = pd.DataFrame(list(zip(rname, link_rel, barcode, lev_dist)), columns=['read_name','linker_relation', 'barcode', 'lev_dist'])
            new_df = pd.concat([new_df, tmp_df])

        linker_type = CategoricalDtype(categories=[1,2,3,4], ordered=True)
        new_df.linker_relation = new_df.linker_relation.astype(linker_type)
        new_df.linker_relation = new_df.linker_relation.cat.rename_categories({1:'Correctable\nlinkers', 2:'Uncorrectable\nlinker 1', 3:'Uncorrectable\nlinker 2', 4:'Uncorrectable\nlinkers'})

        bc_type = CategoricalDtype(categories=['BC3','BC2','BC1'], ordered=True)
        new_df.barcode = new_df.barcode.astype(bc_type)

        lev_type = CategoricalDtype(categories=[0,1,2,4], ordered=True)
        new_df.lev_dist = new_df.lev_dist.astype(lev_type)
        new_df.lev_dist = new_df.lev_dist.cat.rename_categories({0:'0', 1:'1', 2:'2', 4:'>2'})

        x_df = new_df.groupby(by=['linker_relation', 'barcode'])['lev_dist'].value_counts(normalize=True).unstack()
        file = re.sub('_corrected_bcs.tsv','_lev_dist_per_bc.tsv',file)
        x_df.to_csv(file, sep='\t', index=False)

    ##########################################
    ### plot_B - BC pct Dotplot           ####
    ##########################################
    #barcode Dotplot
    plot_B = subfignest[0].subplots(1, 5, sharey=True, gridspec_kw={'wspace':0, 'width_ratios':[1,1,1,1,0.25]}, )
    subfignest[0].suptitle('B', x=0.05, ha='left', fontweight='bold', fontsize=18)
    plot_B[4].axis('off')

    b = 0
    for file in flist:
        file = file + '_corrected_bcs.tsv'
        df = pd.read_csv(file, sep='\t', index_col=None, usecols=[1, 17, 24, 25, 26])
        BC = [1, 2, 3]
        new_df = pd.DataFrame()
        for bc in BC:
            rname = df.read_name
            link_rel = df.linker_relation
            barcode = ['BC' + str(bc)] * len(df.index)
            is_expected = df.iloc[:, (1 + bc)]
            tmp_df = pd.DataFrame(list(zip(rname, link_rel, barcode, is_expected)), columns=['read_name', 'linker_relation', 'barcode', 'is_expected_bc'])
            new_df = pd.concat([new_df, tmp_df])

        linker_type = CategoricalDtype(categories=[1,2,3,4], ordered=True)
        new_df.linker_relation = new_df.linker_relation.astype(linker_type)
        new_df.linker_relation = new_df.linker_relation.cat.rename_categories({1:'Correctable\nlinkers', 2:'Uncorrectable\nlinker 1', 3:'Uncorrectable\nlinker 2', 4:'Uncorrectable\nlinkers'})

        bc_type = CategoricalDtype(categories=['BC3','BC2','BC1'], ordered=True)
        new_df.barcode = new_df.barcode.astype(bc_type)

        new_df.is_expected_bc = new_df.is_expected_bc.astype("bool")

        x_df = new_df.groupby(by=['linker_relation', 'barcode']).agg({'is_expected_bc':'sum'})
        y_df = new_df.value_counts(subset=['linker_relation', 'barcode']).sort_index()
        new_df_pct = pd.DataFrame(100 * (x_df.is_expected_bc/y_df.values))

        sns.scatterplot(data=new_df_pct, x="barcode", y="linker_relation", size="is_expected_bc", hue="is_expected_bc", sizes=(10, 750), palette='RdYlGn', ax=plot_B[b])
        plot_B[b].margins(x=0.5, y=0.25)
        if b == 3:
            plot_B[b].get_yaxis().set_visible(False)
            plot_B[b].spines['top'].set_visible(False)
            plot_B[b].spines['right'].set_visible(False)
            plot_B[b].spines['left'].set_visible(False)
            plot_B[b].spines['bottom'].set_visible(False)
            plot_B[b].legend(bbox_to_anchor=(1, .6), title='Percent\nTrue BC', borderpad=1.2)
            plot_B[b].set(xlabel=None)
            plot_B[b].set_title('{}'.format(tlist[b]))
        else:
            plot_B[b].spines['top'].set_visible(False)
            plot_B[b].spines['right'].set_visible(False)
            plot_B[b].spines['bottom'].set_visible(False)
            plot_B[b].spines['left'].set_visible(False)
            #plot_B[b].spines['left'].set_position(('axes', 0.25))
            plot_B[b].set(xlabel=None, ylabel=None)
            plot_B[b].get_legend().remove()
            plot_B[b].set_title('{}'.format(tlist[b]))

        if b > 0:
            plot_B[b].tick_params(axis='y', which='both', left=False, right=False)
        b += 1

        ############################################
        ### plot_C - linker halve scores Violin  ###
        ############################################
        plot_C = subfignest[1].subfigures(4, 1, wspace=0)
        subfignest[1].suptitle('C', x=0.05, ha='left', fontweight='bold', fontsize=18)

        x = 0
        for file in flist:
            file = file + '_seq_linker_halves_alignments.tsv'
            df = pd.read_csv(file, sep='\t', index_col=None, usecols=[17,18,22,27,31])

            cols = ['l1_half1_score','l1_half2_score','l2_half1_score','l2_half2_score']
            new_df = pd.DataFrame()
            for i in range(1,len(cols)+1):
                link_rel = df.linker_relation
                linker = [math.ceil(i/2)] * len(df.index)
                if (i % 2) == 1:
                    half = 1
                else:
                    half = 2
                link_half = [half] * len(df.index)
                score = df[cols[i-1]]
                tmp_df = pd.DataFrame(list(zip(link_rel, linker, link_half, score)), columns=['linker_relation', 'linker', 'linker_half', 'score'])
                new_df = pd.concat([new_df, tmp_df])

            linker_type = CategoricalDtype(categories=[1, 2, 3, 4], ordered=True)
            new_df.linker_relation = new_df.linker_relation.astype(linker_type)
            new_df.linker_relation = new_df.linker_relation.cat.rename_categories(
                {1: 'Correctable\nlinkers', 2: 'Uncorrectable\nlinker 1', 3: 'Uncorrectable\nlinker 2',
                 4: 'Uncorrectable\nlinkers'})

            link_type = CategoricalDtype(categories=[2, 1], ordered=True)
            new_df.linker = new_df.linker.astype(link_type)

            half_type = CategoricalDtype(categories=[1, 2], ordered=True)
            new_df.linker_half = new_df.linker_half.astype(half_type)
            new_df.linker_half = new_df.linker_half.cat.rename_categories({1: 'first half', 2: 'second half'})

            new_df.reset_index(inplace=True)

            ax = plot_C[x].subplots(1, 4, sharex=True, sharey=True, gridspec_kw={'wspace':0, 'hspace':0})
            labs = ["Correctable\nlinkers", "Uncorrectable\nlinker 1", "Uncorrectable\nlinker 2", "Uncorrectable\nlinkers"]
            h = []
            l = ['First half', 'Second half']
            for i in range(0, len(ax)):
                sns.violinplot(data=new_df[new_df.linker_relation == labs[i]], x='linker', y='score', hue='linker_half',
                   color=pal[i], split=True, ax=ax[i])

                handles, labels = ax[i].get_legend_handles_labels()
                h.extend(handles)
                l.extend(['',''])

                if x == 0:
                    ax[i].set_title(labs[i])
                else:
                    ax[i].set_title(None)
                if i == 0:
                    ax[i].spines['right'].set_visible(False)
                    ax[i].set(ylabel='{}\n Alignment score'.format(tlist[x]))
                elif i == 3:
                    ax[i].get_yaxis().set_visible(False)
                    ax[i].spines['left'].set_visible(False)
                    ax[i].get_legend().remove()
                else:
                    ax[i].get_yaxis().set_visible(False)
                    ax[i].spines['right'].set_visible(False)
                    ax[i].spines['left'].set_visible(False)
                    ax[i].get_legend().remove()
                if x < 3:
                    ax[i].get_xaxis().set_visible(False)

                ax[0].legend(h, l, loc='lower left', title=None, ncol=4, markerfirst=False, columnspacing=0)

            x+=1

    plt.savefig("/media/draco/lucask/linkercomparison/Figure_y.png")

if __name__ == '__main__': main()
