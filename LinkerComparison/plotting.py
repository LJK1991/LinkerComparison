import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import seaborn as sns
import pandas as pd
import math
import numpy as np
from cycler import cycler

mplstyle.use('fast')

def plot_hist(x, **kwargs):
	ax = sns.histplot(x, binwidth=1)

def plot_mm_lines(x, l1_len, l2_len, mm,**kwargs):
	col = x.unique().tolist()[0]
	# linker 1
	if 'l1' in col:
		mismatch_lines = list(range(l1_len-math.ceil(mm*l1_len), l1_len+1))
	# linker 2
	elif 'l2' in col:
		mismatch_lines = list(range(l2_len-math.ceil(mm*l2_len), l2_len+1))
	for l in mismatch_lines:
		plt.axvline(l, color='k', linestyle='-', linewidth=1)

def plot_linker_scores(df, oprefix, l1_len, l2_len, mm):

	val_vars = ['l2_score', 'l1_score']
	cols = ['index'] + val_vars
	temp = df[cols].melt(id_vars='index', value_vars=val_vars)

	g = sns.FacetGrid(temp, col='variable')
	g.map(plot_hist, 'value')
	g.map(plot_mm_lines, 'variable', l1_len=l1_len, l2_len=l2_len, mm=mm)

	fname = oprefix+'_linker_score_dists.png'
	plt.savefig(fname)

	plt.clf()

# plot heatmap of number of reads recovered with different linker
# mismatch allowances
def plot_linker_heatmap(df, oprefix, l1_len, l2_len, mm):
	# proportional allowances
	p = [0, 1, 2, 3, 4, 5]
	data = [[0 for i in range(len(p))] for j in range(len(p))]
	p_df = pd.DataFrame(data, index=p, columns=p)
	for i in p: # l1
		for j in p: # l2
			l1_min = l1_len - i
			l2_min = l2_len - j
			fwd_df = df.loc[(df.l1_score>=l1_min)&(df.l2_score>=l2_min)]
			one_dir_reads = df.loc[list(set(fwd_df.index))]
			p_df.at[i, j] = len(one_dir_reads.index)

	ax = sns.heatmap(p_df, annot=True)
	_ = ax.set(xlabel='Percent mismatches/indels allowed in l2', ylabel='Percent mismatches/indels allowed in l1', title='Reads recovered')

	fname = '{}_proportion_linker_mismatch_heatmap.png'.format(oprefix)
	fig = plt.gcf()
	fig.set_size_inches(8, 8)
	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def plot_linker_relation(df, oprefix, verbosity=1):
	if verbosity > 1:
		print('Start plotting')

	#colpal = {1:'green', 2:'yellow', 3:'blue', 4:'red'}
	# 1: Correctable linker 2: Uncorrectable linker 1 3: Uncorrectable linker 2 4: Uncorrectable linker
	#colpal = {"Correctable linker": 'green', "Uncorrectable linker 1": 'yellow', "Uncorrectable linker 2": 'blue', "Uncorrectable linker": 'red'}
	labs = ["Correctable linker", "Uncorrectable linker 1", "Uncorrectable linker 2", "Uncorrectable linker"]
	#labs = {1:"Correctable linker", 2:"Uncorrectable linker 1", 3:"Uncorrectable linker 2", 4:"Uncorrectable linker"}
	fig, axes = plt.subplots(1, 2, figsize=(10,5), sharey=True, layout='tight')
	sns.histplot(ax=axes[1], data=df, x='l1_score', multiple='stack', hue='linker_relation', hue_order=[1,2,3,4], legend=False, element='bars', binwidth=1, binrange=(0,30), palette='tab10')
	sns.histplot(ax=axes[0], data=df, x='l2_score', multiple='stack', hue='linker_relation', hue_order=[1,2,3,4], legend=True, element='bars', binwidth=1, binrange=(0,30), palette='tab10')
	legend = axes[0].get_legend()
	handles = legend.legendHandles
	legend.remove()
	axes[0].legend(loc='upper left', title='Linker state', handles=handles, labels=labs)
	#axes[1].legend(loc='center right', bbox_to_anchor=(1.25,0.25,0.5,0.5), title='Linker state', handles=handles, labels=labs)

	fname = oprefix+'_linker_rel_dists.png'
	plt.savefig(fname)

	if verbosity > 1:
		print('Finished plotting ')
	plt.clf()

	rel_counts = df['linker_relation'].value_counts()
	fname = oprefix + '_linker_rel_counts.tsv'
	rel_counts.to_csv(fname, sep='\t')

def linker_analysis(df, oprefix, l1, l2, l1_len, l2_len, verbosity=1):
	colpal = {1: 'green', 2: 'yellow', 3: 'blue', 4: 'red'}
	# 1: Correctable linker 2: Uncorrectable linker 1 3: Uncorrectable linker 2 4: Uncorrectable linker

	df['l1_start_dif'] = np.nan
	df['l1_stop_dif'] = np.nan
	df['l2_start_dif'] = np.nan
	df['l2_stop_dif'] = np.nan
	df['l1_pb'] = np.nan
	df['l2_pb'] = np.nan
	df['l1_pb'] = df['l1_pb'].astype(object)
	df['l2_pb'] = df['l2_pb'].astype(object)

	if verbosity > 1:
		print('starting per base calculation')

	for i in df.index:
		l1_start = df.l1_start[i]
		#l1_stop = df.l1_stop[i]
		l2_start = df.l2_start[i]
		#l2_stop = df.l2_stop[i]

		l1_seq = df.seq[i][l1_start:l1_start + l1_len]
		l2_seq = df.seq[i][l2_start:l2_start + l2_len]

		score = []
		# if shorter | starts before or after certain position disregard read
		#if len(l1_seq) < 15 or l1_start < (l1_start - 8) or l1_start > (l1_start + 8):
		#	same_base = np.nan
		#else:
		#for base in reference linker
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

		# if shorter | starts before or after certain position disregard read
		#if len(l2_seq) < 15 or l2_start < (l2_start - 8) or l2_start > (l2_start + 8):
		#	same_base = np.nan
		#else:
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

	if verbosity > 1:
		print('finished per base calculation, start plotting')

	labs = ["Correctable\nlinker", "Uncorrectable\nlinker 1", "Uncorrectable\nlinker 2", "Uncorrectable\nlinker"]
	plt.rc('axes', prop_cycle=cycler('color',['green', 'yellow', 'blue', 'red']))
	#here extract all score's per catagory and make 1 big plot.
	fig = plt.figure(figsize=(12,8), layout='tight')
	#subfigs = fig.subfigures(3,1, height_ratios=[1,1,0.2])
	ax1 = fig.subplots(1,2, sharey=True)
	#ax2 = subfigs[1].subplots(1,6)
	rel = [1,2,3,4]

	for state in rel:
		tmp_df_l1 = pd.DataFrame.from_dict(dict(df.loc[df['linker_relation'] == state, 'l1_pb'])).T
		l1_prop = tmp_df_l1.sum(axis=0)/len(tmp_df_l1.index)
		sns.lineplot(ax=ax1[1], data=l1_prop, legend=False)

		tmp_df_l2 = pd.DataFrame.from_dict(dict(df.loc[df['linker_relation'] == state, 'l2_pb'])).T
		l2_prop = tmp_df_l2.sum(axis=0) / len(tmp_df_l2.index)
		sns.lineplot(ax=ax1[0], data=l2_prop, legend=True)

	ax1[0].set_title('linker 2')
	ax1[0].set(xlabel=None)
	ax1[0].set_ylabel('proportion')
	ax1[0].set_xlabel('base position')
	ax1[1].set_title('linker 1')
	ax1[1].set(xlabel=None, ylabel=None)
	ax1[1].get_yaxis().set_visible(False)
	ax1[1].set_xlabel('base position')
	ax1[1].legend(loc='center right', bbox_to_anchor=(0.9, 0.25, 0.5, 0.5), labels=labs, title='Linker state')

	fname = oprefix + '_linker_perbase_analysis.png'
	plt.savefig(fname)
	plt.clf()

	fig = plt.figure(figsize=(12,8) ,layout='constrained')
	subfig = fig.subfigures(1,2, wspace=0.0, width_ratios=[1,2])
	ax2 = subfig[0].subplots(1,1)
	ax3 = subfig[1].subplots(1,2, sharey=True)
	#ax2 = fig.subplots(1,3)
	#sns.histplot(ax=ax2[0], x=df['l2_len'], hue=df['linker_relation'], stat='percent', multiple='layer', element='poly', binwidth=1)
	sns.violinplot(ax=ax2, x=df['linker_relation'], y=df['l2_len'], label=labs, palette=colpal, scale='width')
	sns.violinplot(ax=ax3[0], x=df['linker_relation'], y=df['l2_start'], label=labs, palette=colpal, scale='width')
	sns.violinplot(ax=ax3[1], x=df['linker_relation'], y=df['l2_stop'], label=labs, palette=colpal, scale='width')

	ax2.set_title('linker 2 length')
	ax2.tick_params(labelrotation=45)
	ax2.set(ylabel=None, xlabel=None)
	ax2.set_xticklabels(labs)
	ax3[0].set_title('linker 2 start position')
	ax3[0].set(ylabel=None, xlabel=None)
	ax3[0].tick_params(labelrotation=45)
	ax3[0].set_yticks([0,18,48,56,86,100])
	ax3[0].set_xticklabels(labs)
	ax3[1].set_title('linker 2 stop position')
	ax3[1].set(ylabel=None, xlabel=None)
	ax3[1].get_yaxis().set_visible(False)
	ax3[1].tick_params(labelrotation=45)
	ax3[1].set_xticklabels(labs)

	fname = oprefix + '_linker2_stat_analysis.png'
	plt.savefig(fname)
	plt.clf()

	fig = plt.figure(figsize=(12,8) ,layout='constrained')
	subfig = fig.subfigures(1, 2, wspace=0.0, width_ratios=[1, 2])
	ax4 = subfig[0].subplots(1, 1)
	ax5 = subfig[1].subplots(1, 2, sharey=True)

	sns.violinplot(ax=ax4, x=df['linker_relation'], y=df['l1_len'], label=labs, palette=colpal, scale='width')
	sns.violinplot(ax=ax5[0], x=df['linker_relation'], y=df['l1_start'], label=labs, palette=colpal, scale='width')
	sns.violinplot(ax=ax5[1], x=df['linker_relation'], y=df['l1_stop'], label=labs, palette=colpal, scale='width')

	ax4.set_title('linker 1 length')
	ax4.set(ylabel=None, xlabel=None)
	ax4.tick_params(labelrotation=45)
	ax4.set_xticklabels(labs)
	ax5[0].set_title('linker 1 start position')
	ax5[0].set(ylabel=None, xlabel=None)
	ax5[0].tick_params(labelrotation=45)
	ax5[0].set_yticks([0, 18, 48, 56, 86, 100])
	ax5[0].set_xticklabels(labs)
	ax5[1].set_title('linker 1 stop position')
	ax5[1].set(ylabel=None, xlabel=None)
	ax5[1].get_yaxis().set_visible(False)
	ax5[1].tick_params(labelrotation=45)
	ax5[1].set_xticklabels(labs)

	fname = oprefix + '_linker1_stat_analysis.png'
	plt.savefig(fname)
	plt.clf()
	if verbosity > 1:
		print('Finished plotting')




