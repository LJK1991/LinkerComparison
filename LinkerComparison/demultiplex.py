import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandarallel import pandarallel
from Bio import pairwise2
import time
import math
import sys
import argparse
import os

from utils import *
from plotting import *
'''
This code is heavily adapted from https://github.com/fairliereese/LR-splitpipe
Made the code run in more compartementilized blocks, this way it was easier to check for any errors that could/would occur.
'''

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mode', dest='mode', help='Which functions you want to run. options: fastq_to_df, score_linkers, align_linkers, align_linker_halves, get_lev_dist, plot')
	parser.add_argument('-f', '--infile', dest='fastq', help='FASTQ file you want to analyze containing SPLiT-seq reads')
	parser.add_argument('-o', '--output', dest='output', help='Output file path/prefix')
	parser.add_argument('--chunksize', dest='chunksize', default=10**5, help='Number of lines to read in at any given time')
	parser.add_argument('--verbosity', dest='verbosity', default=1, help='Verbosity setting. Higher number = more messages')
	parser.add_argument('-t','--threads', dest='threads', default=4, help='The amount of core that should be used during parallelizing. default is 4')
	parser.add_argument('--delete_input', dest='delete_input', action='store_true', help='Delete temporary files', default=False)
	parser.add_argument('-l', '--linker', dest='linker', default='split', help='which linker sequences should be used as reference, default is split. options: split, GtoT, 2as1, share, parse')
	parser.add_argument('--mm', dest='mm', default=0.1, help='Percent of mismatches allowed in linker sequences (based on length of reference sequences)')
	parser.add_argument('-p', '--plot_type', dest='plot', default='all', help='Which plots to make. default = all. options:\n'
																 'score - plot the linker scores\n'
																 'heatmap - plots the amount of reads its find per score in a heatmap\n'
																 'link_rel - includes the relation between linkers within the score plot\n'
																 'link_align - plots linkers against referenec sequence and shows percentage of correct bases per position\n'
																 'all - creates all plots')
	parser.add_argument('-e', '--edit_dist', dest='edit_dist', default=2, help='the maximum levenhstein distance with which barcodes can be corrected.')
	args = parser.parse_args()
	return args

def format_chunksize(c):
	if '**' in str(c):
		i, j = c.split('**')
		i = float(i)
		j = float(j)
		c = i**j
	c = int(c)
	return c

def main():
	args = get_args()
	mode = args.mode
	oprefix = args.output
	t = int(args.threads)
	v = int(args.verbosity)
	delete_input = args.delete_input
	chunksize = format_chunksize(args.chunksize)
	linker = args.linker
	fname = args.fastq
	mm = float(args.mm)
	plot = args.plot
	bc_edit_dist = args.edit_dist

	if mode == 'fastq_to_df' or mode == 'all':
		fname = fastq_to_df(fname, oprefix, chunksize, v)

	if mode == 'score_linkers' or mode == 'all':
		if fname.endswith('.fastq'):
			print('file is a fastq file. please start with the fastq_to_df function.')
		fname = score_linkers(fname, oprefix, t=t, chunksize=chunksize, verbose=v, delete_input=False, linker=linker)

	if mode == 'align_linkers' or mode == 'all':
		if fname == None:
			fname = str(oprefix) + '_table.tsv'
		fname = align_linkers(fname, oprefix, t, chunksize, mm, v, delete_input, linker)

	if mode == 'align_linker_halves' or mode == 'all':
		if fname == None:
			fname = str(oprefix) + '_seq_linker_alignments.tsv'
		fname = align_halves(fname, oprefix, t, v, chunksize, linker)

	if mode == 'get_lev_dist' or mode == 'all':
		if fname == None:
			fname = str(oprefix) + '_seq_linker_alignments.tsv'
		fname = get_lev_dist(fname, oprefix, bc_edit_dist, chunksize, t, v)

	#change the plots that can be made
	#and test
	if mode == 'plot' or mode == 'all':
		plotting_graphs(oprefix, linker, plot, v, mm)

if __name__ == '__main__': main()