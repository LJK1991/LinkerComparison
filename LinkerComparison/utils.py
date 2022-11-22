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
import os.path
from os import path

from plotting import *

###################################################################################
############################### Helper functions ##################################
###################################################################################

def get_linkers(linker_type, verbose=0):
    """
    Get the SPLiT-seq linker sequences
    """

    linker = linker_type
    # linker sequence between barcodes 1 and 2
    if linker == 'split':
        l1 = 'ATCCACGTGCTTGAGAGGCCAGAGCATTCG'
    elif linker == 'GtoT':
        l1 = 'ATCCACGTGCTTGATAGGCCAGAGCATTCG'
    elif linker == '2as1':
        l1 = 'ACGACTAACGTGATGTGGCCGATGTTTCGCATCGGCGTACGACTAACGTGAT'
        l2 = 'AACGTGATGTGGCCGATGTTTCGCATCGGCGTACGACTAACGTGATGTGGCC'
    elif linker == 'parse':
        l1 = 'ATCCACGTGCTTGAGACTGTGG'
    elif linker == 'share':
        l1 = 'ATCCACGTGCTTGAGCGCGCTGCATACTTG'
    # linker sequence between barcodes 2 and 3
    if linker != '2as1':
        l2 = 'GTGGCCGATGTTTCGCATCGGCGTACGACT'

    if verbose > 1:
        print('linker type used: ' + str(linker))
        print('linker sequence 1: ' + str(l1))
        print('linker sequence 2: ' + str(l2))

    return l1, l2

def load_barcodes_set():
    """
	Load the barcodes. Adapted from the Parse biosciences pipeline.
	"""
    pkg_path = os.path.dirname(__file__)
    pkg_path = '/'.join(pkg_path.split('/')[:-1])

    with open(pkg_path + '/barcodes/bc_dict_v1.pkl', 'rb') as f:
        edit_dict_v1 = pickle.load(f)

    # read in barcode sequences
    bc_8nt_v1 = pd.read_csv(pkg_path + '/barcodes/bc_8nt_v1.csv',names=['barcode'],index_col=0).barcode.values
    bc1_edit_dict = edit_dict_v1
    bc2_edit_dict = edit_dict_v1
    bc3_edit_dict = edit_dict_v1

    bc_8nt_set_dict = {}
    bc_8nt_set_dict['bc1'] = set(bc_8nt_v1)
    bc_8nt_set_dict['bc2'] = set(bc_8nt_v1)
    bc_8nt_set_dict['bc3'] = set(bc_8nt_v1)

    bc_8nt_bc1 = bc_8nt_set_dict['bc1']
    bc_8nt_bc2 = bc_8nt_set_dict['bc2']
    bc_8nt_bc3 = bc_8nt_set_dict['bc3']
    return list(bc_8nt_bc1), list(bc_8nt_bc2), list(bc_8nt_bc3)

def load_barcodes():
    """
	Load the barcodes. Adapted from the Parse biosciences pipeline.

	Returns:
		bc#_edit_dict (dict): Dict for barcode<1,2,3> with
			key: query bc
			item: corrected bc
	"""
    pkg_path = os.path.dirname(__file__)
    pkg_path = '/'.join(pkg_path.split('/')[:-1])
    with open(pkg_path + '/Linker_comparison/barcodes/bc_dict_v1.pkl', 'rb') as f:
        edit_dict_v1 = pickle.load(f)

    bc1_edit_dict = edit_dict_v1
    bc2_edit_dict = edit_dict_v1
    bc3_edit_dict = edit_dict_v1

    return bc1_edit_dict, bc2_edit_dict, bc3_edit_dict

def get_min_edit_dists(bc, edit_dict, max_d=3):
    """
	Returns a list of nearest edit dist seqs.
	Adapted from Parse biosciences

	Parameters:
		bc (str): 8nt barcode
		edit_dict (dict): Dict of bcs within
			edit distance of bc
		max_d (int): Edit distance

	Returns:
		bc_matches
		edit_dist
	"""
    bc_matches = edit_dict[0][bc]
    edit_dist = 0
    if (len(bc_matches)==0) and (max_d>=1):
        edit_dist+=1
        bc_matches = edit_dict[1][bc]
    if (len(bc_matches)==0) and (max_d>=2):
        edit_dist+=1
        bc_matches = edit_dict[2][bc]
    if (len(bc_matches)==0) and (max_d>=3):
        edit_dist+=1
        bc_matches = edit_dict[3][bc]
    return bc_matches, edit_dist

###################################################################################
############################### Lambda functions ##################################
###################################################################################

def score_linkers_x(x, l_seqs, l_prefs):
    """
    Function to apply across rows of a dataframe

    Parameters:
		x (pandas Series): Row from parent dataframe
		l_seqs (list of str): Linker sequences
		l_prefs (list of str): Names of linkers (l1, l2)

	Returns:
		entry (pandas Series): Row with linker alignment scores
	"""

    entry = {}

    # +1 for match
    # -1 for mismatch
    # -1 for gap open
    # -1 for gap extend
    # this scoring schema allows us to know exactly
    # how many errors are in each best alignment
    for l_seq, pref in zip(l_seqs, l_prefs):
        a = pairwise2.align.localms(x.seq, l_seq, 1, -1, -1, -1, one_alignment_only=True, score_only=True)
        score = a
        entry['{}_score'.format(pref)] = score
    return entry

def align_linkers_x(x, l1, l2):
    """
    Function to find inds of linkers across rows of a dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe
		l1 (str): Linker 1
		l2 (str): Linker 2

	Returns:
		entry (pandas Series): Row with linker alignment indices
	"""
    entry = {}

    # compute alignments for both linkers in the correct orientation
    # use the default alignment
    l1_a = pairwise2.align.localms(x.seq, l1, 1, -1, -1, -1, one_alignment_only=True)
    l2_a = pairwise2.align.localms(x.seq, l2, 1, -1, -1, -1, one_alignment_only=True)

    # grab start and end of each linker
    try:
        l1_start = l1_a[0].start
    except:
        print("couldn't align reads")
        print(x.read_name)
        raise ValueError()
    l1_stop = l1_a[0].end
    l1_score = l1_a[0].score
    l2_start = l2_a[0].start
    l2_stop = l2_a[0].end
    l2_score = l2_a[0].score

    # calculate some metrics
    l1_len = l1_stop - l1_start
    l2_len = l2_stop - l2_start
    bc2_len = l1_start - l2_stop

    # construct an entry
    entry['l1_score'] = l1_score
    entry['l2_score'] = l2_score
    entry['l1_start'] = l1_start
    entry['l1_stop'] = l1_stop
    entry['l2_start'] = l2_start
    entry['l2_stop'] = l2_stop
    entry['l1_len'] = l1_len
    entry['l2_len'] = l2_len
    entry['bc2_len'] = bc2_len

    #for quicker UMI demultiplexing
    umi = x.seq[l2_start - 18:l2_start - 8]
    entry['umi'] = umi
    entry['umi_len'] = len(umi)
    bc1 = x.seq[l1_stop:l1_stop + 8]
    bc2 = x.seq[l2_stop:l2_stop + 8]
    bc3 = x.seq[l2_start - 8:l2_start]
    entry['raw_bc1'] = bc1
    entry['raw_bc2'] = bc2
    entry['raw_bc3'] = bc3
    entry['bc'] = '_'.join([bc3, bc2, bc1])

    return entry

def align_linker_halves_x(x, l1_half1, l1_half2, l2_half1, l2_half2):
    l1_start = x.l1_start
    l1_stop = x.l1_stop
    l2_start = x.l2_start
    l2_stop = x.l2_stop

    entry = {}

    l1_h1 = pairwise2.align.localms(x.seq, l1_half1, 1, -1, -1, -1, one_alignment_only=True)
    l1_h2 = pairwise2.align.localms(x.seq, l1_half2, 1, -1, -1, -1, one_alignment_only=True)
    l2_h1 = pairwise2.align.localms(x.seq, l2_half1, 1, -1, -1, -1, one_alignment_only=True)
    l2_h2 = pairwise2.align.localms(x.seq, l2_half2, 1, -1, -1, -1, one_alignment_only=True)

    try:
        l1_h1_start = l1_h1[0].start
    except:
        print("couldn't align reads")
        print(x.read_name)
        raise ValueError()
    l1_h1_stop = l1_h1[0].end
    try:
        l1_h2_start = l1_h2[0].start
    except:
        print("couldn't align reads")
        print(x.read_name)
        raise ValueError()
    l1_h2_stop = l1_h2[0].end
    try:
        l2_h1_start = l2_h1[0].start
    except:
        print("couldn't align reads")
        print(x.read_name)
        raise ValueError()
    l2_h1_stop = l2_h1[0].end
    try:
        l2_h2_start = l2_h2[0].start
    except:
        print("couldn't align reads")
        print(x.read_name)
        raise ValueError()
    l2_h2_stop = l2_h2[0].end

    l1_h1_len = l1_h1_stop - l1_h1_start
    l1_h2_len = l1_h2_stop - l1_h2_start
    l2_h1_len = l2_h1_stop - l2_h1_start
    l2_h2_len = l2_h2_stop - l2_h2_start

    entry['l1_half1_score'] = l1_h1[0].score
    entry['l1_half1_start'] = l1_h1_start
    entry['l1_half1_stop'] = l1_h1_stop
    entry['l1_half1_len'] = l1_h1_len
    entry['l1_half2_score'] = l1_h2[0].score
    entry['l1_half2_start'] = l1_h2_start
    entry['l1_half2_stop'] = l1_h2_stop
    entry['l1_half2_len'] = l1_h2_len
    entry['l1_delta_halves'] = l1_h1_stop - l1_h2_start
    entry['l2_half1_score'] = l2_h1[0].score
    entry['l2_half1_start'] = l2_h1_start
    entry['l2_half1_stop'] = l2_h1_stop
    entry['l2_half1_len'] = l2_h1_len
    entry['l2_half2_score'] = l2_h2[0].score
    entry['l2_half2_start'] = l2_h2_start
    entry['l2_half2_stop'] = l2_h2_stop
    entry['l2_half2_len'] = l2_h2_len
    entry['l2_delta_halves'] = l2_h1_stop - l2_h2_start
    entry['l1_halfVfull_start'] = bool(l1_h1_start == l1_start)
    entry['l1_halfVfull_stop'] = bool(l1_h2_stop == l1_stop)
    entry['l2_halfVfull_start'] = bool(l2_h1_start == l2_start)
    entry['l2_halfVfull_stop'] = bool(l2_h2_stop == l2_stop)

    return entry

def correct_bcs_x(x, bc_edit_dist, bc1_dict, bc2_dict, bc3_dict):
    """
	Function to correct bcs across rows of a dataframe.
	Adapted from Parse biosciences

	Parameters:
		x (pandas Series): Row from parent dataframe
		counts (pandas DataFrame): Output from get_perfect_bc_counts
			with number of reads observed for each barcode
		count_thresh (int): Minimum number of reads
		bc_edit_dist (int): Maximum edit distance of bc to correct
		bc#_edit_dict (dict): Dict for barcode<1,2,3> with
			key: query bc
			item: corrected bc

	Returns:
		entry (pandas Series): Row with corrected bcs
	"""
    bc1 = x.raw_bc1
    bc2 = x.raw_bc2
    bc3 = x.raw_bc3

    debug = True
    if debug:
        print('Attempting to correct barcodes : {}, {}, {}'.format(bc1, bc2, bc3))

    bc1_matches, edit_dist1 = get_min_edit_dists(bc1, bc1_dict, max_d=bc_edit_dist)
    bc2_matches, edit_dist2 = get_min_edit_dists(bc2, bc2_dict, max_d=bc_edit_dist)
    bc3_matches, edit_dist3 = get_min_edit_dists(bc3, bc3_dict, max_d=bc_edit_dist)

    entry = {}

    if len(bc1_matches) > 1 or len(bc1_matches) == 0:
        entry['bc1_edit_dist'] = 4
        entry['bc1_corrected'] = ''
    if len(bc2_matches) > 1 or len(bc2_matches) == 0:
        entry['bc2_edit_dist'] = 4
        entry['bc2_corrected'] = ''
    if len(bc3_matches) > 1 or len(bc3_matches) == 0:
        entry['bc3_edit_dist'] = 4
        entry['bc3_corrected'] = ''

    if len(bc1_matches) == 1:
        entry['bc1_edit_dist'] = edit_dist1
        entry['bc1_corrected'] = bc1_matches[0]
    if len(bc2_matches) == 1:
        entry['bc2_edit_dist'] = edit_dist2
        entry['bc2_corrected'] = bc2_matches[0]
    if len(bc3_matches) == 1:
        entry['bc3_edit_dist'] = edit_dist3
        entry['bc3_corrected'] = bc3_matches[0]

    #for all barcodes
    expected_BC1 = "AACGTGAT"
    #barcode 3 of the parse set, because i was an idiot.
    expected_BC2 = "ACTCGTAA"

    if entry['bc1_corrected'] == expected_BC1 or entry['bc1_corrected'] == expected_BC2:
        entry['is_true_BC1'] = True
    else:
        entry['is_true_BC1'] = False
    if entry['bc2_corrected'] == expected_BC1 or entry['bc2_corrected'] == expected_BC2:
        entry['is_true_BC2'] = True
    else:
        entry['is_true_BC2'] = False
    if entry['bc3_corrected'] == expected_BC1 or entry['bc3_corrected'] == expected_BC2:
        entry['is_true_BC3'] = True
    else:
        entry['is_true_BC3'] = False

    return entry

def is_duplicate_x(x, bc_dict):
    bc_set = bc_dict.keys()
    bc = x.bc
    umi = x.umi
    entry = None
    if bc in bc_set:
        umi_set = bc_dict[bc]
        if umi in umi_set:
            entry = True
        else:
            entry = False
            #bc_dict[bc].append(umi)
    else:
        entry = False
        #bc_dict.update({bc:[umi]})

    return entry

def linker_type_x(x, l1_min, l2_min):
    entry = None
    if (x.l1_score >= l1_min) & (x.l2_score >= l2_min):
        entry = 1
    elif (x.l1_score < l1_min) & (x.l2_score >= l2_min):
        entry = 2
    elif (x.l1_score >= l1_min) & (x.l2_score < l2_min):
        entry = 3
    elif (x.l1_score < l1_min) & (x.l2_score < l2_min):
        entry = 4

    return entry
###################################################################################
############################### Processing steps ##################################
###################################################################################

def fastq_to_df(fname, oprefix, chunksize, verbose=1):
    """
	Save the input fastq file into a table format into a file called
	oprefix_table.tsv.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress

	Returns:
		ofile (str): Name of output file
	"""
    # get the sequences from each read
    seqs = []
    read_names = []
    strands = []
    i = 1
    ofile = '{}_table.tsv'.format(oprefix)
    with open(fname, 'r') as f:
        while True:
            read_name = f.readline().strip()
            read_name = read_name[1:]
            if len(read_name)==0:
                break
            seq = f.readline().strip()
            strand = f.readline().strip()
            qual = f.readline()
            seqs.append(seq)
            strands.append(strand)
            read_names.append(read_name)

            # print out notification and write to file
            if i % chunksize == 0 and i != 1 :
                if verbose == 2:
                    print('Processed {} reads'.format(i))

                # pack everything into a dataframe
                df = pd.DataFrame(seqs)
                df.columns = ['seq']
                df['read_name'] = read_names
                df['strand'] = strands

                # first write
                if i == chunksize:
                    df.to_csv(ofile, sep='\t', index=False)
                    read_names = []
                    strands = []
                    seqs = []
                else:
                    df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')
                    read_names = []
                    strands = []
                    seqs = []
            i += 1

    # cleanup
    if len(seqs) > 0:
        df = pd.DataFrame(seqs)
        df.columns = ['seq']
        df['read_name'] = read_names
        df['strand'] = strands
        df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')
    return ofile

def score_linkers(fname, oprefix, t=1, verbose=1, chunksize=10**5, delete_input=False , linker='split'):
    """
	Find and report highest linker scores in each read.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
    """

    l1, l2 = get_linkers(linker, verbose)
    l_seqs = [l1, l2]
    l_prefs = ['l1', 'l2']

    # loop through chunks of the file
    i = 0
    ofile = '{}_seq_linker_alignment_scores.tsv'.format(oprefix)
    print(chunksize)
    for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
        if t == 1:
            l_df = df.apply(score_linkers_x, args=(l_seqs, l_prefs), axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            l_df = df.parallel_apply(score_linkers_x, args=(l_seqs, l_prefs), axis=1, result_type='expand')

        df = pd.concat([df, l_df], axis=1)

        # first write
        if i == 0:
            df.to_csv(ofile, sep='\t', index=False)
        else:
            df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

        i += chunksize
        if verbose > 1:
            print('Found linker scores for {} reads'.format(i))

    # delete input file
    if delete_input:
        os.remove(fname)

    return ofile

def align_linkers(fname, oprefix, t=1, chunksize=10**5, mm=0.1, verbose=1, delete_input=False, linker='split'):
    """
	Find indices of highest-scoring linker in each read.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on

		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
    """

    l1, l2 = get_linkers(linker, verbose)

    # calculate lowest scores viable for l1 and l2
    l1_min = len(l1) - math.ceil(len(l1)*mm)
    l2_min = len(l2) - math.ceil(len(l2)*mm)

    # loop through chunks of the file
    i = 0
    # create a umi/cell dict
    cell_umi_dict = {}
    ofile = '{}_seq_linker_alignments.tsv'.format(oprefix)
    for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
        # find the start and end of each alignment for the valid reads
        if t == 1:
            l_df = df.apply(align_linkers_x, args=(l1, l2), axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            l_df = df.parallel_apply(align_linkers_x, args=(l1, l2), axis=1, result_type='expand')
        df = pd.concat([df, l_df], axis=1)

        #remove all umi's that are not 10 in length, and then drop duplicates keeping first
        df = df[df['umi_len'] == 10]
        df.drop_duplicates(subset=['umi','bc'], keep='first', inplace=True)
        #check if the umi's already exist in the previous dataframes, if so remove them.

        if t == 1:
            df['duplicate'] = df.apply(is_duplicate_x, args=[cell_umi_dict], axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            df['duplicate'] = df.parallel_apply(is_duplicate_x, args=[cell_umi_dict], axis=1, result_type='expand')

        #df['duplicate'] = df.apply(is_duplicate_x, args=[cell_umi_dict], axis=1, result_type='expand')
        df = df[df['duplicate'] == False]
        df.reset_index(drop=True, inplace=True)
        #in the chance that dataframe becomes empty, continue to the next cycle
        if len(df.index) == 0:
            continue

        #adding to cell_umi_dict
        for j in df.index:
            if df.loc[j,'bc'] in cell_umi_dict.keys():
                cell_umi_dict[df.loc[j,'bc']].append(df.loc[j,'umi'])
            elif df.loc[j,'bc'] not in cell_umi_dict.keys():
                cell_umi_dict[df.loc[j,'bc']] = [df.loc[j,'umi']]

        #drop useless columns
        df.drop(['duplicate', 'umi_len'], axis=1, inplace=True)

        df.reset_index(drop=True, inplace=True)

        if t == 1:
            df['linker_relation'] = df.apply(linker_type_x, args=(l1_min, l2_min), axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            df['linker_relation'] = df.parallel_apply(linker_type_x, args=(l1_min, l2_min), axis=1, result_type='expand')

        #link_relation = []
        #for j in df.index:
        #    if (df.loc[j, "l1_score"] >= l1_min) & (df.loc[j, "l2_score"] >= l2_min):
        #        link_relation.append(1)
        #    elif (df.loc[j, "l1_score"] < l1_min) & (df.loc[j, "l2_score"] >= l2_min):
        #        link_relation.append(2)
        #    elif (df.loc[j, "l1_score"] >= l1_min) & (df.loc[j, "l2_score"] < l2_min):
        #        link_relation.append(3)
        #    elif (df.loc[j, "l1_score"] < l1_min) & (df.loc[j, "l2_score"] < l2_min):
        #        link_relation.append(4)
        #df['linker_relation'] = link_relation

        df.reset_index(drop=True, inplace=True)

        # first write
        if i == 0:
            df.to_csv(ofile, sep='\t', index=False)
        else:
            df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

        i += chunksize
        if verbose > 1:
            print('Found linkers for {} reads'.format(i))

    # delete input file
    if delete_input:
        os.remove(fname)

    return ofile

def align_halves(fname, oprefix, t=1, verbose=1, chunksize=10**5, linker='split'):
    l1, l2 = get_linkers(linker, verbose)
    l1_half = int(len(l1) / 2)
    l2_half = int(len(l2) / 2)
    l1_half1, l1_half2 = l1[:l1_half], l1[l1_half:]
    l2_half1, l2_half2 = l2[:l2_half], l2[l2_half:]

    # loop through chunks of the file
    i = 0
    ofile = '{}_seq_linker_halves_alignments.tsv'.format(oprefix)
    for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
        # find the start and end of each alignment for the valid reads
        if t == 1:
            l_df = df.apply(align_linker_halves_x, args=(l1_half1, l1_half2, l2_half1, l2_half2), axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            l_df = df.parallel_apply(align_linker_halves_x, args=(l1_half1, l1_half2, l2_half1, l2_half2), axis=1, result_type='expand')
        df = pd.concat([df, l_df], axis=1)

        df.reset_index(drop=True, inplace=True)

        if i == 0:
            df.to_csv(ofile, sep='\t', index=False)
        else:
            df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

        i += chunksize
        if verbose > 1:
            print('Found linkers for {} reads'.format(i))

    return ofile

def get_lev_dist(fname, oprefix, bc_edit_dist, chunksize, t, verbose):
    bc1_edit_dict, bc2_edit_dict, bc3_edit_dict = load_barcodes()

    i=0
    for df in pd.read_csv(fname, chunksize=chunksize, sep="\t"):
        if t == 1:
            tmp_df = df.apply(correct_bcs_x, args=(bc_edit_dist, bc1_edit_dict, bc2_edit_dict, bc3_edit_dict), axis=1, result_type='expand')
        else:
            pandarallel.initialize(nb_workers=t, verbose=1, use_memory_fs=False)
            tmp_df = df.parallel_apply(correct_bcs_x, args=(bc_edit_dist, bc1_edit_dict, bc2_edit_dict, bc3_edit_dict), axis=1, result_type='expand')
        df = pd.concat([df, tmp_df], axis=1)

        df.reset_index(drop=True, inplace=True)

        ofile = oprefix + "_corrected_bcs.tsv"
        if i == 0:
            df.to_csv(ofile, sep='\t', index=False)
        else:
            df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

        i += chunksize
        if verbose > 1:
            print('corrected {} barcodes'.format(i))

def plotting_graphs(oprefix, linker, plot, verbose, mm):
    l1, l2 = get_linkers(linker, verbose)
    l1_len = len(l1)
    l2_len = len(l2)
    if plot == 'score' or plot == 'all':
        if path.exists('{}_seq_linker_alignments.tsv'.format(oprefix)):
            df = pd.read_csv('{}_seq_linker_alignments.tsv'.format(oprefix), sep='\t')
        else:
            df = pd.read_csv('{}_seq_linker_alignment_scores.tsv'.format(oprefix), sep='\t')
        df.reset_index(inplace=True)
        plot_linker_scores(df, oprefix, l1_len, l2_len, mm)

    if plot == 'heatmap' or plot == 'all':
        if path.exists('{}_seq_linker_alignments.tsv'.format(oprefix)):
            df = pd.read_csv('{}_seq_linker_alignments.tsv'.format(oprefix), sep='\t')
        else:
            df = pd.read_csv('{}_seq_linker_alignment_scores.tsv'.format(oprefix), sep='\t')
        df.reset_index(inplace=True)
        plot_linker_heatmap(df, oprefix, l1_len, l2_len, mm)

    if plot == 'link_rel' or plot == 'all':
        if path.exists('{}_seq_linker_alignments.tsv'.format(oprefix)):
            df = pd.read_csv('{}_seq_linker_alignments.tsv'.format(oprefix), sep='\t', usecols=[3, 4, 17])
        else:
            df = pd.read_csv('{}_seq_linker_alignment_scores.tsv'.format(oprefix), sep='\t', usecols=[3, 4, 5])
        df.reset_index(inplace=True)
        plot_linker_relation(df, oprefix, verbose)

    if plot == 'link_anal' or plot == 'all':
        df = pd.read_csv('{}_seq_linker_alignments.tsv'.format(oprefix), sep='\t', usecols=[0, 3, 4, 5, 6, 7, 8, 9, 10, 17])
        df.reset_index(inplace=True)
        linker_analysis(df, oprefix, l1, l2, l1_len, l2_len, verbose)
