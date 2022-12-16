# LinkerComparison
Thesis Chapter 1 - Comparing post-sequencing processing pipelines, and linker sequences of the SPLiT-seq technique

## PipelineComparison
Here six computational pipelines are used on a single dataset from [Rosenberg et. al.](https://www.science.org/doi/10.1126/science.aam8999?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) to compare their performance.
For each pipeline a bash.sh script was created that creates a timestamp on start and stop of the bash script. Maximum RAM usage was set to 64 GB with the linux ulimit command.
Described in more detail [here](link to paper if it is ever published)

## LinkerComparison
This folder contains all the scripts used to analyze the linker designs of a [SPLiT-seq](https://www.science.org/doi/10.1126/science.aam8999?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) experiment.

Code is heavily based/inspired from [LR-splitpipe](https://github.com/fairliereese/LR-splitpipe)

### Running analysis
To perform the linker analysis use demultiplex.py
```
Usage: python ./Linkercomparison/demultiplex.py [options]

Options:
  -m --mode     Which mode the script should be run in available modes are fastq_to_df, score_linkers, align_linkers, align_linker_halves, get_lev_dist, plot
  -f --infile   The reads.fastq or table.tsv that should be used as input file.
  -o --output   The location and pretext of the output. e.g. /location/of/output will create an output file in the folder '/location/of/' called 'output_ofthemodeused.file'
  --chunksize   The amount of rows a read in at any given point in time during some of the steps of the script. default is 100000. This is to prevent memory overloading.
  --verbosity   The amount of text is written to the screen as the scripts runs.
  -t --threads  The amount of threads or cores used.
  --delete_input  whether to delete the input used in the script.
  -l --linker   which linkertype is used in the script, default is split. options availaible: split, GtoT, 2as1, share, parse
  --mm          percent mismatches allowed in each linker depending on the linker length. default is 0.1
  -p --plot_type  Which plot should be created. options are scores, heatmap, link_rel, link_align, all. 
  -e --edit_dist  edit distance used to correct the barcodes. default is 2.
```
