# LinkerComparison
Thesis Chapter 1 - Comparing post-sequencing processing pipelines, and linker sequences of the SPLiT-seq technique

## PipelineComparison
Here six computational pipelines are used on a single dataset from [Rosenberg et. al.](https://www.science.org/doi/10.1126/science.aam8999?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) to compare their performance.
For each pipeline a bash.sh script was created that creates a timestamp on start and stop of the bash script. Maximum RAM usage was set to 64 GB with the linux ulimit command.
the venndiagram.py is used to create the UsePlot which compares the Cell Barcode output of the pipelines compared.
Described in more detail [here](link to paper if it is ever published)

## LinkerComparison
This folder contains all the scripts used to analyze the linker designs of a [SPLiT-seq](https://www.science.org/doi/10.1126/science.aam8999?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) experiment.

Code is heavily based/inspired from [LR-splitpipe](https://github.com/fairliereese/LR-splitpipe)

### Running analysis
To perform the linker analysis use demultiplex.py
```
Usage: python ./LinkerComparison/demultiplex.py [options]

Options:
  -m --mode     Which mode the script should be run in available modes are fastq_to_df, score_linkers, align_linkers, align_halves, get_lev_dist, plot
  -f --infile   The reads.fastq or table.tsv that should be used as input file.
  -o --output   The location and pretext of the output. e.g. /location/of/output will create an output file in the folder '/location/of/' called 'output_ofthemodeused.file'
  --chunksize   The amount of rows a read in at any given point in time during some of the steps of the script. default is 100000. This is to prevent memory overloading.
  --verbosity   The amount of text is written to the screen as the scripts runs.
  -t --threads  The amount of threads or cores used.
  --delete_input  whether to delete the input used in the script.
  -l --linker   which linkertype is used in the script, default is split. options availaible: split, GtoT, 2as1, share, parse
  --mm          percent mismatches allowed in each linker depending on the linker length. default is 0.1
  -p --plot_type  Which plot should be created. options are score, heatmap, link_rel, link_align, all. 
  -e --edit_dist  edit distance used to correct the barcodes. default is 2.
```

### Basic pipline steps
#### fastq to df
This function is copied from [LR-splitpipe](https://github.com/fairliereese/LR-splitpipe) and turns the .fastq file into a table.tsv where each row represent a read.
```
Example:
python ./LinkerComparison/demultiplex.py -m fatq_to_df -f /path/to/test.fastq -o /path/to/test --chunksize 100 --verbosity 2
```
#### score_linkers
Function is largely copied as well from [LR-splitpipe](https://github.com/fairliereese/LR-splitpipe).
It creates file with scores of your linkers, does not try to obtain any extra data. If you want a 'quick' look this works best.
To visualize you can use the plot options to display score or heatmap
```
Example:
python ./LinkerComparison/demultiplex.py -m score_linkers -f /path/to/test_table.tsv -o /path/to/test -t n --chunksize 100 --verbosity 2 --linker split
python ./LinkerComparison/demultiplexe.py -m plot -p score -o /path/to/test --verbosity 2 --linker split
```
#### align_linkers
This functions works the same as score_linkers, however instead of only extracting the Alignment Score it also get other parameters.
which are [Start,Stop,length, linker_relation, umi, umi_len, bc1, bc2, bc3, CB]
This function will also perform deduplicating of UMIs increasing the runtime significantly, especially with large datasets.
It can be followed by plotting score, heatmap, link_rel or link_anal
```
Example:
python ./LinkerComparison/demultiplex.py -m align_linkers -f /path/to/test_table.tsv -o /path/to/test -t n --chunksize 100 --verbosity 2 --linker split
python ./LinkerComparison/demultiplexe.py -m plot -p link_rel -o /path/to/test --verbosity 2 --linker split
```

#### align_linker_halves
Works the same as align_linker, except instead of aligning the full linkers against the sequence will slice the linker in half at the middle and perform alignments with each seperate half. Extracts the same paramters as align_linkers but then from their respective halves.
```
Example:
python ./LinkerComparison/demultiplex.py -m align_halves -f /path/to/test_table.tsv -o /path/to/test -t n --chunksize 100 --verbosity 2 --linker split
python ./LinkerComparison/demultiplexe.py -m plot -p link_rel -o /path/to/test --verbosity 2 --linker split
```

#### get_lev_dist
This function will perform barcode correction on the extracted BC from align_linkers, returning the Levensthein distance and the sequence of the corrected barcode.
```
Example:
python ./LinkerComparison/demultiplex.py -m get_lev_dist -f /path/to/test_table.tsv -o /path/to/test -t n --chunksize 100 --verbosity 2 -e 2
```
#### plot
as described in previous modes, -m plot together with -p 'plot_type' can be used to make graphical images of the tables generated.
the script should provide an error when input is incompatible. e.g. using an score_linkers output to create link_rel or link_anal plots.
This is a short description of the plots generated.
score:    a histogram of the Alignment scores of the linkers per read.
heatmap:  a heatmap of the Alignment scores of the linkers
link_rel: a histogram of the Alignment scores of the linkers per read. In color are the four different linker states that an occur.
link_anal: a lineplot of the percentage of correct bases per persition compared to the reference. additionally violin plots of the linker length, start and stop positions.

### figure_x, y and supp_figure_y
These scripts are used to create the figures in this [paper](link to paper)
Four different linker designs were tested to create these plots and are currently required to run the scripts. 
