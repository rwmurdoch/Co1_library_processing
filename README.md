# Co1_library_processing
A set of Qiime2 and R scripts that process Co1 datasets, assuming that a Co1 classifier is on-hand

This set of scripts has been applied to a Co1 amplicon library dataset, consisting of over 400 samples.  I applied my experimental BOLD classifier, build from all species Co1 records from the BOLD database using an as-of-yet unannounced process.  If you would like to experiment with this classifier, please contact me at rmurdoch@utk.edu.  

Some sort of classifier is necessary for the script to function.  If you just want to run the script and produce ASVs and feature table, you can use a dummy classifer, such as the Silva classfier distributed by the Qiime2 group.  Just ignore the taxonomic placements.

Qiime2 must be installed via the conda distribution.  Clear instructions are available at the Qiime2 website.

From a practical perspective,:
1. make a folder for the project of interest; all scripts will be run from within this folder. Make sure that the accessory scripts "size.filter.seqs.q2.R", "create.OTU.table.R", and "extract.readfiles.sh" are all in the master project folder. 
2. create a folder entitled "reads" and copy in all relevant read-files. create an empty folder entitled "exports".
3. create a sample-by-sample metadata.csv file according to qiime2 metadata file format (https://docs.qiime2.org/2019.1/tutorials/metadata/)
4. if your reads are within subfolders in the "reads" folder, run the script "extract.readfiles.sh", included in this repository.
5. within the qiime2Co1.sh script, find and replace all instances of "willcox" with the project name of your choosing. The other scripts will interact with the master script without any modification.
6. manually scan the script for other required changes, such as classifier, dada2 trim lengths, frequency cutoffs for the feature-table, etc.
7. if desired, change the two sequence length limits in the "size.filter.seqs.q2.R" script.  Currently, it set to allow only ASVs withing 147 to 167 bp.
8. a file called "results.tar.gz" contains many of the key output files.  The "exports" folder

The general outline of the processing pipeline, from a data processing perspective:

1. Amplicon library data was processed using the Qiime2-2018.11 environment.
2. Primers were trimmed from reads using cutadapt.
3. Actual sequence variants were determined using dada2 with forward and reverse truncation length of 140bp
4. ASVs shorter than 147bp or longer than 167bp were removed.
5. ASVs occurring fewer than twice across all samples were removed.
6. ASVs were taxnomically placed using a custom sklearn feature classifier was trained using all BOLD Co1 database entries, downloaded in June, 2018.


