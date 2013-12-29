Script descriptions.

1. process-reads-fasta.py
Pre-process a list of uncollapsed fasta files. Collapse the identical reads and
rename the reads. The new identifier for each read in the output files has the
following format: SampleName_rA_xB. Here A is an number identifies the read,
which an sequentially increasing number assigned to the read by the script. B is
the total count of the read (depth of the read). 'SampleName' is the name of the
sample/tissue/library that the input fasta file represents.

For example, if the sample name is 'root', and the read occurred 120 times in
the library, then it's identifier could be root_r23_x120. Here '23' is just a
number that donates the order of the reads when processing it, it has no use
otherwise.

Input:
1. samplenamlist: a file in which the sample names are listed. Each
sample name should be on one line. Sample name should be one word without any
blank characters. The number of sample names in the file should be the same as
the number of input fasta files. Each name should be unique.

2. list of fasta files: a list of uncollapsed fasta files, each represents a
RNAseq dataset from a sample/library/tissue/etc. The order of the files should
be consistent with the names listed in the samplenamlist file.


The script assumes each read sequence occupies ONLY one line in the input file.

The output files are in the same folder as the input files, with ".processed "
as suffix.


2. convert-mirdeep2-fasta.py
Pre-precess a list of collapsed fasta format read files. The input read files
following the miRDeep2 naming convention: the identifier of each read ends with
"_xN", where 'N' is the depth of the read in the library. The new identifier for
each read in the output files has the following format: SampleName_rA_xB. Here A
is an number identifies the read, which an sequentially increasing number
assigned to the read by the script. B is the total count of the read (depth of
the read). 'SampleName' is the name of the sample/tissue/library that the input
fasta file represents.

The script assumes each read sequence occupies ONLY one line in the input file.

The output files are in the same folder as the input files, with ".processed "
as suffix.


3. convert-readcount-file.py
Pre-precess a list of read count files. In a read count file, each line contains
two columns: the first column is a read sequence, the second column is the depth
of the read. The two columns can be seperated by spaces or tabs.  The new
identifier for each read in the output files has the following format:
SampleName_rA_xB. Here A is an number identifies the read, which an sequentially
increasing number assigned to the read by the script. B is the total count of
the read (depth of the read). 'SampleName' is the name of the
sample/tissue/library that the input fasta file represents.

The output files are in the same folder as the input files, with ".processed "
as suffix.


4. bowtie-align-reads.py
Align fasta format reads to genome using bowtie.

The input fasta format should follow the output of
'process-reads-fasta.py'. That is, the read should be collapsed; the read ID
must be in format of 'samplename_rA_xN'. Here samplename is the name of the
sample from with small RNASeq library came from, 'A' is an unique number to
identify the read, 'N' is the depth of the read.

The output of the script are SAM format alignment files, which can be then used
in mir_PREFeR.py to do miRNA prediction. The output are in the same folder as
the input files, with suffix ".sam".

NOTE: Pair-end reads are not supported currently. Not all bowtie options are
exposed in the script. To use more bowtie options, please use bowtie
directly to align the reads.

To use the script:
python bowtie-align-reads.py [options] <read fasta1> ... <read fastaN>

Example:
python bowtie-align-reads.py -p 2 -m 20 -f -r TAIR10.fa SAMPLE1.fa.processed

This command maps reads in SAMPLE1.fa.processed to TAIR10.fa using two
threads (-p). Reads that mapped to more than 20 positions are not
reported (-m). Unmapped alignments are filtered using SAMtools (-f).

To see the available options, using:
python bowtie-align-reads.py -h
