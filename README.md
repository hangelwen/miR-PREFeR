# miR-PREFeR: <strong>mi</strong>cro<strong>R</strong>NA <strong>PRE</strong>diction <strong>F</strong>rom small <strong>R</strong>NAseq data #


**The miR-PREFeR pipeline is still in active development, to be able to use the newest features from the pipeline, it's best to check the  <https://github.com/hangelwen/miR-PREFeR> page and obtain the newest version. The current version is only tested under Python 2.6.7, Python 2.7.2 and Python 2.7.3 and should work under Python 2.6.* and Python 2.7.\* (Tested platforms: Linux, Mac OS). It does NOT work under Python 3.0 currently. If you find any problem, please contact the author.**

## 1. Required programs ##
To run the miR-PREFeR pipeline, the ViennaRNA package(version 1.8.5 or 2.1.x) and samtools(0.1.15 or later. Tested under 0.1.18) should be installed on the system. miR-PREFeR uses samtools commands to manipulate SAM and BAM alignment files, and uses RNALfold from the ViennaRNA package to do RNA secondary structure folding. The miR-PREFeR pipeline takes SAM alignment files as input, so an aligner such as Bowtie is also required (Not required by the pipeline, but needed for preparing the input data).

### ViennaRNA package ###
The ViennaRNA package can be downloaded from <http://www.tbi.univie.ac.at/~ronny/RNA/index.html.>. The website provides precompiled  packages for some platforms. For these platforms, download the corresponding package and install it according to the platform package management system. If no precompiled package is provided for your platform, you can install it from source code.
To install the package from source code, simply decompress the source code archive and cd into the folder,
then:

    ./configure
    make

and (as root):
    make install

If the user does not have root permission. then the package can be installed to
user specified locations:

    ./configure --prefix="/path/to/install" --without-perl
    make
    make install

After finish installing the package, **add the bin directory to the PATH environment variable.** Instructions on how to add a path to the PATH environment variable can be found at <http://www.cyberciti.biz/faq/unix-linux-adding-path/> or just google "`add directory to PATH`".


**NOTE: Because that RNALfold from the ViennaRNA package version 2.0.4 has a bug (If the input sequence has no valid folding structure, the program produces a segmentation fault), please make sure to use ViennaRNA package version 1.8.x or the newest 2.1.2.** The 1.8.* version is easier to install than the 2.1.2 version. Because the new version needs the `gengetopt` package (<http://www.gnu.org/s/gengetopt/gengetopt.html>), which may not be available on your system.

### Samtools ###
Samtools can be downloaded from <http://samtools.sourceforge.net/>. Please follow the instructions from the package to install it. Please note the version 0.1.15 or later is needed (The pipeline uses the `samtools depth` command, which was introduced from version 0.1.15).


## 2. Obtain and install the pipeline ##

The miR-PREFeR pipeline is hosted on Github (https://github.com/hangelwen/miR-PREFeR). To obtain the program, open an terminal window and execute the following command:

    git clone https://github.com/hangelwen/miR-PREFeR.git

This makes an directory named 'miR-PREFeR' in the current directory, and clones the newest version of the miR-PREFeR pipeline code into the directory. The pipeline is ready to use.

If you do not have `git` on your system, simply go to <https://github.com/hangelwen/miR-PREFeR/releases>, download the zip file or the tar.gz file of the latest version.

**The miR-PREFeR pipeline is still in active development, to be enable to use the newest features from the pipeline, it's better to check the  <https://github.com/hangelwen/miR-PREFeR> page and obtain the newest version.**

## 4. Test the pipeline. ###

The package contains a small dataset for testing whether the pipeline works after you download it. Please refer to [HOW\_TO\_RUN_EXAMPLE.txt](HOW_TO_RUN_EXAMPLE.txt) file to see more infomation about testing.

## 4. How to run the pipeline. ##

### a. Prepare input data for the pipeline. ###

The miR-PREFeR pipeline needs the following input:

1. A fasta file, which contains the gnome sequences of the species under study.
2. one or more SAM files which contains the alignments of small RNAseq data with the gnome.
3. (Optional) An GFF (<http://www.sanger.ac.uk/resources/software/gff/spec.html>) file which lists regions in the gnome sequences that should be ignored from miRNA analysis.

#### a). Genome fasta file ####
Fasta format specification can be found at <http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml>. In miR-PREFeR, when an identifier of a genomic sequence is needed, the first part that does not contain any white space characters (whitespace, tab, etc) is used. For example, for the following sequence, 'ath-MIR773a' is used as the identifier of the seqeunce. **Thus, please ensure that all the sequences in the FASTA files have different identifiers.**

    >ath-MIR773a MI0005103
    AGGAGGCAAUAGCUUGAGCAAAUAAUUGAUUGCAGAAGUCCAUCGACUAAAGCUGUCACCUGUUUGCUUCCAGCUUUUGUCUCCU

#### b). SAM alignment files ####
The miR-PREFeR pipeline takes SAM format alignment files. SAM alignment files can be generated by many aligners. Here we use Bowtie (<http://bowtie-bio.sourceforge.net/index.shtml>) as an example.

**Prepare RNAseq fasta files for Bowtie**

Because miR-PREFeR is able to use multiple small RNAseq data as input and utilizes information from multiple RNAseq samples (expression pattern consistency, differential expression, etc) to increase the performance of the miRNA loci prediction. **Thus, It's an requirement that the RNAseq reads data (fasta files that contain the RNAseq reads) be preprocessed by the provided `process-reads-fasta.py` script.** The script takes an list of **uncollapsed fasta files** as input, and renames the reads in the input files, and produces a list of collapsed fasta files. Here `uncollpased fasta` means that identical reads in the fasta file have multiple entries.

For example, if you have 3 small RNAseq samples in uncollpased fasta format with names `SAMPLE1.fasta`, `SAMPLE2.fasta`, and `SAMPLE3.fasta`, then by executing:

    python process-reads-fasta.py SAMPLE1.fasta SAMPLE2.fasta SAMPLE3.fasta

Three new files with name `SAMPLE1.fasta.processed`, `SAMPLE2.fasta.processed`, and `SAMPLE3.fasta.processed` are produced. These three files are then aligned to the genome fasta file using Bowtie.

**Align the RNA-seq fasta files with Bowtie:**

a1. indexing the genome sequences:

    bowtie-build -q -f TAIR10.fas bowtie-index/TAIR10

a2. alignment

    bowtie -a -v 0 -p 8 -m 10 -S bowtie-index/TAIR10  -f SAMPLE1.fasta.processed > SAMPLE1.sam 2> SMAPLE1.log
    bowtie -a -v 0 -p 8 -m 10 -S bowtie-index/TAIR10  -f SAMPLE2.fasta.processed > SAMPLE2.sam 2> SMAPLE2.log
    bowtie -a -v 0 -p 8 -m 10 -S bowtie-index/TAIR10  -f SAMPLE3.fasta.processed > SAMPLE3.sam 2> SMAPLE3.log

The -v 0 option allows no mismatch in the alignment
The -m 10 option discards reads that can be mapped to more than 10 positions.
The -p 8 option uses 8 threads to do the alignment


#### c). The (optional) gff file  ####

For a genome that already has known annotations, some regions in the genome can be ignored when doing the miRNA analysis. For example, many species have protein coding sequence (CDS) annotations, so there is no need to do run miR-PREFeR on those regions. These regions should be listed in a GFF file (<http://www.sanger.ac.uk/resources/software/gff/spec.html>).

### b. Prepare a configuration file for the pipeline. ###

The miR\_PREFeR.py script takes a configuration file as input. The configuration file lists all the information (such as where are the SAM files, genome fasta file, etc) needed to run the pipeline. An [example configuration file](example/config.example) with detailed descriptions of each option is in the `example` folder. The options are also explained in the following section.

1. PIPELINE_PATH: The path of the miR-PREFeR pipeline.
2. FASTA_FILE: The path of the genome sequences in fasta format. Absolute path preferred.
3. ALIGNMENT_FILE: The path of the RNAseq alignment files in SAM format. Absolute path preferred.
4. GFF_FILE: The path of the optional GFF file. Absolute path preferred.
5. PRECURSOR_LEN: The max length of a miRNA precursor. The default is 300
6. READS_DEPTH\_CUTOFF: The first step of the pipeline is to identify candidate regions of the miRNA loci. If READS\_DEPTH\_CUTOFF = N, then genomic position that the mapped depth is smaller than N is not considered. Default value is 20.
7. NUM\_OF\_CORE: Number of CPUS/Cores avalible for this computation. If commented out or leave blank, 1 is used.
8. OUTFOLDER: Outputfolder. If not specified, use the current working directory.
9. NAME\_PREFIX: Prefix for naming the output files. For portability, please DO NOT contain any spaces and special characters. The prefered includes 'A-Z', 'a-z', '0-9', and underscore '_'.
10. MAX_GAP: Maximum gap length between two contigs to form a candidate region. Default value is 100.

### c. Run the pipeline. ###

With the configuration file ready, the pipeline can be used to predict miRNA loci. The miR-PREFeR pipeline can be run as following:

    python miR_PREFeR.py [options] command configfile

Currently, the following options are avalible:

1. **-h**: Show help information.
2. **-l**: Generate a log file. From the log file you can find the status during the running of the pipeline. **Recommend to always use this option.**
3. **-k**: After finish the whole pipeline, do not remove the temporary folder that contains the intermediate files. This will save disk space. If it's not specified, you can delete the temporary folder after getting the result.

`command` could be one of `check`, `pipeline`, `prepare`, `candidate`, `fold`, `predict`, and `recover`. (Run `python miR_PREFeR.py -h` to see help on options and commands.)

1. **check:** Check the presence of the depended programs (RNALfold and samtools), and check the recovery information (Shows which stage the previous computation on the same data was ceased. See the `recover` command in 6.).
2. **pipeline:** Run the whole pipeline. That is, run `prepare`, `candidate`,`fold` and `predict` sequentially. This is the normal way to run miR-PREFeR.
3. **prepare:**  Run the first step of the pipeline. This step prepares some data files needed in the following steps.
4. **candidate:**  Identify possible candidate regions. This step can ONLY be run if the `prepare` step has been finished on the configfile file.
5. **fold:** Fold the candidate regions. This step can ONLY be run if the `prepare` and `candidate` steps have been finished on the configfile file.
6. **predict:** Predict miRNA loci. This step can ONLY be run if the `prepare`, `candidate` and `fold` steps have been finished on the configfile file.
7. **recover:** Tries its best to recover an unfinished job and continues to run from the unfinished step it was ceased. This is designed to easily continue a job other than re-run the whole pipeline from start. For example, if a job was kill half way (This happens sometime. For example, the job takes too long time and the user killed it. Or, if the job was run on an Cluster and was killed halfway because some resources exceed the max values.) In these cases, one does not need to run the pipeline from start, the job can be continued by using the `recover` command:

    `python miR_PREFeR.py -l recover configfile`

**NOTE on job recovery:** To recover and continue a job, the intermediate output files in the temporary folder (See the Output section) should not be removed. All files needed to do recovery are in this folder.


## 4. Output ##

The predicted miRNA loci are written to a GFF file named NAME\_PREFIX + "_miRNA.gff3" in the output folder, where 'NAME\_PREFIX' is the prefix for naming output files which is specified in the configuration file. The mature sequences and the stem loop sequences of the predicted miRNAs are written to two fasta files, named NAME\_PREFIX + "\_miRNA.mature.fa" and NAME\_PREFIX + "\_miRNA.stemloop.fa".

During the running of the pipeline, there is a folder named NAME\_PREFIX + "_tmp" in the output folder. The folder contains intermediate files generated by the pipeline. If the job is finished without interruption, the folder is deleted by default (Specify the `-k` option can keep the folder, but this is usually not needed.).  **If the job is interrupted half way, then the folder is not removed by default, so that one can run the job in the `recover` mode to continue the computation. (Because all files needed to do recovery are in this temp folder, so if one wants to recover an interrupted job, please do not remove the folder. The folder will be removed automatically once the whole pipeline is successfully run.).**
