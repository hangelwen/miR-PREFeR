import sys
import subprocess
import multiprocessing
import os
import re
import tempfile
import distutils
import distutils.spawn

"""Align fasta format reads to genome using bowtie.

The input fasta format should follow the output of
'process-reads-fasta.py'. That is, the read should be collapsed; the read ID
must be in format of 'samplename_rA_xN'. Here samplename is the name of the
sample from with small RNASeq library came from, 'A' is an unique number to
identify the read, 'N' is the depth of the read.

The output of the script are SAM format alignment files, which can be then used
in mir_PREFeR.py to do miRNA prediction. The output are in the same folder as
the input files.

To use the script:
python bowtie-align-reads.py [options] <fasta1> <fasta2> ... <fastaN>

To see the available options, using:
python bowtie-align-reads.py -h

"""


def parse_option_optparse():
    from optparse import OptionParser
    helpstr = """python bowtie-align-reads.py [options] <fasta1> ... <fastaN>

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
    python bowtie-align-reads.py -p 2 -k 20 -f -r TAIR10.fa SAMPLE1.fa.processed

    This command maps reads in SAMPLE1.fa.processed to TAIR10.fa using two
    threads (-p). Reads that mapped to more than 20 positions are not
    reported (-k). Unmapped alignments are filtered using SAMtools (-f).

    To see the available options, using:
    python bowtie-align-reads.py -h

    """

    parser = OptionParser(helpstr)

    parser.add_option("-r", "--reference",
                      help="Reference genome in fasta format. If you have multiple reference files, please use multipe -r options. If you have aready index the reference sequences, you should use the -i option. ", action="append")

    parser.add_option("-i", "--index",
                      help="Use the supplied bowtie index file, instead of indexing the genome in the script. \
                      If your genome index files are in folder /user/home/index, with names TAIR10.1.ebwt, TAIR10.2.ebwt, etc, you must use -i /user/home/index/TAIR10 for the -i option.")
    parser.add_option("-t", "--temp", help="Temporary folder to hold the bowtie index files. If not supplied, the current directory is used. Only used with -r.")
    parser.add_option("-v", "--allowedmismatch", type=int, default=0, help="-v option in bowtie. Number of mismatches allowed. Default is 0.")
    parser.add_option("-k", "--multialignment", type=int, default=20, help="-k option in bowtie.  Report up to <int> vaild alignment. Default is 20.")
    parser.add_option("-p", "--processor", type=int, default=1,
                      help="Use multiple threads to do alignment.")
    parser.add_option("-f", "--filterunmapped", action="store_true", help="Filter out unmapped alignments in the output.")

    (options, args) = parser.parse_args()
    if len(args) < 1:
         parser.error("incorrect number of arguments. Run the script with -h option to see help.")

    if options.reference and options.index:
        parser.error("Options -r and -i are mutually exclusive. Please use only one of them.")

    if not options.reference and not options.index:
        parser.error("Either option -r or -i should be provided.")

    if options.index and options.temp:
        parser.error("Option -t is not needed for option '-r'")

    if options.reference:
        for name in options.reference:
            absname = os.path.abspath(os.path.expanduser(name))
            if not os.path.exists(absname):
                parser.error("File " + name + " in option -r does not exist!!")

    if options.index:
        suffix = ['1.ebwt','2.ebwt','3.ebwt','4.ebwt','rev.1.ebwt','rev.2.ebwt']
        for s in suffix:
            name = options.index + "." + s
            absname = os.path.abspath(os.path.expanduser(name))
            if not os.path.exists(absname):
                parser.error("Index file " + name + " does not exist!! Please use the -r option instead.")

    for name in args:
        absname = os.path.abspath(os.path.expanduser(name))
        if not os.path.exists(absname):
            parser.error("File " + name + " does not exist!!")

    return options, args


def check_Bowtie():
    ret = distutils.spawn.find_executable("bowtie")
    if ret:
        return True
    else:
        return False

def check_Bowtiebuild():
    ret = distutils.spawn.find_executable("bowtie-build")
    if ret:
        return True
    else:
        return False


def check_samtools():
    ret = distutils.spawn.find_executable("samtools")
    if ret: #check the depth, faidx, view command exist.
        return True
    return False


def check_readid(readname):
    """Only check the first 2000 sequences. If all the first 2000 sequences have
    the right ID formart, then the script assumes all the read IDs are in the
    right format.

    """
    count = 0
    pattern=r"^>\S+_r[0-9]+_x[0-9]+$"
    allgood = True
    with open(readname) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
                if count>= 2000:
                    break
                m = re.match(pattern, line.strip())
                if not m:  # seqid right
                    allgood = False
                    break
    return allgood


def run_bowtie_build(fastalist, tempfolder):
    command = "bowtie-build -f "
    for fasta in fastalist:
        command = command + fasta + " "

    command = command + " " + os.path.join(tempfolder, "genomeindex")
    sys.stdout.write("\nIndexing reference genomes, command:\n")
    sys.stdout.write(command+'\n')
    outerr = ""
    try:
        build_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = build_process.communicate()
        return os.path.join(tempfolder, "genomeindex")
    except Exception:
        sys.stderr.write("Error occurred when indexing reference sequences.\n")
        if outerr:
            sys.stderr.write(outerr+'\n')
        os.exit(-1)


def run_bowite(fastaname, outname, indexbase, options):
    #bowtie -a -v 0 -p 8 -m 30 -S bowtie-index/TAIR10  -f SAMPLE1.fasta.processed > SAMPLE1.sam 2> SMAPLE1.log
    command = "bowtie -p " + str(options.processor)
    command = command + " -v " + str(options.allowedmismatch)
    command = command + " --best --strata -k " + str(options.multialignment) + " -S "
    command = command + indexbase + " " + " -f " + fastaname
    fout = open(outname, 'w')
    try:
        sys.stdout.write("\nMapping file " + fastaname + ", command:\n")
        sys.stdout.write(command+"\n")
        bowtie_process = subprocess.check_call(command.split(), stdout=fout)
        return outname
    except Exception:
        sys.stderr.write("Error occurred when mapping reads.\n")
        sys.stderr.write("Command that was running: " + command + "\n")
        sys.exit(-1)


def filter_unmapped_alignment(samname, newsamname):
    command = "samtools view -Sh -F 4 " + samname + " -o " + newsamname
    try:
        sys.stdout.write("\nFiltering unmapped alignments, command:\n")
        sys.stdout.write(command+"\n")
        filter_process = subprocess.check_call(command.split())
        return newsamname
    except Exception:
        sys.stderr.write("Error occurred when filtering unmapped alignments")
        sys.exit(-1)

if __name__ == '__main__':
    options, args = parse_option_optparse()
    if not check_Bowtie():
        sys.stderr.write("bowtie is required but not installed or not in the PATH.\n")
        sys.exit(-1)
    if not check_Bowtiebuild():
        sys.stderr.write("bowtie-build command is required but not in the PATH.\n")
        sys.exit(-1)
    if options.filterunmapped:
        if not check_samtools():
            sys.stderr.write("SAMtools is required for the -f option but not in the PATH.\n")
            sys.exit(-1)
    # check the read ID.
    allgood = True
    for name in args:
        if not check_readid(name):
            sys.stderr.write("ERROR: The format of the read IDs in file " + name +" is not right.\n\n")

            sys.stderr.write("The ID for each read in the output files has the following format: SampleName_rA_xN. Here A is an number uniquely identifies the read, N is the total count of the read (depth of the read). 'SampleName' is the name of the sample/tissue/library that the input fasta file represents. For example, if the sample name is 'root', and the read occurred 120 times in the library, then it's identifier could be root_r23_x120. Here '23' is just a number that donates the order of the reads when processing it, it has no use otherwise.\n")
            sys.stderr.write("Please use the provided scripts(process-reads-fasta.py, convert-mirdeep2-fasta.py, convert-readcount-file.py) to preprocess the read fasta files. Refer to the README file for more information.\n\n")
            allgood = False

    if not allgood:
        exit(-1)
    indexbase = ""
    if not options.index:
        indexfolder = ""
        if options.temp:
            tempfolder = os.path.expanduser(options.temp)
            if not os.path.exists(tempfolder):
                os.makedirs(tempfolder)
            indexfolder = tempfolder
        else:
            indexfolder = tempfile.mkdtemp(prefix='bowtie-index')
        run_bowtie_build(options.reference, indexfolder)
        indexbase = os.path.join(indexfolder, "genomeindex")
    else:
        indexbase = options.index

    outnames = []
    for name in args:
        finaloutname = name+".sam"
        if options.filterunmapped:
            outf = tempfile.NamedTemporaryFile(delete=False)
            run_bowite(name, outf.name, indexbase, options)
            outnames.append(filter_unmapped_alignment(outf.name, finaloutname))
            os.unlink(outf.name)
        else:
            outnames.append(run_bowite(name, finaloutname, indexbase, options))

    sys.stdout.write("===============================================\nDONE\n")
    sys.stdout.write("Output SAM files can be found at:\n")
    for name in outnames:
        sys.stdout.write(name + "\n")
