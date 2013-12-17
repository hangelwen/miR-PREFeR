import sys
import os
import multiprocessing

'''Pre-process a list of uncollapsed fasta files. Collapse the identical reads and
rename the reads. The new identifier for each read in the output files has the
following format: SampleName_rA_xB. Here A is an number identifies the read,
which an sequentially increasing number assigned to the read by the script. B is
the total count of the read (depth of the read). 'SampleName' is the name of the
sample/tissue/library that the input fasta file represents.

Input:
1. samplenamlist: a file in which the sample names are listed. Each
sample name should be on one line. Sample name should be one word without any
blank characters. The number of sample names in the file should be the same as
the number of input fasta files. Each name should be unique.

2. list of fasta files: a list of uncollapsed fasta files, each represents a
RNAseq dataset from a sample/library/tissue/etc. The order of the files should
be consistent with the names listed in the samplenamlist file.

For example, if the sample name is 'root', and the read occurred 120 times in
the library, then it's identifier could be root_r23_x120. Here '23' is just a
number that donates the order of the reads when processing it, it has no use
otherwise.

The script assumes each read sequence occupies ONLY one line in the input file.

The output files are in the same folder as the input files, with ".processed "
as suffix.

'''
sys.stderr.write("========================================================================================= \n")
sys.stderr.write("Pre-processing RNAseq fasta files\n")
sys.stderr.write("Please make sure the input fasta files are uncollpased fasta files.\n")
sys.stderr.write("Please make sure the sequence of each read only has one line in the input file.\n")
sys.stderr.write("The output files are in the same folders as the input files, with suffix '.processed'. \n")
sys.stderr.write("========================================================================================= \n")
if len(sys.argv) < 3:
    sys.stderr.write("Usage: python preocess-reads-fasta.py <samplenamelist> <fasta1> <fasta2> .. <fastaN>\n")
    sys.exit(-1)


prefix = []
with open(sys.argv[1]) as f:
    for line in f:
        if line.strip():
            prefix.append(line.strip())

if len(prefix) != len(sys.argv[2:]):
    sys.stderr.write("Error: number of sample/tissue names in the samplenamelist file must be the same as the number of input fasta files.\n")
    sys.exit(-1)

for idx, name in enumerate(sys.argv[2:]):
    if not os.path.exists(name):
        sys.stderr.write("Error: file " +name + " does not exist!!!\n")
        sys.exit(-1)


def process_one_file(name, prefix, queue):
    sys.stdout.write("Start processing file "+name+"\n")
    dict_reads = {}
    outf = open(name + ".processed", 'w')
    with open(name) as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                key = line.strip()
                if key in dict_reads:
                    dict_reads[key] += 1
                else:
                    dict_reads[key] = 1
    for cnt, r in enumerate(dict_reads):
        identifier = ">"+prefix+"_r"+str(cnt)+"_x"+str(dict_reads[r])
        outf.write(identifier+"\n")
        outf.write(r+"\n")
    outf.close()
    sys.stdout.write("Finish file "+name+"\n")
    queue.put((name, len(dict_reads)))

jobs = []
inforqueue = multiprocessing.Queue()
for idx, name in enumerate(sys.argv[2:]):
    p = multiprocessing.Process(target=process_one_file, args=(name, prefix[idx]
                                                               , inforqueue))
    p.start()
    jobs.append(p)

for job in jobs:
    job.join()
    info = inforqueue.get()
    sys.stdout.write("File " + info[0]+ " has "+ str(info[1]) + " unique reads\n")
sys.stdout.write("DONE\n\n")
