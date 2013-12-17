import sys
import os
import multiprocessing

'''

Pre-precess a list of collapsed fasta format read files. The input read files
following the miRDeep2 naming convention: the identifier of each read ends with
"_xN", where 'N' is the depth of the read in the library. The new identifier for
each read in the output files has the following format: SampleName_rA_xB. Here A
is an number identifies the read, which an sequentially increasing number
assigned to the read by the script. B is the total count of the read (depth of
the read). 'SampleName' is the name of the sample/tissue/library that the input
fasta file represents. For more information, see the doc in the
process-reads-fasta.py

The script assumes each read sequence occupies ONLY one line in the input file.

The output files are in the same folder as the input files, with ".processed "
as suffix.

'''

sys.stderr.write("========================================================================================= \n")
sys.stderr.write("Pre-processing RNAseq fasta files (mirdeep2 format).\n")
sys.stderr.write("Please make sure the sequence of each read only has one line in the input file.\n")
sys.stderr.write("The output files are in the same folders as the input files, with suffix '.processed'. \n")
sys.stderr.write("========================================================================================= \n")

if len(sys.argv) < 3:
    sys.stderr.write("Usage: python process-reads-fasta.py <samplenamelist> <fasta1> <fasta2> .. <fastaN>\n")
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
    outf = open(name + ".processed", 'w')
    count = 1
    with open(name) as f:
        for line in f:
            if line.startswith(">"):
                sp = line.split("x")
                seqid = ">"+prefix+"_r"+str(count)+"_x"+sp[-1]
                outf.write(seqid)
                count += 1
            else:
                outf.write(line)
    outf.close()
    sys.stdout.write("Finish file "+name+"\n")

jobs = []
inforqueue = multiprocessing.Queue()
for idx, name in enumerate(sys.argv[2:]):
    p = multiprocessing.Process(target=process_one_file, args=(name, prefix[idx]
                                                               , inforqueue))
    p.start()
    jobs.append(p)

for job in jobs:
    job.join()

sys.stdout.write("DONE\n\n")
