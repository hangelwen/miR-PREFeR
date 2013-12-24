import sys
import os
import multiprocessing

'''
Pre-precess a list of read count files. In a read count file, each line contains
two columns: the first column is a read sequence, the second column is the depth
of the read. The two columns can be seperated by spaces or tabs.  The new
identifier for each read in the output files has the following format:
SampleName_rA_xB. Here A is an number identifies the read, which an sequentially
increasing number assigned to the read by the script. B is the total count of
the read (depth of the read). 'SampleName' is the name of the
sample/tissue/library that the input fasta file represents. For more
information, see the doc in the process-reads-fasta.py

The output files are in the same folder as the input files, with ".processed "
as suffix.

'''

sys.stderr.write("========================================================================================= \n")
sys.stderr.write("Pre-processing read-count format files.\n")
sys.stderr.write("The output files are in the same folders as the input files, with suffix '.processed'. \n")
sys.stderr.write("========================================================================================= \n")

if len(sys.argv) < 3:
    sys.stderr.write("Usage: python process-reads-fasta.py <samplenamelist> <rc1> <rc2> .. <rcN>\n")
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
            if not line.strip():
                continue
            sp = line.split()
            seqid = ">"+prefix+"_r"+str(count)+"_x"+sp[-1].strip()
            outf.write(seqid+'\n')
            outf.write(sp[0]+'\n')
            count += 1
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
