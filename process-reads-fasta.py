import sys
import os
import multiprocessing
import gc

'''
Pre-process a list of uncollpased fasta files. Collapse the identical reads and
rename the reads. The identifier for each read in the output files are:
SN-A-B. Here N is the order of the input fasta file, A is an number identify the
read, B is the total count of the read. For example, if in the second input
file, a read occurred 1200 times, then it's description line looks like:
>S2-23-1200. Here '23' is just a number that donates the order of the reads when
processing it, it has no use otherwise.

The script assumes each read sequence occupies ONLY one line in the input file.

The output files are in the same folder as the input files, with ".processed "
as suffix.
'''

sys.stderr.write("Pre-processing RNAseq fasta files\n")
sys.stderr.write("Please make sure the input fasta files are uncollpased fasta files.\n")
if len(sys.argv) == 1:
    sys.stderr.write("Usage: python preocess-reads-fasta.py <fasta1> <fasta2> .. <fastaN>\n")
    sys.exit(-1)



prefix = []
for idx, name in enumerate(sys.argv[1:]):
    if not os.path.exists(name):
        sys.stderr.write("Error: file " +name + " does not exist!!!\n")
        sys.exit(-1)
    prefix.append("S"+str(idx+1))


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
        identifier = ">"+prefix+"-"+str(cnt)+"-"+str(dict_reads[r])
        outf.write(identifier+"\n")
        outf.write(r+"\n")
    outf.close()
    sys.stdout.write("Finish file "+name+"\n")
    queue.put((name, len(dict_reads)))

jobs = []
inforqueue = multiprocessing.Queue()
for idx, name in enumerate(sys.argv[1:]):
    p = multiprocessing.Process(target=process_one_file, args=(name, prefix[idx]
                                                               , inforqueue))
    p.start()
    jobs.append(p)

for job in jobs:
    job.join()
    info = inforqueue.get()
    sys.stdout.write("File " + info[0]+ " has "+ str(info[1]) + " unique reads\n")
sys.stdout.write("DONE\n\n")
