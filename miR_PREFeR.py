import sys
import os.path
import multiprocessing
import Queue
import cPickle
import time
import re
import imp
import subprocess
import os
import distutils.spawn
import tempfile
import string
import gc
import logging
import shutil

dict_failure_reasons = {1: "Structure is not stemloop",
                       2: "Mature-star duplex has too many bugles/interior loops",
                       3: "Mature-star duplex has large bugles/interior loops",
                       4: "Star sequence not expression in the RNA-seq data",
                       5: "Expression pattern not good"
}

def parse_option():
    import argparse
    helpstr = """    check = Check the dependency and the config file only (default).
    prepare = Prepare data.
    candidate = Generate candidate regions.
    fold = Fold the candidate regions.
    predict = Predict miRNAs.
    pipeline = Run the whole pipeline. This is the same as running 'check', 'prepare', 'candidate', 'fold', 'predict' sequentially.
    recover = Recover a unfinished job. By default, miR-PREFeR makes checkpoint of the results of each stage. Thus, an unfinished job can be started from where it was checkpointed to save time.
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("action",
                        choices=('check','prepare','candidate','fold','predict','pipeline', 'recover'),
                        default = 'check', help = helpstr)
    parser.add_argument("configfile",
                        help="    Configure file.")
    args = parser.parse_args()
    action = args.action
    configfile = args.configfile

    dict_option = parse_configfile(configfile)
    dict_option['ACTION'] = action
    return dict_option

def parse_option_optparse():
    from optparse import OptionParser
    helpstr = """python mir_PREFeR.py [options] command configfile

command could be one of the following:
    check = Check the dependency and the config file only (default).
    prepare = Prepare data.
    candidate = Generate candidate regions.
    fold = Fold the candidate regions.
    predict = Predict miRNAs.
    pipeline = Run the whole pipeline. This is the same as running 'check', 'prepare', 'candidate', 'fold', and 'predict' sequentially.
    recover = Recover a unfinished job. By default, miR-PREFeR makes checkpoint of the results of each stage. Thus, an unfinished job can be started from where it was checkpointed to save time.

configfile: configuration file"""

    parser = OptionParser(helpstr)
    parser.add_option("-l", "--log", action="store_true", dest="log",
                      help="Generate a log file.")
    parser.add_option("-k", "--keep-tmp", action="store_true", dest="keeptmp",
                      help="After finish the whole pipeline, do not remove the temporary folder that contains the intermidate files.")
    actions = ['check','prepare','candidate','fold', 'predict', 'pipeline', 'recover']
    (options, args) = parser.parse_args()
    if len(args) != 2:
         parser.error("incorrect number of arguments. Run the script with -h option to see help.")
    if args[0] not in actions:
        parser.error("unknow command")
    dict_option = parse_configfile(args[1])
    dict_option['ACTION'] = args[0]
    dict_option["LOG"] = options.log
    dict_option["KEEPTMP"] = options.keeptmp
    return dict_option

def parse_configfile(configfile):
    if not os.path.exists(configfile):
        sys.stderr.write("Configuration file " + configfile + " does not exist!!\n")
        sys.exit(-1)
    dict_option = {
        "CONFIG_FILE":configfile,
        "FASTA_FILE":"",
        "ALIGNMENT_FILE":[],
        "GFF_FILE":"",
        "PRECURSOR_LEN":300,
        "READS_DEPTH_CUTOFF":20,
        "MAX_GAP":100,
        "NUM_OF_CORE":1,
        "OUTFOLDER":"./",
        "NAME_PREFIX":"",
        "PIPELINE_PATH":"",
        "ONLY_SEQ":[],
        "DELETE_IF_SUCCESS":"Y"
    }
    with open(configfile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            sp = line.strip().split("=")
            if len(sp)>1 and sp[1]:
                key = sp[0].strip()
                if key == "ALIGNMENT_FILE":
                    names = sp[1].split(",")
                    for name in names:
                        if not os.path.exists(name.strip()):
                            sys.stderr.write("File " + name.strip() +
                                             " does not exist!!\n")
                            sys.exit(-1)
                        dict_option[key].append(name.strip())
                    continue
                if key == "ONLY_SEQ":
                    names = sp[1].split(",")
                    for name in names:
                        dict_option[key].append(name.strip())
                if key == "GFF_FILE" or key == "FASTA_FILE":
                    if not os.path.exists(sp[1].strip()):
                        sys.stderr.write("File " + sp[1].strip() +
                                          " does not exist!!\n")
                        sys.exit(-1)
                    dict_option[key] = sp[1].strip()
                    continue
                if key == "NUM_OF_CORE":
                    cpucount = multiprocessing.cpu_count()
                    cpu_to_use = int(sp[1].strip())
                    if cpucount < cpu_to_use:
                        sys.stderr.write(
                            "Warnning: NUM_OF_CORE is larger than CPUS/Cores on" +
                            " the machine. Use "+str(cpucount)+" instead.\n")
                        cpu_to_use = cpucount
                    dict_option[key] = cpu_to_use
                    continue
                if key == "PRECURSOR_LEN" or key == "READS_DEPTH_CUTOFF" or key == 'MAX_GAP':
                    dict_option[key] = int(sp[1].strip())
                    continue
                if key =="OUTFOLDER" or key == "NAME_PREFIX" or key == "DELETE_IF_SUCCESS":
                    dict_option[key] = sp[1].strip()
                    continue
                if key == "PIPELINE_PATH":
                    if not os.path.exists(sp[1].strip()):
                        sys.stderr.write("miR-PREFeR path " + sp[1].strip() +
                        " does not exist!!\n")
                        sys.exit(-1)
                    dict_option[key] = sp[1].strip()
                    continue
    return dict_option


def get_current_local_time():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())


def write_formatted_string(message, leftpad_len, f, padchar=""):
    outstr = "{0:"+padchar+"<"+str(leftpad_len)+"}"+message
    f.write(outstr.format(""))
    f.write("\n")
    f.flush()


def write_formatted_string_withtime(message, leftpad_len, f, padchar=""):
    t = get_current_local_time()
    outstr = "{0:"+padchar+"<"+str(leftpad_len)+"}"+message
    f.write(outstr.format(t))
    f.write("\n")
    f.flush()


def write_gff_line(seqid, start, end, strand, ID, name, score=".", source="miR-PREFeR",
                   feature="miRNA", frame=".", other = "", fout=sys.stdout):
    '''
    Write a gff line to fout.
    The line is:
    seqid source feature start end score strand frame ID=id;NAME=name;Other=other
    '''
    outstr = "\t".join([seqid, source, feature, str(start), str(end), str(score),
                        strand, frame, "ID="+ID+";NAME="+name+";Other="+str(other)])
    fout.write(outstr)
    fout.write("\n")


def get_complement(seq):
    #TODO if the sequence contains other characters, i.e. IUPAC?
    trans_table = string.maketrans("ATGCU","UACGA")
    return seq.translate(trans_table)


def get_reverse_complement(seq):
    return get_complement(seq)[::-1]


def compute_RPKM(transcript_len, reads_in_transcript, reads_total):
    return float(reads_in_transcript)*1000000000/(reads_total * transcript_len)


def check_Bowtie():
    ret = distutils.spawn.find_executable("bowtie")
    if ret:
        return True
    else:
        return False


def check_samtools():
    ret = distutils.spawn.find_executable("samtools")
    if ret:
        return True
    else:
        return False


def check_RNALfold():
    ret = distutils.spawn.find_executable("RNALfold")
    if ret:
        return True
    else:
        return False

def get_RNALfold_version():
    command = "RNALfold -V"
    try:
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
    except Exception as e:
        return "UnknowVersion"
    if outmessage == "":
        return "UnknowVersion"
    else:
        return str(outmessage).strip()

def is_bug_RNALfold(version):
    '''
    RNALfold in Vienna 2.0.4 has a bug: when the input sequence has no valid
    secondary structure, it produces a segmentation fault and can not continue.
    '''
    if version.strip() == "RNALfold 2.0.4":
        return True
    return False

def get_length_from_sam(samfile):
    dict_len = {}
    with open(samfile) as f:
        for line in f:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    sp = line.split()
                    dict_len[sp[1].split(":")[1]] = int(sp[2].split(":")[1])
            else:
                return dict_len


def index_genome(fastaname):
    '''
    Run 'samtools faidx fastafile' command to index the genome fasta file.
    '''
    command = "samtools faidx "+fastaname
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when indexing the genome file\n")
        sys.exit(-1)


def gen_keep_regions_from_gff(gffname, tmpdir, dict_len, minlen):
    '''
    Generate a BED format file which contains all the regions that are not
    overlap with the regions of the features in the gff file.

    Note that gff file is 0 based, and the ending position is inclusive. BED
    file is 1 based, and the ending position is exclusive.
    '''
    def overlap(r1, r2):
        '''
        If r1 and r2 overlaps, return the combined region, else, return None.
        '''
        if r1[1] < r2[0] or r1[0] > r2[1]:
            return None
        else:
            return (min(r1[0],r2[0]), max(r1[1],r2[1]))

    #sort the gff file
    tempgffname = os.path.join(tmpdir, "temp.remove.gff")
    tempbedname = os.path.join(tmpdir, "temp.keepregion.bed")
    command = "sort -k1,1 -k4,4n " + gffname + " -o " + tempgffname
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when sorting the gff file\n")
        sys.exit(-1)
    foutput = open(tempbedname, 'w')
    seqids = dict_len.keys()
    with open(tempgffname) as f:
        seqid = "NOTKNOW"
        region = None
        for line in f:
            if line.startswith("#"):
                continue
            sp = line.split()
            if sp[0] not in seqids:
                continue
            cur_seqid = sp[0]
            if cur_seqid != seqid:
                if seqid != "NOTKNOW":
                    if dict_len[seqid] - region[1] >= minlen:
                        foutput.write(seqid+"\t"+str(region[1]+1)+"\t"+str(dict_len[seqid])+"\n")
                if int(sp[3]) >=minlen:
                    foutput.write(cur_seqid+"\t1"+"\t"+sp[3]+"\n")
                seqid = cur_seqid
                region = (int(sp[3]), int(sp[4]))
            else:
                ret = overlap(region, (int(sp[3]), int(sp[4])))
                if ret:
                    region = ret
                else:
                    if int(sp[3])-region[1] >=minlen:
                        foutput.write(cur_seqid+"\t"+str(region[1]+1)+"\t"+sp[3]+"\n")
                    region = (int(sp[3]), int(sp[4]))
                    seqid = cur_seqid
            #write the last one
        if dict_len[seqid] - region[1] >= minlen:
            foutput.write(seqid+"\t"+str(region[1]+1)+"\t"+str(dict_len[seqid])+"\n")
        foutput.close()
    return tempbedname


def sam2bam(samfile, bamfile):
    command = "samtools view -bS -o" + bamfile + " " + samfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when converting SAM to BAM\n")
        sys.exit(-1)


def combine_bamfiles(headersamfile, outbamfile, *bamfiles):
    '''
    Combine BAM files sample as one file.
    '''
    command = "samtools cat -h " + headersamfile + " -o " + outbamfile
    for bamname in bamfiles:
        command = command + " "+bamname
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when combining BAM files\n")
        sys.exit(-1)
    return outbamfile


def gen_keep_regions_sort_bam(bamfile, bedfile, outbamprefix):
    '''
    Only keep alignments that overlap with the regions specified in the BED
    file, and then sort the bam file.

    To achieve this, the following samtools command can be used:
    samtools view -L bedfile -F 4 bamfile -b -o outbamfile
    '''
    #TODO samtools view generate reads overlap with the region, not only reads
    #in the region. Should this be a problem here??
    command1 = "samtools view -L " + bedfile + " -F 4 -b -o " + outbamprefix+".bam" + " " +bamfile
    command2 = "samtools sort "+ outbamprefix+".bam " + outbamprefix + ".sort"
    command3 = "samtools index " + outbamprefix+".sort.bam"
    try:
        subprocess.check_call(command1.split())
        subprocess.check_call(command2.split())
        subprocess.check_call(command3.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when filtering regions and sorting the BAM file\n")
        sys.exit(-1)
    return outbamprefix+".sort.bam"

def sort_index_bam(bamfile, outbamprefix):
    command2 = "samtools sort "+ bamfile + " " + outbamprefix + ".sort"
    command3 = "samtools index " + outbamprefix+".sort.bam"
    try:
        subprocess.check_call(command2.split())
        subprocess.check_call(command3.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when filtering regions and sorting the BAM file\n")
        sys.exit(-1)
    return outbamprefix+".sort.bam"


def expand_bamfile(bamfile, maxdepth, outputbamfile, outputsamfile):
    tempsamfile = tempfile.NamedTemporaryFile(mode='w',prefix = str(os.getpid())+"tempsam", suffix=".sam", delete=False)
    command = "samtools view -h -o " + tempsamfile.name + " " + bamfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when converting BAM to SAM\n")
        sys.exit(-1)
    tempsamfile.close()
    outf = open(outputsamfile, 'w')
    with open(tempsamfile.name) as f:
        for line in f:
            if line.startswith("@"):
                outf.write(line)
                continue
            sp = line.split()
            depth = int(sp[0].split("-")[-1])
            if depth > maxdepth:
                depth = maxdepth
            for i in xrange(depth):
                outf.write(line)
        outf.close()
    command = "samtools view -bS -o" + outputbamfile + " " + outputsamfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when converting SAM to BAM\n")
        sys.exit(-1)
    return outputbamfile, outputsamfile


def index_bam(bamfile, outindexfile):
    command = "samtools index " + bamfile + " "+outindexfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when indexing the BAM file\n")
        sys.exit(-1)
    return outindexfile


def filter_bam_by_flag(bamfile, flag, outbamname, keep=True):
    keepflag = " -F " + str(flag)
    if keep:
        keepflag = " -f " + str(flag)
    command = "samtools view -b "+keepflag + " -o " +outbamname + " " + bamfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when filtering BAM file using flags.\n")
        sys.exit(-1)
    return outbamname


def prepare_data(dict_option, outtempfolder, logger):

    if logger:
        logger.info("Getting genomic sequence lengths.")
    #get the length of all the genomic sequences in the fasta/alignment files
    dict_len = get_length_from_sam(dict_option["ALIGNMENT_FILE"][0])

    #if the fasta file is not index, then create a index
    if not os.path.exists(dict_option["FASTA_FILE"]): #index the genome
        index_command = "samtools faidx " + dict_option["FASTA_FILE"]
        try:
            subprocess.check_call(index_command.split())
        except Exception as e:
            if logger:
                logger.error("Error occurred when indexing the genome file. "+
                             "Command: "+index_command)
            sys.stderr.write(e)
            sys.stderr.write("Error occurred when indexing the genome file\n")
            sys.exit(-1)

    #convert the SAM files(generated by Bowtie) to BAM files
    if logger:
        logger.info("Converting SAM files to BAM files using samtools.")
    bamfiles = []
    for samname in dict_option["ALIGNMENT_FILE"]:
        if logger:
            logger.info("Converting "+samname)
        bamname = os.path.join(outtempfolder, os.path.basename(samname)+".bam")
        sam2bam(samname, bamname)
        bamfiles.append(bamname)

    #combine multiple BAM files from multiple sample together
    if logger:
        logger.info("Combining multiple BAM files from multiple samples together.")
    combinedbamname = os.path.join(outtempfolder, "combined.bam")
    combine_bamfiles(dict_option["ALIGNMENT_FILE"][0], combinedbamname, *bamfiles)

    #removing reads that are overlapped with features in the gff file, if provided.
    if os.path.exists(dict_option["GFF_FILE"]):
        #TODO: minlen should be user adjustable, not fix here.
        tempkeepregion = gen_keep_regions_from_gff(dict_option["GFF_FILE"], outtempfolder, dict_len, 60)
        if logger:
            logger.info("Removing reads that are overlapped with features in the gff file.")
        combinedbamname = gen_keep_regions_sort_bam(combinedbamname, tempkeepregion, os.path.join(outtempfolder,"combined.filtered"))
    else:
        combinedbamname = sort_index_bam(combinedbamname, os.path.join(outtempfolder,"combined.filtered"))

    expandedsamname = os.path.join(outtempfolder, "expanded.sam")
    expandedbamname = os.path.join(outtempfolder, "expanded.bam")
    expandedbam_plus = os.path.join(outtempfolder, "expanded.plus.bam")
    expandedbam_minus = os.path.join(outtempfolder, "expanded.minus.bam")
    if logger:
        logger.info("Generating expanded BAM and SAM files")
    expand_bamfile(combinedbamname, dict_option["READS_DEPTH_CUTOFF"], expandedbamname, expandedsamname)
    if logger:
        logger.info("Generating "+expandedbam_plus)
    filter_bam_by_flag(expandedbamname, 16, expandedbam_plus, keep=False)
    if logger:
        logger.info("Generating "+expandedbam_minus)
    filter_bam_by_flag(expandedbamname, 16, expandedbam_minus, keep=True)
    return combinedbamname, expandedsamname, expandedbamname, expandedbam_plus, expandedbam_minus


def gen_contig_typeA(expandedbam_plus, expandedbam_minus, dict_option,
                     contig_minlen, logger):
    def get_next_non_zero_region(name):
        with open(name) as f:
            region = []
            line = f.readline()
            sp = line.split()
            seperate_depth = [int(x) for x in sp[2:]]
            depth = sum(seperate_depth)
            pre_seqid = sp[0]
            pre_pos = int(sp[1])
            start_pos = int(sp[1])

            plus_depth = seperate_depth[0]
            minus_depth = seperate_depth[1]

            for line in f:
                sp = line.split()
                cur_seqid = sp[0]
                pos = int(sp[1])
                seperate_depth = [int(x) for x in sp[2:]]
                depth = sum(seperate_depth)

                if cur_seqid != pre_seqid:
                    if plus_depth > minus_depth:
                        yield [pre_seqid, start_pos, pre_pos+1, "+"]
                    else:
                        yield [pre_seqid, start_pos, pre_pos+1, "-"]
                    plus_depth = seperate_depth[0]
                    minus_depth = seperate_depth[1]
                    start_pos = pos
                    pre_pos = pos
                    pre_seqid = cur_seqid
                    region = []

                if pos - pre_pos > 1:
                    region.append(cur_seqid)
                    region.append(start_pos)
                    region.append(pre_pos+1)
                    if plus_depth > minus_depth:
                        region.append("+")
                    else:
                        region.append("-")
                    yield region
                    start_pos = pos
                    pre_pos = pos
                    region = []
                    plus_depth = seperate_depth[0]
                    minus_depth = seperate_depth[1]
                else:
                    pre_pos = pos
                    plus_depth = plus_depth + seperate_depth[0]
                    minus_depth = minus_depth + seperate_depth[1]
            if plus_depth > minus_depth:
                yield [cur_seqid, start_pos, pre_pos+1, "+"]
            else:
                yield [cur_seqid, start_pos, pre_pos+1, "-"]
    try:
        if logger:
            logger.info("Generating the depth file using samtools.")
        samtools_process = subprocess.Popen(["samtools","depth",expandedbam_plus, expandedbam_minus], stdout=subprocess.PIPE)
        awk_rule = "$3+$4>"+str(dict_option["READS_DEPTH_CUTOFF"])
        awk_process = subprocess.Popen(["awk",awk_rule], stdout=subprocess.PIPE, stdin=samtools_process.stdout)
        samtools_process.stdout.close()
        output = awk_process.communicate()[0]
    except Exception as e:
            sys.stderr.write(e)
            sys.stderr.write("Error occurred when generating contigs(typeA).\n")
            sys.exit(-1)
    depthfilename = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_tmp", "bam.depth.cut"+str(dict_option["READS_DEPTH_CUTOFF"]))
    f=open(depthfilename,'w')
    f.write(output.decode())
    f.close()
    if logger:
        logger.info("Generating contigs.")
    dict_contigs = {}
    for region in get_next_non_zero_region(depthfilename):
        if region[2] - region[1] < contig_minlen:
            continue
        dict_contigs.setdefault(region[0], []).append((region[1],region[2],region[3]))
    return depthfilename, dict_contigs


def dump_loci_seqs_samtool(dict_loci, fastaname, outputprefix, num_of_proc):
    '''
    Generate num_of_proc fasta files that contains all the loci/extended region sequence.

    Return the file names and the number of sequences in each file.
    '''

    seqids = sorted(dict_loci.keys())
    num_seq = 0
    total_loci = 0
    for seqid in seqids:
        total_loci += len(dict_loci[seqid])
    # How many loci in one file?
    num_in_one_file = int(total_loci/num_of_proc)
    cur_loci = 0
    cur_fout = 0
    fname = outputprefix+"_0.fa"
    fout  = open(fname, 'w')
    ret_list = []  # a list of (name, num_of_seq) tuples
    for seqid in seqids:
        for loci in dict_loci[seqid]:
            if cur_loci < num_in_one_file * (cur_fout+1):
                pass
            else:
                cur_fout += 1
                fout.close()
                ret_list.append((fname, num_seq))
                num_seq = 0
                fname = outputprefix+"_"+str(cur_fout)+".fa"
                fout = open(outputprefix+"_"+str(cur_fout)+".fa", 'w')
            cur_loci += 1
            plus_seq = []
            minus_seq = []
            plus_tag = []
            minus_tag = []
            loci_pos = str(loci[0][0]) + "-" + str(loci[0][1])
            for idx, extendregion in enumerate(loci[1:]):
                region = seqid+":"+str(extendregion[0][0])+"-"+str(extendregion[0][1]-1)
                command = "samtools faidx "+fastaname+" "+region
                seq = ""
                try:
                    samtools_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    outmessage, outerr = samtools_process.communicate()
                    seq = str("".join(outmessage.decode().split("\n")[1:]))
                except Exception as e:
                    sys.stderr.write(e)
                    sys.stderr.write("Error occurred when cutting loci sequence.\n")
                    sys.exit(-1)
                strands = set()
                otherinfo = ""
                for peak in extendregion[1]:
                    strands.add(peak[2])
                    otherinfo = otherinfo + str(peak[0])+","+str(peak[1])+","+str(peak[2])+";"
                otherinfo = otherinfo.rstrip(";")
                seqtag = ""
                # a tag indicates which extend region this sequence belongs to.
                # Possible values:
                # "0": This loci has only one extend region, and 'seq' is the sequence for the region
                # "L": This loci has two extend regions, and 'seq' is the sequence for the left region
                # "R": This loci has two extend regions, and 'seq' is the sequence for the right region
                extendregiontag = " 0 "
                if len(loci)==3 and idx == 0:  # Left side region
                    extendregiontag = " L "
                elif len(loci)==3 and idx == 1:  # Right side region
                    extendregiontag = " R "
                if len(strands) == 1:
                    if "+" in strands:
                        seqtag = ">"+region+" + "+ loci_pos + extendregiontag + otherinfo+"\n"  # >ExtendRegion strand loci 0/L/R peaks
                    else:
                        seqtag = ">"+region+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                        seq = get_reverse_complement(seq)
                    fout.write(seqtag)
                    fout.write(seq+"\n")
                    num_seq += 1
                    continue
                if len(strands) == 2:
                    seqtag = ">"+region+" + "+ loci_pos + extendregiontag + otherinfo+"\n"
                    plus_tag.append(seqtag)
                    plus_seq.append(seq)
                    #fout.write(seqtag)
                    #fout.write(seq+"\n")
                    num_seq += 1
                    seq = get_reverse_complement(seq)
                    seqtag = ">"+region+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                    minus_tag.append(seqtag)
                    minus_seq.append(seq)
                    #fout.write(seqtag)
                    #fout.write(seq+"\n")
                    num_seq += 1
            if len(plus_tag):
                for which, tag in enumerate(plus_tag):
                    fout.write(tag)
                    fout.write(plus_seq[which]+"\n")
                for which, tag in enumerate(minus_tag):
                    fout.write(minus_tag[which])
                    fout.write(minus_seq[which]+"\n")
    ret_list.append((fname, num_seq))
    fout.close()
    return ret_list

def dump_loci_seqs_and_alignment(dict_loci, sortedbamname, fastaname, outputseqprefix, outputalnprefix, num_of_proc):
    '''
    Generate num_of_proc fasta files that contains all the loci/extended region sequence.
    Generate num_of_proc dump files that contains all the loci/extended region alignment infomation.

    The entries in the dump file are the same order as in the output fasta files.
    '''

    seqids = sorted(dict_loci.keys())
    num_loci = 0
    num_seq = 0
    total_loci = 0
    for seqid in seqids:
        total_loci += len(dict_loci[seqid])
    # How many loci in one file?
    num_in_one_file = int(total_loci/num_of_proc) + 1
    cur_loci = 0
    cur_fout = 0
    fname = outputseqprefix+"_0.fa"
    fout  = open(fname, 'w')
    fnamedump = outputalnprefix+"_0.dump"
    foutdump  = open(fnamedump, 'wb')
    ret_list = []  # a list of (name, num_of_seq) tuples

    dump_key = None
    dump_value = []
    for seqid in seqids:
        for loci in dict_loci[seqid]:
            num_loci += 1
            if cur_loci < num_in_one_file * (cur_fout+1):
                pass
            else:
                cur_fout += 1
                fout.close()
                #ret_list.append((fname, num_seq))
                ret_list.append((fname, fnamedump))
                fname = outputseqprefix+"_"+str(cur_fout)+".fa"
                fout = open(fname, 'w')
                fnamedump = outputalnprefix+"_"+str(cur_fout)+".dump"
                foutdump  = open(fnamedump, 'wb')
            cur_loci += 1
            plus_seq = []
            minus_seq = []
            plus_tag = []
            minus_tag = []

            plus_dumpinfo = []
            minus_dumpinfo = []
            loci_pos = str(loci[0][0]) + "-" + str(loci[0][1])
            for idx, extendregion in enumerate(loci[1:]):
                #extendregion is [), but faidx is []
                region = seqid+":"+str(extendregion[0][0])+"-"+str(extendregion[0][1]-1)
                regionpos = seqid+":"+str(extendregion[0][0])+"-"+str(extendregion[0][1])
                command = "samtools faidx "+fastaname+" "+region
                seq = ""
                try:
                    samtools_process = subprocess.Popen(command.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    outmessage, outerr = samtools_process.communicate()
                    seq = str("".join(outmessage.decode().split("\n")[1:]))
                except Exception as e:
                    sys.stderr.write(e)
                    sys.stderr.write("Error occurred when cutting loci sequence.\n")
                    sys.exit(-1)
                strands = set()
                otherinfo = ""
                for peak in extendregion[1]:
                    strands.add(peak[2])
                    otherinfo = otherinfo + str(peak[0])+","+str(peak[1])+","+str(peak[2])+";"
                otherinfo = otherinfo.rstrip(";")
                seqtag = ""
                # a tag indicates which extend region this sequence belongs to.
                # Possible values:
                # "0": This loci has only one extend region, and 'seq' is the sequence for the region
                # "L": This loci has two extend regions, and 'seq' is the sequence for the left region
                # "R": This loci has two extend regions, and 'seq' is the sequence for the right region
                extendregiontag = " 0 "
                if len(loci)==3 and idx == 0:  # Left side region
                    extendregiontag = " L "
                elif len(loci)==3 and idx == 1:  # Right side region
                    extendregiontag = " R "
                if len(strands) == 1:
                    if "+" in strands:
                        seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo+"\n"  # >ExtendRegion strand loci 0/L/R peaks
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                    else:
                        seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                        seq = get_reverse_complement(seq)
                    fout.write(seqtag)
                    fout.write(seq+"\n")
                    num_seq += 1
                    alignments = samtools_view_region(sortedbamname, seqid,
                                                      extendregion[0][0],
                                                      extendregion[0][1])
                    dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                    matures = gen_possible_matures_loci(dict_loci_info, extendregion)
                    cPickle.dump([dump_key, extendregiontag.strip(), dict_loci_info, matures], foutdump, protocol=2)
                    continue
                if len(strands) == 2:
                    seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo+"\n"
                    plus_tag.append(seqtag)
                    plus_seq.append(seq)
                    #fout.write(seqtag)
                    #fout.write(seq+"\n")
                    num_seq += 1
                    seq = get_reverse_complement(seq)
                    seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                    minus_tag.append(seqtag)
                    minus_seq.append(seq)
                    num_seq += 1

                    alignments = samtools_view_region(sortedbamname, seqid,
                                                      extendregion[0][0],
                                                      extendregion[0][1])
                    dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                    matures = gen_possible_matures_loci(dict_loci_info, extendregion)
                    dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                    plus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures])
                    dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                    minus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures])

                    #fout.write(seqtag)
                    #fout.write(seq+"\n")

            if len(plus_tag):
                for which, tag in enumerate(plus_tag):
                    fout.write(tag)
                    fout.write(plus_seq[which]+"\n")
                    cPickle.dump(plus_dumpinfo[which], foutdump, protocol=2)
                for which, tag in enumerate(minus_tag):
                    fout.write(minus_tag[which])
                    fout.write(minus_seq[which]+"\n")
                    cPickle.dump(minus_dumpinfo[which], foutdump, protocol=2)
    ret_list.append((fname, fnamedump))
    fout.close()
    return num_loci, num_seq, ret_list

def dump_loci_seqs_and_alignment_multiprocess(dict_loci, piece_info_list,
                                              sortedbamname, fastaname,
                                              outputseqprefix, outputalnprefix,
                                              logger):
    def dump_piece(dict_loci, piece_info_list,sortedbamname, fastaname,
                   outputfastaname, outputalnname, queue):
        fout = open(outputfastaname,'w')
        foutdump = open(outputalnname,'w')
        dump_key = None
        num_seq = 0
        num_loci = 0
        for seqid in piece_info_list[0]:
            start = 0
            end = len(dict_loci[seqid])
            if piece_info_list[0][0] == seqid:  # the first seqid
                start = piece_info_list[1][0]
            if piece_info_list[0][-1] == seqid:  # the last seqid
                end = piece_info_list[1][1]
            for loci in dict_loci[seqid][start:end]:
                num_loci += 1
                plus_seq = []
                minus_seq = []
                plus_tag = []
                minus_tag = []

                plus_dumpinfo = []
                minus_dumpinfo = []
                loci_pos = str(loci[0][0]) + "-" + str(loci[0][1])
                for idx, extendregion in enumerate(loci[1:]):
                    #extendregion is [), but faidx is []
                    region = seqid+":"+str(extendregion[0][0])+"-"+str(extendregion[0][1]-1)
                    regionpos = seqid+":"+str(extendregion[0][0])+"-"+str(extendregion[0][1])
                    command = "samtools faidx "+fastaname+" "+region
                    seq = ""
                    try:
                        samtools_process = subprocess.Popen(command.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        outmessage, outerr = samtools_process.communicate()
                        seq = str("".join(outmessage.decode().split("\n")[1:]))
                    except Exception as e:
                        sys.stderr.write(e)
                        sys.stderr.write("Error occurred when cutting loci sequence.\n")
                        sys.exit(-1)
                    strands = set()
                    otherinfo = ""
                    for peak in extendregion[1]:
                        strands.add(peak[2])
                        otherinfo = otherinfo + str(peak[0])+","+str(peak[1])+","+str(peak[2])+";"
                    otherinfo = otherinfo.rstrip(";")
                    seqtag = ""
                    extendregiontag = " 0 "
                    if len(loci)==3 and idx == 0:  # Left side region
                        extendregiontag = " L "
                    elif len(loci)==3 and idx == 1:  # Right side region
                        extendregiontag = " R "
                    if len(strands) == 1:
                        if "+" in strands:
                            seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo+"\n"  # >ExtendRegion strand loci 0/L/R peaks
                            dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                        else:
                            seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                            dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                            seq = get_reverse_complement(seq)
                        fout.write(seqtag)
                        fout.write(seq+"\n")
                        num_seq += 1
                        alignments = samtools_view_region(sortedbamname, seqid,
                                                          extendregion[0][0],
                                                          extendregion[0][1])
                        dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                        matures = gen_possible_matures_loci(dict_loci_info, extendregion)
                        cPickle.dump([dump_key, extendregiontag.strip(), dict_loci_info, matures], foutdump, protocol=2)
                        continue
                    if len(strands) == 2:
                        seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo+"\n"
                        plus_tag.append(seqtag)
                        plus_seq.append(seq)
                        num_seq += 1
                        seq = get_reverse_complement(seq)
                        seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo+"\n"
                        minus_tag.append(seqtag)
                        minus_seq.append(seq)
                        num_seq += 1
                        alignments = samtools_view_region(sortedbamname, seqid,
                                                          extendregion[0][0],
                                                          extendregion[0][1])
                        dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                        matures = gen_possible_matures_loci(dict_loci_info, extendregion)
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                        plus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures])
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                        minus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures])
                if len(plus_tag):
                    for which, tag in enumerate(plus_tag):
                        fout.write(tag)
                        fout.write(plus_seq[which]+"\n")
                        cPickle.dump(plus_dumpinfo[which], foutdump, protocol=2)
                    for which, tag in enumerate(minus_tag):
                        fout.write(minus_tag[which])
                        fout.write(minus_seq[which]+"\n")
                        cPickle.dump(minus_dumpinfo[which], foutdump, protocol=2)
        fout.close()
        foutdump.close()
        queue.put(((outputfastaname, outputalnname), num_loci, num_seq))
        queue.close()

    if logger:
        logger.info("Generating candidate loci fasta sequences and loci reads alignment information. "
                    +str(len(piece_info_list))+" parallel processes.")
    inforqueue = multiprocessing.Queue()
    jobs = []
    finalresult = []
    for i in range(len(piece_info_list)):
        fname = outputseqprefix+"_"+str(i)+".fa"
        fnamedump = outputalnprefix+"_"+str(i)+".dump"
        p = multiprocessing.Process(target = dump_piece, args=(dict_loci,
                                                               piece_info_list[i],
                                                               sortedbamname,
                                                               fastaname,
                                                               fname,
                                                               fnamedump,
                                                               inforqueue))
        p.start()
        jobs.append(p)
    total_loci = 0
    total_seq = 0
    for job in jobs:
        job.join()
        info = inforqueue.get()
        total_loci += info[1]
        total_seq += info[2]
        finalresult.append(info[0])
    finalresult.sort()
    return total_loci, total_seq, finalresult

def gen_candidate_region_typeA(dict_contigs, dict_len, dict_option, tmpdir,
                               logger):
    '''
    Connect contigs that have a distance smaller than MAX_GAP.

    dict_loci structure:
    key: seqid
    value: a list of regioninfo. A regioninfo is:
    [region, (extendedregion1, peaks1), (extendedregion2, peaks2)]
    '''
    def next_region_typeA(contiglist, maxgap):
        if len(contiglist) == 0:
            raise StopIteration
        r_now = contiglist[0]
        r_peak = [contiglist[0]]
        for i in xrange(1,len(contiglist)):
            r_next = contiglist[i]
            if r_next[0] - r_now[1] < maxgap:
                r_now = (r_now[0], r_next[1])
                r_peak.append(contiglist[i])
            else:
                yield r_now, r_peak
                r_now = r_next
                r_peak = [contiglist[i]]
        yield r_now, r_peak

    def extend_region(region, max_transcript, seqlen):
        length = region[1] - region[0]
        if length > max_transcript + 50: #if the region is too long, ignore it.
            return []
        if length > max_transcript:
            return [region]
        if length < 60: #shortest precursor is about 60
            leftstart = region[0] - (max_transcript-length-25)-25
            if leftstart<0:
                leftstart=0
            leftend = region[1] + 25
            if leftend > seqlen:
                leftend = seqlen
            rightstart = region[0] - 25
            if rightstart < 0:
                rightstart = 0
            rightend = region[1] + (max_transcript-length-25) + 25
            if rightend > seqlen:
                rightend = seqlen
            return [(leftstart, leftend), (rightstart, rightend)]
        else:
            left = region[0] - 25
            right = region[1] + 25
            if left < 0:
                left = 0
            if right > seqlen:
                right = seqlen
            return [(left, right)]

    dict_loci = {}
    total_loci = 0  # the number of loci
    if logger:
        logger.info("Connect contigs that have a distance smaller than MAX_GAP.")
    for seqID in sorted(dict_contigs): #regions in the dict are already sorted.
        dict_loci[seqID] = []  # changed to a list, this preserves the order.
        for region, peaks in next_region_typeA(dict_contigs[seqID], dict_option["MAX_GAP"]):
            regioninfo = [region]  # [region, (extendedregion1, peaks1), (extendedregion2, peaks2)]
            for extendedregion in extend_region(region, dict_option["PRECURSOR_LEN"], dict_len[seqID]):
                regioninfo.append((extendedregion, peaks))
            if len(regioninfo) > 1:  #  this means the previous for loop was executed at least once.
                dict_loci[seqID].append(regioninfo)
                total_loci += 1
    #  generate information for using multiple processes to dumping loci info in
    #  the next step. For each piece, we record the seq IDs the piece have, and
    #  the start index and end index in the first seqid and the last seqid. Then
    #  the process for the piece will process loci from first_seqid[start] to
    #  last_seqid[end]. The info is recorded in a list. The list has
    #  NUM_OF_CORE values. Each element of the list is an two element (A and B)
    #  list. A is a list of all seq IDs for this piece (ordered). B is a 2-tuple
    #  contains the start and the end index for this piece
    num_each_piece = int(total_loci/dict_option["NUM_OF_CORE"]) + 1
    piece_info = []
    cur_piece = 0
    cur_loci = 0
    cur_pieceinfo = [[], [0,0]]
    cur_start = 0
    cur_end = 0
    cur_id = 0
    for seqID in sorted(dict_loci):
        cur_id = seqID
        cur_pieceinfo[0].append(seqID)
        for idx, loci in enumerate(dict_loci[seqID]):
            cur_loci += 1
            if cur_loci < num_each_piece * (cur_piece + 1):
                pass
            else:
                cur_end = idx
                cur_pieceinfo[1][0] = cur_start
                cur_pieceinfo[1][1] = cur_end
                piece_info.append(cur_pieceinfo)
                cur_start = idx
                cur_pieceinfo = [[seqID], [cur_end,0]]
                cur_piece += 1
    if len(piece_info)<dict_option["NUM_OF_CORE"]:  # last piece
        cur_pieceinfo[1][1] = len(dict_loci[cur_id])
        piece_info.append(cur_pieceinfo)

    ###############debug###############
    tempgff = os.path.join(tmpdir,dict_option["NAME_PREFIX"]+"_ExRegionA.gff3")
    fgff = open(tempgff,'w')
    for seqID in sorted(dict_contigs): #regions in the dict are already sorted.
        cnt = 0
        for region, peaks in next_region_typeA(dict_contigs[seqID], dict_option["MAX_GAP"]):
            otherinfo = ""
            for peak in peaks:
                otherinfo = otherinfo + ":".join([str(x) for x in peak]) + "|"
            for extendedregion in extend_region(region, dict_option["PRECURSOR_LEN"], dict_len[seqID]):
                name = "ExRegionA_"+str(cnt)
                cnt += 1
                write_gff_line(seqID,extendedregion[0],extendedregion[1],"+",name,name,feature="ExRegionA",other=otherinfo,fout=fgff)
    fgff.close()
    ###################################
    return dict_loci, piece_info


def samtools_view_region(sortedbamname, seqid, start, end):
    #TODO samtools view generate reads overlap with the region, not only reads
    #in the region. Should this be a problem here??
    region = seqid + ":" + str(start) + "-" + str(end)
    command = "samtools view " + sortedbamname + " "+region
    lines = []
    ret= ""
    try:
        samtools_process = subprocess.Popen(command.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when viewing a region from bamfile.\n")
        sys.exit(-1)
    alignments = outmessage.split("\n")
    if alignments[0]:
        alignments.pop()
        return alignments
    else:  # no alignment overlap with this region.
        return []


def gen_loci_alignment_info(alignments, seqid, loci):

    '''
    Generate a dict contains all the alignment information of a loci.

    The input 'alignments' are result from samtools_view_region for the loci.
    The input 'loci' is an extended region (A dict_loci[seqid][loci_start] value)


    '''
    # alignments for a position in a loci . The format is:
    # key=strand
    # value=[[length, count, total-count][(length1,count1),(length2, count2),...]]
    # length: the length of the read start from this position that are most
    # abundant. count: the count of that read. total-count: the total number of
    # reads that start from this position. length1, length2,..lengthN: the
    # lengths of reads starts from this position.
    # total-count = count1 + count2 + ... + countN
    # dict_pos_info = {"+":[[0, 0, 0],[]],
    #                   "-":[[0, 0, 0],[]]
    # }
    # key is position, value is a dict_pos_info
    dict_loci_info = {"+": 0, # total reads count of plus
                      "-": 0, # total reads count of plus
                      0: {} # detailed reads info
                  }

    for aln in alignments:
        sp = aln.split()
        startpos = 0
        try:  # skip unmapped alignment
            startpos = int(sp[3])
        except ValueError as e:
            continue
        if startpos < loci[0][0] or startpos >loci[0][1]:
            continue
        depth = int(sp[0].split("-")[-1])  # read names must be in "name-depth" format
        strand = "+"
        if int(sp[1]) & 16:  # on minus strand
            strand = "-"
        cur_seqid = sp[2]
        if cur_seqid != seqid:  # This should not happen
            continue
        readlen = len(sp[9])
        value = (readlen, depth)
        if startpos in dict_loci_info[0]:
            if strand in dict_loci_info[0][startpos]:
                dict_loci_info[0][startpos][strand][1].append(value)
                dict_loci_info[0][startpos][strand][0][2] += depth
                dict_loci_info[strand] += depth
                if depth > dict_loci_info[0][startpos][strand][0][1]:
                    dict_loci_info[0][startpos][strand][0][0] = readlen
                    dict_loci_info[0][startpos][strand][0][1] = depth
            else:  # There is an entry for this position, but not for the strand
                dict_loci_info[0][startpos][strand] = [[readlen, depth, depth],[(readlen, depth)]]
                dict_loci_info[strand] = +depth
        else:  # There is no entry for this position, create one.
            dict_pos_info = {strand: [[readlen, depth, depth],[(readlen, depth)]]
                         }
            dict_loci_info[0][startpos] = dict_pos_info
            dict_loci_info[strand] += depth
    return dict_loci_info


def gen_possible_matures_loci(dict_loci_info, loci):
    ret_matures = []
    for start, end, strand in loci[1]:
        mature_depth = 0
        mature_length = 0
        mature_pos = 0
        for pos in xrange(start,end):
            if pos in dict_loci_info[0]:
                if strand in dict_loci_info[0][pos]:
                    if dict_loci_info[0][pos][strand][0][1] > mature_depth:
                        mature_depth = dict_loci_info[0][pos][strand][0][1]
                        mature_length = dict_loci_info[0][pos][strand][0][0]
                        mature_pos = pos
        mature = (mature_pos, mature_pos+mature_length, strand, mature_depth)
        ret_matures.append(mature)
    return ret_matures



########################################################################
##  test mature region
########################################################################
def test_mature_region(dict_loci, dict_option, sortedbamname, outputname):
    outf = open(outputname, 'w')
    count = 0
    for seqid in dict_loci:
        for startpos in dict_loci[seqid]:
            loci_start = dict_loci[seqid][startpos][0][0]
            loci_end = dict_loci[seqid][startpos][0][1]
            alignments = samtools_view_region(sortedbamname, seqid, loci_start, loci_end)
            dict_loci_info = gen_loci_alignment_info(alignments, seqid, dict_loci[seqid][startpos])
            matures = gen_possible_matures_loci(dict_loci_info, dict_loci[seqid][startpos])
            for start, end, strand, depth in matures:
                name = "mature_"+str(count)
                count += 1
                otherinfo = "depth="+str(depth)
                write_gff_line(seqid, start, end, strand, name, name, feature="Mature", other=otherinfo, fout=outf)
            gc.collect()
    outf.close()
    return outputname


########################################################################
#  Structure related functions
########################################################################

def get_structures_next_extendregion(rnalfoldoutname, minlen, filter=False,
                                     minloop=3):
    '''
    Parse the RNALfold output file, return a list of (0/L/R, structure,
    norm_energy) tuples for an extend region each time (generator).

    Structure whose length is shorter than minlen are filtered.

    If filter=True, then structures which are not stem loop, or the loop size is
    smaller than minloop, are filtered.

    '''

    energypattern = re.compile(r'\(\s*(-?[0-9]+[\.]?[0-9]*)\s*\)')
    structures = []
    which = 0
    first = True
    sp=[]
    with open(rnalfoldoutname) as f:
        for line in f:
            sp = line.strip().split()
            if line.startswith(">"):
                if not first:
                    yield (which, structures)
                structures = []
                which = sp[3]  # 0, L, or R
                first = False
            else:
                if len(sp) >=3:  # structure line
                    if len(sp[0]) < minlen:
                        continue
                    elif (not is_stem_loop(sp[0],minloop)) & filter:
                        continue
                    norm_energy = energypattern.search(line).group(1)
                    norm_energy = float(norm_energy)/len(sp[0])
                    structures.append((norm_energy,int(sp[-1]), sp[0]))
                else:
                    continue
        yield (which, structures)


def is_stem_loop(ss, minloopsize):
    last_open = ss.rfind("(")
    first_close = ss.find(")")
    if first_close - last_open-1 >= minloopsize:
        return True
    else:
        return False

def stemloop_or_onebifurcation(ss):
    if is_stem_loop(ss, 3):
        return True
    # not an stem loop, but a structure with two bifurcation: (()()).
    bifurcation = 0
    stack = []
    last = 'Push'
    for idx, char in enumerate(ss):
        if char == '(':
            if idx !=0:
                if not stack:
                    if last == 'Pop': # ()() like structure
                        return False
                    last = "Push"
                    stack.append(char)
                else: # stack is not empty
                    if last == "Pop":
                        if bifurcation >= 1:
                            return False
                        else:
                            bifurcation = 1
                    stack.append(char)
                    last = "Push"
            else:
                stack.append(char)
                last = "Push"
        elif char == ")":
            last = "Pop"
            stack.pop()
    return True

def pos_genome_2_local(genomestart, genomeend, strand, regionstart, regionend,
                       foldstart, foldend):
    '''
    Convert genome coordinate to local structure coordinate.

    Local coordinate is 0 based, starts from foldstart.
    NOTE: All the regions are [start, end). That is, the start is inclusive and
    the end is exclusive.
     -------------------------------------------------------------------> genome
          |>----------------->------------------->------------>|   extend region
      regionstart                                         regionend (genome coordinate, start at 1)
                |--------------------------------->|   folded region
               foldstart                         foldend (regionstart is 1)
                      |------------>|
                 localstart  localend  (local coordinate, foldstart is 0)
                 genomestart genomeend  (genome coordinate, start at 1)

    If strand is "+":
    genomestart = regionstart + foldstart - 1 + localstart
    genomeend =  regionstart + foldstart -1 + localend

    If strand is "-":

     -------------------------------------------------------------------> genome
          |<-----------------<-------------------<------------<|   extend region
      regionstart                                         regionend (genome coordinate)
                |<---------------------------------<|   folded region
               foldend                         foldstart
                      |<------------|
                 localend  localstart  (local coordinate, foldstart is 0)
                 genomestart genomeend  (genome coordinate)

    genomestart =  regionend-1 - foldstart+1 -localend+1
    genomeend = regionend-1 - foldstart+1 - localstart+1

    Note that genomestart corresponds to regionend, genomeend corresponds to
    regionend when strand is "-".
    '''
    if strand == "+":
        return (genomestart-regionstart-foldstart+1, genomeend-regionstart-foldstart+1)
    else:
        return (regionend-genomeend-foldstart+1, regionend-genomestart-foldstart+1)


def pos_local_2_genome(localstart, localend, strand, regionstart, regionend,
                       foldstart, foldend):
    '''
    Convert local structure coordinate to genome coordinate.

    NOTE: All the regions are [start, end). That is, the start is inclusive and
    the end is exclusive.

     -------------------------------------------------------------------> genome
          |>----------------->------------------->------------>|   extend region
      regionstart                                         regionend (genome coordinate, start at 1)
                |--------------------------------->|   folded region
               foldstart                         foldend (regionstart is 1)
                      |------------>|
                 localstart  localend  (local coordinate, foldstart is 0)
                 genomestart genomeend  (genome coordinate, start at 1)

    If strand is "+":
    genomestart = regionstart + foldstart - 1 + localstart
    genomeend =  regionstart + foldstart -1 + localend


    If strand is "-":
     -------------------------------------------------------------------> genome
          |<-----------------<-------------------<------------<|   extend region
      regionstart                                         regionend (genome coordinate)
                |<---------------------------------<|   folded region
               foldend                         foldstart (regionend-1 is 1)
                      |<------------|
                 localend  localstart  (local coordinate, foldstart is 0)
                 genomestart genomeend  (genome coordinate)

    genomestart =  regionend-1 - foldstart+1 -localend+1
    genomeend = regionend-1 - foldstart+1 - localstart+1

    Note that genomestart corresponds to regionend, genomeend corresponds to
    regionend when strand is "-".
    '''
    if strand == "+":
        return (regionstart+foldstart-1+localstart, regionstart+foldstart-1+localend)
    else:
        return (regionend-foldstart-localend+1, regionend-foldstart-localstart+1)


def get_maturestar_info(ss, mature, foldstart, foldend, regionstart, regionend,
                        strand):

    '''
    mature format is [maturestart, matureend) in GENOME coordinate.

    Make sure the mature is in one arm of a stem.
    '''

    dict_bp = {}  #
    stack = []
    for idx, char in enumerate(ss):
        if char =="(":
            stack.append(idx)
        if char == ")":
            if len(stack) == 0:
                return None
            dict_bp[idx] = stack[-1]
            dict_bp[stack[-1]] = idx
            stack.pop()

    # zero started
    mature_local_pos = pos_genome_2_local(mature[0], mature[1], strand,
                                          regionstart, regionend, foldstart,
                                          foldend)
    foldregion_genome = pos_local_2_genome(0, len(ss), strand,
                                          regionstart, regionend, foldstart,
                                          foldend)
    # mature not in the fold region
    if not (mature[0]>=foldregion_genome[0] and mature[1]<=foldregion_genome[1]):
        return None

    mature_ss = ss[mature_local_pos[0]:mature_local_pos[1]]
    #  Not in a arm of a stem, return None. This should not happen given
    #  that the mature is already in an arm of a stem.
    if mature_ss.find("(")!=-1 and mature_ss.find(")")!=-1:
        return None

    #  local positions
    star_start = 0
    star_end = 0
    star_ss = ""

    firstbp = ss.find("(", mature_local_pos[0], mature_local_pos[1])
    lastbp =  ss.rfind("(", mature_local_pos[0], mature_local_pos[1])
    prime5 = True
    if firstbp == -1:
        firstbp = ss.find(")",mature_local_pos[0], mature_local_pos[1])
        lastbp =  ss.rfind(")",mature_local_pos[0], mature_local_pos[1])
        prime5 = False
    if firstbp == -1:  # all are dots
        return None

    start_unpaired_len = mature_local_pos[1]-1-lastbp
    end_unpaired_len = firstbp-mature_local_pos[0]
    star_start = dict_bp[lastbp] -start_unpaired_len + 2
    star_end = dict_bp[firstbp] + end_unpaired_len + 2 + 1  # +1 because the end pos is exclusive

    #duplex is the region from first paired to the last paired
    mature_duplex = ss[firstbp: lastbp]
    star_duplex = ss[dict_bp[lastbp]: dict_bp[firstbp]]
    total_dots = mature_duplex.count(".") + star_duplex.count(".")
    total_bps = len(mature_duplex) - mature_duplex.count(".")
    if total_bps<14:
        return None
    # if prime5:  # the mature is on 5' arm
    #     print("5PRIME")
    #     start_unpaired_len = mature_local_pos[1]-1-lastbp
    #     end_unpaired_len = firstbp-mature_local_pos[0]
    #     star_start = dict_bp[lastbp] -start_unpaired_len + 2
    #     star_end = dict_bp[first] + end_unpaired_len + 2 + 1  # +1 because the end pos is exclusive
    #     star_ss = ss[star_start:star_end]
    # else:  # the mature is on 3' arm
    #     print("3PRIME")
    #     start_unpaired_len = mature_local_pos[1]-1-lastbp
    #     end_unpaired_len = firstbp-mature_local_pos[0]
    #     star_start = dict_bp[lastbp] - start_unpaired_len + 2
    #     star_end = dict_bp[firstbp] + end_unpaired_len + 2 + 1
    star_ss = ss[star_start: star_end]
    if star_ss.find("(") != -1 and star_ss.find(")")!= -1:  # this could happen if the input ss in not a stem loop
        return None
    # convert local positions to genome positions.
    genome_star_start, genome_star_end = pos_local_2_genome(star_start,
                                                             star_end,
                                                             strand,
                                                             regionstart,
                                                             regionend,
                                                             foldstart,
                                                             foldend)
    #  return:
    #  star start, star end in genome coordinate
    #  foldregion start, end in genome coordinate
    #  star structure
    #  whether the mature is at the 5' arm of the folding sequence
    #  mature structure
    #  total number of unmatched bases in the duplex
    #  total number of basepairs
    return genome_star_start, genome_star_end, foldregion_genome[0], foldregion_genome[1], star_ss, prime5, mature_ss, total_dots, total_bps

def check_expression(ss, dict_extendregion_info, maturepos_genome, mature_depth, starpos_genome, strand):
    '''
    starpos_genome is the genome coordinate of the star sequence position.

    The format of starpos_genome is [start, end)
    The return value is the depth of the star sequence.
    '''
    starlen = starpos_genome[1] - starpos_genome[0]
    total_depth = 0
    star_depth = 0
    for pos in dict_extendregion_info[0]:
        for s in dict_extendregion_info[0][pos]:
            total_depth += dict_extendregion_info[0][pos][s][0][-1]
    if starpos_genome[0] not in dict_extendregion_info[0]:
        star_depth = 0
    elif strand not in dict_extendregion_info[0][starpos_genome[0]]:
        star_depth = 0

    #TODO  Maybe better to store length as key depth as value, easier for lookup.
    else:
        for length, depth in  dict_extendregion_info[0][starpos_genome[0]][strand][1]:
            if length == starlen:
                star_depth = depth

    ratio = (mature_depth+star_depth)/float(total_depth)
    #print("T, M, S, R", total_depth, mature_depth, star_depth, ratio)
    return star_depth, ratio, total_depth


def check_loci(structures, matures, region, dict_aln, which):
    miRNAs = []
    for m0, m1, strand, mdepth in matures:
        # TODO Move this if to previous stage.
        # if mature length is not in [19-24], then ignore this
        if m1-m0 <19 or m1-m0>24:
            continue
        lowest_energy = 0
        lowest_energy = 0
        outputinfo = []
        for energy, foldstart, ss in structures[1]:
            if not stemloop_or_onebifurcation(ss):
                continue
            if energy > lowest_energy:
                continue
            else:
                structinfo = get_maturestar_info(ss, (m0,m1), foldstart,
                                                 foldstart+len(ss),
                                                 region[1][0],
                                                 region[1][1],
                                                 strand)
                if structinfo:
                    exprinfo = check_expression(ss,dict_aln, (m0,m1),
                                                mdepth,
                                                (structinfo[0],structinfo[1]),
                                                strand)
                    if exprinfo[0]>0:  # has star expression
                        if exprinfo[1] < 0.10:
                            continue
                        else:
                            #  The last 'True' means this is an confident miRNA
                            outputinfo = [region[0], structinfo[2],
                                          structinfo[3], m0, m1, structinfo[0],
                                          structinfo[1], ss, strand, True]
                            lowest_energy = energy
                    else:  #  no star expression
                        if exprinfo[1] >=0.8 and exprinfo[2] > 5000:  # but very high expression
                            #  The last 'False' means this is not an confident miRNA
                            outputinfo = [region[0], structinfo[2],
                                          structinfo[3], m0, m1, structinfo[0],
                                          structinfo[1], ss, strand, False]
                            lowest_energy = energy
        if outputinfo:  # the loci contains an miRNA
            miRNAs.append(outputinfo)
    return miRNAs


def filter_next_loci(alndumpname, rnalfoldoutname, minlen=50, filterstemloop=False):
    list_miRNA_loci = []
    alnf = open(alndumpname)
    ss_generator = get_structures_next_extendregion(rnalfoldoutname, minlen,
                                                    filter=filterstemloop)
    while True:
        region = None
        which = None
        dict_aln = None
        matures = None
        try:
            region, which, dict_aln, matures = cPickle.load(alnf)
        except EOFError:
            raise StopIteration
        structures = []
        if which == "0": #only one extend region for this loci
            structures = next(ss_generator)
            miRNAs = check_loci(structures, matures, region, dict_aln, which)
            if miRNAs:
                yield miRNAs
            continue
        else:
            structuresL = next(ss_generator)  # L
            structuresR = next(ss_generator)  # R
            region1, which1, dict_aln1, matures1 = cPickle.load(alnf)
            miRNAsL = check_loci(structuresL, matures, region, dict_aln, which)
            if miRNAsL:
                yield miRNAsL
            else:
                miRNAsR = check_loci(structuresR, matures1, region1, dict_aln1, which1)
                if miRNAsR:
                    yield miRNAsR
                continue

def gen_miRNA_loci_nopredict(alndumpnames, rnalfoldoutnames, minlen, logger,
                             filterstemloop=False):
    def gen_miRNA_loci_local(queue, alndumpname, rnalfoldoutname, minlen,
                             filterstemloop=False):
        mir_generator = filter_next_loci( alndumpname, rnalfoldoutname,
                                          minlen=minlen,
                                          filterstemloop=filterstemloop)
        mirnas = []
        for mir in mir_generator:
            mirnas.append(mir)
        queue.put(mirnas)
        queue.put("done")
        queue.close()
    miRNAqueue = multiprocessing.Queue()
    jobs = []
    finalresult = []
    for i in range(len(alndumpnames)):
        p = multiprocessing.Process(target = gen_miRNA_loci_local, args=(miRNAqueue, alndumpnames[i], rnalfoldoutnames[i],minlen))
        p.start()
        jobs.append(p)

    num_joined = 0
    while True:
        try:
            if num_joined == len(jobs):
                for job in jobs:
                    job.join()
                break
            result = miRNAqueue.get_nowait()
            if result == "done":
                num_joined += 1
            else:
                for mirna in result:
                    finalresult.append(mirna[0])
        except Queue.Empty:
            time.sleep(5)
            continue

    return finalresult


def get_mature_stemloop_seq(seqid, mstart, mend, start, end, fastaname):
    region1 = seqid+":"+str(mstart)+"-"+str(mend)
    region2 = seqid+":"+str(start)+"-"+str(end)
    mature_command = "samtools faidx " + fastaname + " " +region1
    stemloop_command = "samtools faidx " + fastaname + " " +region2
    try:
        samtools_process = subprocess.Popen(mature_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
        mature_seq = str("".join(outmessage.decode().split("\n")[1:]))
        samtools_process = subprocess.Popen(stemloop_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
        stemloop_seq = str("".join(outmessage.decode().split("\n")[1:]))
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when indexing the genome file\n")
        sys.exit(-1)
    return region1, mature_seq, region2, stemloop_seq

def gen_gff_from_result(resultlist, gffname):
    f = open(gffname, 'w')
    resultlist.sort()
    for idx, mirna in enumerate(resultlist):
        maturename = "miRNA_"+ str(idx)
        mirname = "miRNA-stemloop_"+ str(idx)
        write_gff_line(mirna[0],mirna[1],mirna[2]-1,mirna[-2],mirname, mirname,
                       feature="miRNA-stemloop",fout=f)
        write_gff_line(mirna[0],mirna[3],mirna[4]-1,mirna[-2],maturename,
                       maturename, feature="miRNA",fout=f)
    f.close()

def gen_mirna_fasta_from_result(resultlist, maturename, stemloopname, fastaname):
    f1 = open(maturename, 'w')
    f2 = open(stemloopname, 'w')
    for idx, mirna in enumerate(resultlist):
        maturename = "miRNA_" + str(idx)
        mirname = "miRNA-stemloop_" + str(idx)
        seqs = get_mature_stemloop_seq(mirna[0],mirna[3], mirna[4]-1, mirna[1], mirna[2]-1, fastaname)
        matureid = ">" + seqs[0] + " " + mirna[-2]+ " " + mirname
        stemloopid = ">" + seqs[2] + " " + mirna[-2]+ " " + mirname
        matureseq = seqs[1].upper().replace("T", "U")
        stemloopseq = seqs[3].upper().replace("T", "U")
        if mirna[-2] == "-":
            matureseq = get_reverse_complement(matureseq)
            stemloopseq = get_reverse_complement(stemloopseq)
        f1.write(matureid+"\n")
        f1.write(matureseq+"\n")
        f2.write(stemloopid+"\n")
        f2.write(stemloopseq+"\n")
    f1.close()
    f2.close()

def fold_use_RNALfold(inputfastalist, tempfolder, dict_option, maxspan):
    '''
    Fold the sequences using RNALfold and write results to outputname file.
    '''
    inputfastalist.sort()
    def fold(inputname,outputname):
        command = "RNALfold -L " + str(maxspan)
        try:
            outputfile = open(outputname, 'w')
            inputfile = open(inputname)
            subprocess.check_call(command.split(), stdin=inputfile, stdout=outputfile)
            outputfile.close()
        except Exception as e:
            sys.stderr.write(e)
            sys.stderr.write("Error occurred when folding sequences.\n")
            sys.exit(-1)
    outputnames = []
    for i in range(len(inputfastalist)):
        outname = os.path.join(tempfolder, dict_option["NAME_PREFIX"]+"_rnalfoldoutput_"+str(i))
        outputnames.append(outname)
    jobs = []
    for i in range(len(inputfastalist)):
        p = multiprocessing.Process(target = fold, args=(inputfastalist[i], outputnames[i]))
        p.start()
        jobs.append(p)
    for job in jobs:
        job.join()
    return outputnames


def files_all_exist(filelist):
    allexist = True
    for name in filelist:
        if not os.path.exists(name):
            allexist = False
    return allexist

def previous_stage_saved(recovername, stage):
    if not os.path.exists(recovername):
        return False
    recoverfile = open(recovername)
    dict_recover = cPickle.load(recoverfile)
    if stage in dict_recover["finished_stages"]:
        if files_all_exist(dict_recover["files"][stage]):
            return True
    return False

def detect_stage_last_finished(recovername):
    if not os.path.exists(recovername):
        return None
    recoverfile = open(recovername)
    dict_recover = cPickle.load(recoverfile)
    return dict_recover["last_stage"]

def load_recover_file(recovername):
    if not os.path.exists(recovername):
        return None
    recoverfile = open(recovername)
    dict_recover = cPickle.load(recoverfile)
    return dict_recover

def run_check(dict_option, outtempfolder, recovername):
    write_formatted_string_withtime("Checking RNALfold and samtools", 30, sys.stdout)
    if not check_RNALfold():
        write_formatted_string("*** RNALfold is required but not installed or not in the PATH.", 30, sys.stdout)
        write_formatted_string("*** Please refer to the README file.", 30, sys.stdout)
    else:
        RNALfoldversion = get_RNALfold_version()
        if is_bug_RNALfold(RNALfoldversion.strip()):
            write_formatted_string("!!! The version of RNALfold (2.0.4) has a bug. ", 30, sys.stdout)
            write_formatted_string("!!! Please use the latest version or the older version 1.8.x", 30, sys.stdout)
            write_formatted_string("!!! Please refer to the README file.", 30, sys.stdout)
        else:
            write_formatted_string("*** RNALfold is ready.", 30, sys.stdout)

    if not check_samtools():
        write_formatted_string("!!! samtools is required but not installed or not in the PATH.", 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file.", 30, sys.stdout)
    else:
        write_formatted_string("*** SAMtools is ready.\n", 30, sys.stdout)

    last_stage = detect_stage_last_finished(recovername)
    if last_stage == "fold":
        write_formatted_string_withtime("Checkpoint info:", 30, sys.stdout)
        write_formatted_string("*** The pipeline was stopped after stage '" + last_stage + "'.", 30, sys.stdout)
        write_formatted_string("*** No need to run any stage (if all data files are not changed). ", 30, sys.stdout)
    elif last_stage:
        write_formatted_string_withtime("Checkpoint info:", 30, sys.stdout)
        write_formatted_string("*** The pipeline was stopped after stage '" + last_stage + "'.", 30, sys.stdout)
        write_formatted_string("*** You can start the pipeline from next stage (if all data files are not changed and all files generated from previous stages exist). ", 30, sys.stdout)
    else:
        write_formatted_string_withtime("Checkpoint info:", 30, sys.stdout)
        write_formatted_string("*** No checkpoint information found, need to run the whole pipeline.", 30, sys.stdout)

def check_dependency():
    write_formatted_string_withtime("Checking RNALfold and samtools", 30, sys.stdout)
    allgood = True
    if not check_RNALfold():
        write_formatted_string("!!! RNALfold is required but not installed or not in the PATH.", 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file.", 30, sys.stdout)
        allgood = False
    else:
        RNALfoldversion = get_RNALfold_version()
        if is_bug_RNALfold(RNALfoldversion):
            write_formatted_string("!!! The version of RNALfold (2.0.4) has a bug. ", 30, sys.stdout)
            write_formatted_string("!!! Please use the latest version or the older version 1.8", 30, sys.stdout)
            write_formatted_string("!!! Please refer to the README file.", 30, sys.stdout)
            allgood = False
        else:
            write_formatted_string("*** RNALfold is ready.", 30, sys.stdout)

    if not check_samtools():
        write_formatted_string("!!! samtools is required but not installed or not in the PATH.", 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file.", 30, sys.stdout)
        allgood = False
    else:
        write_formatted_string("*** SAMtools is ready.\n", 30, sys.stdout)

    return allgood

def run_prepare(dict_option, outtempfolder, recovername):
    logger = None
    if dict_option["LOG"]:
        logger = dict_option["LOGGER"]
        logger.info("Starting preparing data for the 'candidate' stage.")
    write_formatted_string_withtime("Starting preparing data for the 'candidate' stage.\n", 30, sys.stdout)
    combinedbamname, expandedsamname, expandedbamname, expandedbam_plus, expandedbam_minus = prepare_data(dict_option, outtempfolder, logger)
    if logger:
        logger.info("Start writing recover infomation to the recovery file")
    dict_recover["last_stage"] = "prepare"
    dict_recover["finished_stages"][ "prepare"]={}
    dict_recover["finished_stages"]["prepare"]["combinedbamname"] = combinedbamname
    dict_recover["finished_stages"]["prepare"]["expandedbamname"] = expandedbamname
    dict_recover["finished_stages"]["prepare"]["expandedsamname"] = expandedsamname
    dict_recover["finished_stages"]["prepare"]["expandedbam_plus"] = expandedbam_plus
    dict_recover["finished_stages"]["prepare"]["expandedbam_minus"] = expandedbam_minus
    dict_recover["files"]["prepare"] = []
    dict_recover["files"]["prepare"].append(combinedbamname)
    dict_recover["files"]["prepare"].append(expandedsamname)
    dict_recover["files"]["prepare"].append(expandedbamname)
    dict_recover["files"]["prepare"].append(expandedbam_plus)
    dict_recover["files"]["prepare"].append(expandedbam_minus)
    recoverfile = open(recovername,'w')
    cPickle.dump(dict_recover,recoverfile)
    if logger:
        logger.info("Recovery file successfully updated.")
    write_formatted_string_withtime("Done\n", 30, sys.stdout)
    if logger:
        logger.info("=========================Done (prepare stage)=======================\n\n")


def run_candidate(dict_option, outtempfolder, recovername):
    logger = None
    if dict_option["LOG"]:
        logger = dict_option["LOGGER"]
        logger.info("Starting candidate stage.")
    if not previous_stage_saved(recovername, "prepare"):
        write_formatted_string_withtime("Error: can not start the pipeline from this stage, the files needed are not generated or have been removed/moved. Please run previous stages first, or run the pipeline in the recover mode to automatically continue from where the job was ceased.", 30, sys.stdout)
        exit(-1)
    #now we know we can run this stage. First get file names from the last stage
    dict_recover = load_recover_file(recovername)
    if logger:
        logger.info("Starting identifying candidate regions.")
    write_formatted_string_withtime("Starting identifying candidate regions", 30, sys.stdout)
    rnalfoldinputprefix = os.path.join(outtempfolder, dict_option["NAME_PREFIX"]+".rnalfold.in")
    outputdumpprefix = os.path.join(outtempfolder, dict_option["NAME_PREFIX"]+".alndump")
    loci_dump_name = os.path.join(outtempfolder, dict_option["NAME_PREFIX"]+"_loci_dump.dump")
    dict_len = get_length_from_sam(dict_option["ALIGNMENT_FILE"][0])

    combinedbamname = dict_recover["finished_stages"]["prepare"]["combinedbamname"]
    expandedbamname = dict_recover["finished_stages"]["prepare"]["expandedbamname"]
    expandedsamname = dict_recover["finished_stages"]["prepare"]["expandedsamname"]
    expandedbam_plus = dict_recover["finished_stages"]["prepare"]["expandedbam_plus"]
    expandedbam_minus = dict_recover["finished_stages"]["prepare"]["expandedbam_minus"]
    #if the length of a contig is smaller than min_contig_len, then ignore it
    min_contig_len = 19
    depthfilename, dict_contigs = gen_contig_typeA(expandedbam_plus,expandedbam_minus, dict_option, min_contig_len, logger)
    dict_loci, piece_info = gen_candidate_region_typeA(dict_contigs,dict_len,dict_option,outtempfolder, logger)
    #  a list of (rnafoldinfasta, lociinfodump)  pairs
    num_loci, num_fasta, retnames = dump_loci_seqs_and_alignment_multiprocess(dict_loci,
                                                                              piece_info,
                                                                              combinedbamname,
                                                                              dict_option["FASTA_FILE"],
                                                                              rnalfoldinputprefix,
                                                                              outputdumpprefix,
                                                                              logger)

    loci_dump_file = open(loci_dump_name, 'wb')
    cPickle.dump(dict_loci, loci_dump_file, protocol=2)
    loci_dump_file.close()
    if logger:
        logger.info("Start writing recover infomation to the recovery file")
    dict_recover["last_stage"] = "candidate"
    dict_recover["finished_stages"]["candidate"] = {}
    dict_recover["files"]["candidate"] = []
    dict_recover["finished_stages"]["candidate"]["depthfilename"] = depthfilename
    dict_recover["finished_stages"]["candidate"]["loci_dump_name"] = loci_dump_name
    dict_recover["finished_stages"]["candidate"]["fasta"] = []
    dict_recover["finished_stages"]["candidate"]["infodump"] = []
    dict_recover["finished_stages"]["candidate"]["num_loci"] = num_loci
    dict_recover["finished_stages"]["candidate"]["num_fasta"] = num_fasta
    for fastaname, lociinfoname in retnames:
        dict_recover["finished_stages"]["candidate"]["fasta"].append(fastaname)
        dict_recover["finished_stages"]["candidate"]["infodump"].append(lociinfoname)
        dict_recover["files"]["candidate"].append(fastaname)
        dict_recover["files"]["candidate"].append(lociinfoname)

    recoverfile = open(recovername,'w')
    cPickle.dump(dict_recover,recoverfile)
    if logger:
        logger.info("Recovery file successfully updated.")

    outstr = "{0} candidate loci generated, {1} regions to fold.".format(num_loci, num_fasta)
    write_formatted_string(outstr, 30, sys.stdout)
    write_formatted_string_withtime("Done\n", 30, sys.stdout)
    if logger:
        logger.info(outstr)
    if logger:
        logger.info("=========================Done (candidate stage)=======================\n\n")


def run_fold(dict_option, outtempfolder, recovername):
    logger = None
    if dict_option["LOG"]:
        logger = dict_option["LOGGER"]
        logger.info("Starting folding candidate sequences.")
    if not previous_stage_saved(recovername, "candidate"):
        write_formatted_string_withtime("Error: can not start the pipeline from this stage, the files needed are not generated or have been removed/moved. Please run previous stages first, or run the pipeline in the recover mode to automatically continue from where the job was ceased.", 30, sys.stdout)
        exit(-1)
    write_formatted_string_withtime("Starting folding candidate sequences.", 30, sys.stdout)
    #now we know we can run this stage. First get file names from the last stage
    dict_recover = load_recover_file(recovername)
    if logger:
        logger.info("Folding candidate sequences using RNALfold. "+
                    str(dict_option["NUM_OF_CORE"]) + " parallel processes.")
    foldnames = fold_use_RNALfold(dict_recover["finished_stages"]["candidate"]["fasta"],
                                  outtempfolder,
                                  dict_option,
                                  dict_option["PRECURSOR_LEN"])

    if logger:
        logger.info("Start writing recover infomation to the recovery file")
    dict_recover["last_stage"] = "fold"
    dict_recover["finished_stages"]["fold"] = {}
    dict_recover["files"]["fold"] = []
    dict_recover["finished_stages"]["fold"]["foldnames"] = foldnames
    for name in foldnames:
        dict_recover["files"]["fold"].append(name)
    recoverfile = open(recovername,'w')
    cPickle.dump(dict_recover,recoverfile)
    if logger:
        logger.info("Recovery file successfully updated.")
    write_formatted_string_withtime("Done\n", 30, sys.stdout)
    if logger:
        logger.info("=========================Done (fold stage)=======================\n\n")


def run_predict(dict_option, outtempfolder, recovername):
    logger = None
    if dict_option["LOG"]:
        logger = dict_option["LOGGER"]
        logger.info("Starting predicting miRNAs.")
    if not previous_stage_saved(recovername, "fold"):
        write_formatted_string_withtime("Error: can not start the pipeline from this stage, the files needed are not generated or have been removed/moved. Please run previous stages first, or run the pipeline in the recover mode to automatically continue from where the job was ceased.", 30, sys.stdout)
        exit(-1)
    write_formatted_string_withtime("Starting predicting miRNAs.", 30, sys.stdout)
    #now we know we can run this stage. First get file names from the last stage
    dict_recover = load_recover_file(recovername)
    if logger:
        logger.info("Predicting microRNAs. "+
                    str(dict_option["NUM_OF_CORE"]) + " parallel processes.")
    foldnames = dict_recover["finished_stages"]["fold"]["foldnames"]
    result = gen_miRNA_loci_nopredict(dict_recover["finished_stages"]["candidate"]["infodump"],
                                      foldnames, 60, logger)
    gffname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.gff3")
    maturename = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.mature.fa")
    stemloopname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.stemloop.fa")
    if logger:
        logger.info("Generating a gff file contains all predictions. GFF file name: " + gffname)
    gen_gff_from_result(result,gffname)
    if logger:
        logger.info("Generating two fasta files contains predicted mature sequences and stem-loop sequences. Fasta file names: [mature]: " +
                    maturename + ", [stem-loop]: " + stemloopname)
    gen_mirna_fasta_from_result(result, maturename, stemloopname, dict_option["FASTA_FILE"])
    if logger:
        logger.info("Start writing recover infomation to the recovery file")
    dict_recover["last_stage"] = "predict"
    dict_recover["finished_stages"]["predict"] = {}
    dict_recover["files"]["predict"] = []
    dict_recover["finished_stages"]["predict"]["gffname"] = gffname
    dict_recover["finished_stages"]["predict"]["maturename"] = maturename
    dict_recover["finished_stages"]["predict"]["stemloopname"] = stemloopname
    dict_recover["files"]["predict"].append(gffname)
    dict_recover["files"]["predict"].append(maturename)
    dict_recover["files"]["predict"].append(stemloopname)
    recoverfile = open(recovername,'w')
    cPickle.dump(dict_recover,recoverfile)
    if logger:
        logger.info("Recovery file successfully updated.")
    outstr = "{0} miRNAs identified. Result writes to file {1}".format(len(result), gffname)
    write_formatted_string(outstr, 30, sys.stdout)
    write_formatted_string_withtime("Done\n", 30, sys.stdout)
    if logger:
        logger.info(outstr)
    if logger:
        logger.info("=========================Done (predict stage)=======================\n\n")


def run_removetmp(outtempfolder):
    try:
        write_formatted_string_withtime("Removing the temporary folder.", 30, sys.stdout)
        shutil.rmtree(outtempfolder)
        write_formatted_string("Temporary folder removed.\n", 30, sys.stdout)
        return 0
    except Exception as e:
        sys.stderr.write(e)
        sys.stderr.write("Error occurred when removing the tmp folder.\n")
        sys.stderr.write("You can remove the folder by yourself.\n")
        return 1


if __name__ == '__main__':
    dict_option = parse_option_optparse()

    has_multiple_sample = False
    if len(dict_option['ALIGNMENT_FILE'])>1:
        has_multiple_sample = True

    dict_recover = {"finished_stages": {},
                    "last_stage": "",
                    "files": {}  # value: (stage, filelist)
                }

    #make temp dir to store temporary files
    outfolder = dict_option["OUTFOLDER"]
    outtempfolder = os.path.join(outfolder, dict_option["NAME_PREFIX"]+"_tmp")
    recovername = os.path.join(outtempfolder, dict_option["NAME_PREFIX"]+"_recover")

    if not os.path.exists(outfolder):
        os.mkdir(dict_option["OUTFOLDER"])
    if not os.path.exists(outtempfolder):
        os.mkdir(outtempfolder)

    if dict_option['ACTION'] == 'check':
        run_check(dict_option, outtempfolder, recovername)
        exit(0)

    allgood = check_dependency()
    if not allgood:
        exit(-1)

    logger = None
    if dict_option['LOG']:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        logname = os.path.join(dict_option["OUTFOLDER"],
                               dict_option["NAME_PREFIX"]+".log")
        handler = logging.FileHandler(logname, mode='w')
        fmt = logging.Formatter('%(asctime)20s  %(name)10s:  %(levelname)10s  %(message)s')
        handler.setFormatter(fmt)
        logger.addHandler(handler)
        dict_option["LOGGER"] = logger
    if dict_option['ACTION'] == 'prepare':
        if logger:
            logger.info("ACTION is 'prepare'.")
        run_prepare(dict_option, outtempfolder, recovername)
        exit(0)
    if dict_option['ACTION'] == 'candidate':
        if logger:
            logger.info("ACTION is 'candidate'.")
        run_candidate(dict_option, outtempfolder, recovername)
        exit(0)
    if dict_option['ACTION'] == 'fold':
        if logger:
            logger.info("ACTION is 'fold'.")
        run_fold(dict_option, outtempfolder, recovername)
        exit(0)
    if dict_option['ACTION'] == 'predict':
        if logger:
            logger.info("ACTION is 'predict'.")
        run_predict(dict_option, outtempfolder, recovername)
        if not dict_option["KEEPTMP"]:
            run_removetmp(outtempfolder)
            if logger:
                logger.info("Temporary folder removed.")
        exit(0)
    if dict_option['ACTION'] == 'pipeline':
        if logger:
            logger.info("ACTION is 'pipeline'.")
        run_prepare(dict_option, outtempfolder, recovername)
        run_candidate(dict_option, outtempfolder, recovername)
        run_fold(dict_option, outtempfolder, recovername)
        run_predict(dict_option, outtempfolder, recovername)
        if not dict_option["KEEPTMP"]:
            run_removetmp(outtempfolder)
            if logger:
                logger.info("Temporary folder removed.")
        exit(0)
    if dict_option['ACTION'] == 'recover':
        if logger:
            logger.info("ACTION is 'recover'.")
        last_stage = detect_stage_last_finished(recovername)
        if last_stage:
            write_formatted_string_withtime("The pipeline was stopped after stage '" + last_stage + "'.", 30, sys.stdout)
        if last_stage == "prepare":
            write_formatted_string("*** Starting the pipeline from stage 'candidate'.", 30, sys.stdout)
            run_candidate(dict_option, outtempfolder, recovername)
            run_fold(dict_option, outtempfolder, recovername)
            if not dict_option["KEEPTMP"]:
                run_removetmp(outtempfolder)
                if logger:
                    logger.info("Temporary folder removed.")
        elif last_stage == "candidate":
            write_formatted_string("*** Starting the pipeline from stage 'fold'.", 30, sys.stdout)
            run_fold(dict_option, outtempfolder, recovername)
        elif last_stage == "fold":
            write_formatted_string("*** Starting the pipeline from stage 'predict'.", 30, sys.stdout)
            run_predict(dict_option, outtempfolder, recovername)
            if not dict_option["KEEPTMP"]:
                run_removetmp(outtempfolder)
                if logger:
                    logger.info("Temporary folder removed.")
        elif last_stage == "predict":
            write_formatted_string("*** The pipeline has been finished on the input. No action is performed.\n", 30, sys.stdout)
        else:
            write_formatted_string_withtime("No recovery information found. Please run the pipeline in the 'pipeline' mode.\n", 30, sys.stdout)
        exit(0)
