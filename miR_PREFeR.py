import sys
import os.path
import multiprocessing
import Queue
import cPickle
import time
import re
import subprocess
import os
import distutils.spawn
import tempfile
import string
import gc
import logging
import shutil


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
    parser.add_option("-L", "--log", action="store_true", dest="log",
                      help="Generate a log file.")
    parser.add_option("-k", "--keep-tmp", action="store_true", dest="keeptmp",
                      help="After finish the whole pipeline, do not remove the temporary folder that contains the intermidate files.")
    parser.add_option("-d", "--output-detail-for-debug", action='store_true', dest='debugmode')
    actions = ['check','prepare','candidate','fold', 'predict', 'pipeline', 'recover']
    (options, args) = parser.parse_args()
    if len(args) != 2:
         parser.error("incorrect number of arguments. Run the script with -h option to see help.")
    if args[0] not in actions:
        parser.error("unknown command")
    dict_option = parse_configfile(args[1])
    dict_option['ACTION'] = args[0]
    dict_option["LOG"] = options.log
    dict_option["KEEPTMP"] = options.keeptmp
    dict_option['OUTPUT_DETAILS_FOR_DEBUG'] = False
    if options.debugmode:
        dict_option['OUTPUT_DETAILS_FOR_DEBUG'] = options.debugmode
    return dict_option

def parse_configfile(configfile):
    if not os.path.exists(os.path.expanduser(configfile)):
        sys.stderr.write("Configuration file " + configfile + " does not exist!!\n")
        sys.exit(-1)
    dict_option = {
        "CONFIG_FILE":configfile,
        "FASTA_FILE":"",
        "ALIGNMENT_FILE":[],
        "GFF_FILE_EXCLUDE": "",
        "GFF_FILE_INCLUDE": "",
        "PRECURSOR_LEN":300,
        "READS_DEPTH_CUTOFF":10,
        "MAX_GAP":100,
        "NUM_OF_CORE":1,
        "OUTFOLDER":"./",
        "TMPFOLDER": "",
        "NAME_PREFIX":"",
        "PIPELINE_PATH":"",
        "DELETE_IF_SUCCESS":"Y",
        "CHECKPOINT_SIZE":3000,  #checkpoint every 3000 sequences
        "MIN_MATURE_LEN": 18,
        "MAX_MATURE_LEN": 23,
        "ALLOW_3NT_OVERHANG": False,
        "ALLOW_NO_STAR_EXPRESSION": True,
        "OUTPUT_DETAILS_FOR_DEBUG": False
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
                        if not os.path.exists(os.path.expanduser(name.strip())):
                            sys.stderr.write("File " + name.strip() +
                                             " does not exist!!\n")
                            sys.exit(-1)
                        dict_option[key].append(os.path.expanduser(name.strip()))
                    continue
                if key == "GFF_FILE_EXCLUDE" or key == "GFF_FILE_INCLUDE" or key == "FASTA_FILE":
                    if not os.path.exists(os.path.expanduser(sp[1].strip())):
                        sys.stderr.write("File " + sp[1].strip() +
                                          " does not exist!!\n")
                        sys.exit(-1)
                    dict_option[key] = os.path.expanduser(sp[1].strip())
                    continue
                if key == "NUM_OF_CORE":
                    cpucount = multiprocessing.cpu_count()
                    cpu_to_use = int(sp[1].strip())
                    if 2*cpucount < cpu_to_use:
                        sys.stderr.write(
                            "Warnning: 2*NUM_OF_CORE is larger than CPUS/Cores on" +
                            " the machine. Use "+str(2*cpucount)+" instead.\n")
                        cpu_to_use = 2*cpucount
                    dict_option[key] = cpu_to_use
                    continue
                if key == "PRECURSOR_LEN" or key == "MAX_MATURE_LEN" or key == "MIN_MATURE_LEN":
                    dict_option[key] = int(sp[1].strip())
                    continue
                if key == "READS_DEPTH_CUTOFF" or key == 'MAX_GAP' or key == 'CHECKPOINT_SIZE':
                    dict_option[key] = int(sp[1].strip())
                    continue
                if key =="OUTFOLDER" or key== "TMPFOLDER" or key == "NAME_PREFIX" or key == "DELETE_IF_SUCCESS":
                    dict_option[key] = sp[1].strip()
                    continue
                if key == 'ALLOW_NO_STAR_EXPRESSION' or key == 'ALLOW_3NT_OVERHANG':
                    value = sp[1].strip().lower()
                    if value != 'y' and value != 'n':
                        sys.stderr.write("Value for ALLOW_3NT_OVERHANG/ALLOW_NO_STAR_EXPRESSION must be one of 'Y/y/N/n'.\n")
                        sys.exit(-1)
                    if value == 'y':
                        dict_option[key] = True
                    else:
                        dict_option[key] = False
                    continue
                if key == "PIPELINE_PATH":
                    if not os.path.exists(os.path.expanduser(sp[1].strip())):
                        sys.stderr.write("miR-PREFeR path " + sp[1].strip() +
                        " does not exist!!\n")
                        sys.exit(-1)
                    dict_option[key] = os.path.expanduser(sp[1].strip())
                    continue
    allgood = True
    if dict_option["PRECURSOR_LEN"] < 60 or dict_option["PRECURSOR_LEN"] > 3000:
        sys.stderr.write("Error: allowed precursor range: 60-3000\n")
        allgood = False
    if dict_option['MIN_MATURE_LEN'] > dict_option['MAX_MATURE_LEN']:
        sys.stderr.write("Error: MIN_MATURE_LEN is greater than MAX_MATURE_LEN.\n")
        allgood = False
    if dict_option["READS_DEPTH_CUTOFF"] < 2:
        sys.stderr.write("Error: READS_DEPTH_CUTOFF should >=2.\n")
        allgood = False
    if dict_option["CHECKPOINT_SIZE"] < 10:
        sys.stderr.write("Error: CHECKPOINT_SIZE should >=10.\n")
        allgood = False
    if dict_option['GFF_FILE_INCLUDE'] and dict_option['GFF_FILE_EXCLUDE']:
        sys.stderr.write("Error: GFF_FILE_EXCLUDE and GFF_FILE_INCLUDE are mutual exclusive, please remove one of them.\n")
        allgood = False
    if not allgood:
        exit(-1)

    return dict_option


def display_dict_option(dict_option):
    write_formatted_string_withtime("", 30, sys.stdout)
    sys.stdout.write("============================================================================\n")
    sys.stdout.write("===============Configurations for this experiment:==========================\n")
    sys.stdout.write("============================================================================\n")
    for k in sorted(dict_option.keys()):
        print(k+": "+str(dict_option[k]))
    sys.stdout.write("============================================================================\n")
    sys.stdout.write("============================================================================\n\n")


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


def get_read_depth_fromID_as_string(seqid):
    '''
    Read ID should be in 'samplename_rA_xN' format. Where A is an number
    identify the read, and N is the depth of the read.
    '''

    pattern=r"^\S+_x([0-9]+)"
    m = re.match(pattern, seqid)
    if m:
        return  m.groups()[0]
    else:
        raise Exception('Read Id format is not right. Read id must be in "samplename_rA_xN" format. See readme section "Prepare input data for the pipeline."\n')



def check_gff(gffname):
    '''
    Check the input gff file is in the right format.
    We did not check all the requirements made by the GFF3 format, but only
    those fields that are essential for our pipeline.

    Note that fields on each line MUST be separated with one TAB, not SPACE. This
    is a requirement of the GFF3 format.

    The return value is a tuple with two elements. The first one is True/False.
    if True, the second is a set containing all the IDs of the sequences; if
    False, the second is a error message.
    '''

    set_ids = set()
    good_count = 0
    with open(gffname) as f:
        for idx, line in enumerate(f):
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    if good_count == 0:
                        return (False, "No GFF entry in the given gff file.")
                    else:
                        return (True, set_ids)
                else:
                    continue
            else:
                if not line.strip():  #empty lines
                    continue
                sp = line.strip().split('\t')
                if len(sp) != 9:
                    return (False, "Line "+str(idx+1) +" does not have " +
                            str(len(sp))+
                            " fields, 9 required. Make sure fields on one line is separated by one TAB, not SPACEs.")
                else:
                    try:
                        start = int(sp[3])
                    except ValueError:
                        return (False, "Column 4 on line "+str(idx+1) +" should be an integer.")
                    try:
                        end = int(sp[4])
                    except ValueError:
                        return (False, "Column 5 on line "+str(idx+1) +" should be an integer.")
                    if end < start:
                        return (False, "End position is smaller than start position on line "+str(idx+1) +".")
                    else:
                        set_ids.add(sp[0])
                        good_count += 1
    return (True, set_ids)


def check_reference(refname):

    '''
    Check the reference file. The reference file should be a multi-fasta file.

    The return value is a tuple with two elements. The first one is True/False.
    if True, the second is a set containing all the IDs of the sequences; if
    False, the second is a error message.

    Because miR-PREFeR uses samtools, and samtools requires all sequence lines
    of an entry have the same length except the last line, all sequence ID
    lines do not contain white space. So we check the following:
    1.  The ID of each sequence is unique. Each ID line starts with ">". The ID
    should no white space characters. Everything after
    the first blank character is ignored(ID only count to the first non-blank
    character).
    2.  There is no black lines in between the sequence line of each entry.
    3.  Each sequence line of the same entry is of fixed length, except the last
    line of each entry.
    '''

    command = "samtools faidx " + refname
    try:
        subprocess.check_call(command.split())
    except subprocess.CalledProcessError:
        return (False, 0)

    # idpattern = r'^>\s*([-_a-zA-Z0-9]+)\s*.*$'
    idpattern = r'^>\s*(\S+)\s*.*$'
    dict_refID = {}
    with open(refname) as f:
        refID = ""
        for idx, line in enumerate(f):
            if line.startswith(">"):
                m = re.match(idpattern, line.rstrip())
                if m:
                    refID = m.groups()[0]
                    if refID in dict_refID:  #duplicate ID
                        return (False, "More than one sequence have the same sequence ID: " + refID)
                    else:
                        dict_refID[refID] = 1
                else:
                    return (False, "Sequence ID at line " + str(idx+1) + " is not correct. Please make sure that the sequence ID line starts with '>'. The sequence ID starts from the first non-blank character after '>' and counts until the first blank character. The sequence ID only contains a-z, A-Z, 0-9, underscore(_), and dash(-).")
            else:
                continue
    return (True, dict_refID)


def check_sam_format(samname):
    '''
    This function does the following:
    1. Check the read IDs are in right format (see get_read_depth_fromID_as_string)
    2. Check whether the SAM header is present. The pipeline requires the SAM
    file has the SAM header.
    3. Check the flag column of the SAM file. The pipeline requires the second
    column (the flag column) of the file is in 10-based integer format. SAMtools
    can output SAM file which has a Hex or string format flag column (-X and -x
    options in SAMtools view command).

    This only checks the first 10000 alignment lines of the SAM file. If all
    10000 lines are right, then it assumes the file is well-formed.

    '''

    def check_sam_flag(flag):
        if re.match(r'^[0-9]+$', flag):
            return True
        return False

    pattern=r"^(\S+)_r[0-9]+_x([0-9]+)$"
    has_header = False
    seqid_right = True
    flag_right = True
    counter = 0
    samplename = ""
    with open(samname) as f:
        line = f.readline()
        if line.startswith("@"):
            has_header = True
        for idx, line in enumerate(f):
            if line.startswith("@"):
                continue
            else:
                if not samplename:
                    samplename = "_".join(line.split()[0].split("_")[0:-2])
                counter += 1
                if counter==10000:
                    break
                sp = line.split()
                m = re.match(pattern, sp[0])
                if m:  # seqid right
                    if m.groups()[0] == samplename:
                        continue
                    else:
                        write_formatted_string("ERROR: ReadID not right on line " + str(idx+1) + ", ReadID=" +sp[0] + ", samplename=" + samplename , 30, sys.stderr)
                        seqid_right = False
                else:  # seqid wrong
                    seqid_right = False
                if check_sam_flag(sp[1]):  # flag right
                    continue
                else:  #
                    flag_right = False
    return (has_header, seqid_right, flag_right)


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
    if ret: #check the depth, faidx, view command exist.
        command = "samtools view"
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
        if outerr.find("unrecognized") !=-1:
            message = "Samtools view command does not exist. Please make sure samtools is correctly installed and the version is > 0.1.15. Refer to the README for more infomation"
            return False, message
        command = "samtools depth"
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
        if outerr.find("unrecognized") !=-1:
            message = "SAMtools depth command does not exist. Please make sure samtools is correctly installed and the version is > 0.1.15. Refer to the README for more infomation"
            return False, message
        command = "samtools faidx"
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
        if outerr.find("unrecognized") !=-1:
            message = "SAMtools faidx command does not exist. Please make sure samtools is correctly installed and the version is > 0.1.15. Refer to the README for more infomation"
            return False, message
        command = "samtools index"
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
        if outerr.find("unrecognized") !=-1:
            message = "SAMtools index command does not exist. Please make sure samtools is correctly installed and the version is > 0.1.15. Refer to the README for more infomation"
            return False, message

        command = "samtools sort"
        check_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = check_process.communicate()
        if outerr.find("unrecognized") !=-1:
            message = "SAMtools sort command does not exist. Please make sure samtools is correctly installed and the version is > 0.1.15. Refer to the README for more infomation"
            return False, message

        command = "samtools"
        version_num = "Unknown"
        version_process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = version_process.communicate()
        m=re.search(r'([Vv]ersion.+)\n',outerr)
        if m:
            version_num = m.groups()[0]
        return True, version_num
    else:
        return False, "SAMtools not installed or not in the PATH."


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
        return "Unknown"
    if outmessage == "":
        return "Unknown"
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
        sys.stderr.write("Error occurred when indexing the genome file\n")
        sys.exit(-1)


def gen_temp_gff(gffname, tmpdir):
    tempgffname = os.path.join(tmpdir, "temp.remove.unsort.gff")
    fout = open(tempgffname, 'w')
    with open(gffname) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                else:
                    continue
            else:
                if not line.strip():
                    continue
                else:
                    fout.write(line)
    fout.close()
    return tempgffname


def gen_keep_regions_from_exclude_gff(gffname, tmpdir, dict_len, minlen):
    '''
    Generate a BED format file which contains all the regions that are not
    overlap with the regions of the features in the gff file.

    Note that gff file is 1 based, and the ending position is inclusive.
    BED file is 0 based, and the ending position is exclusive.
    '''

    def overlap(r1, r2):
        '''
        If r1 and r2 overlaps, return the combined region, else, return None.
        '''
        if r1[1] < r2[0] or r1[0] > r2[1]:
            return None
        else:
            return (min(r1[0],r2[0]), max(r1[1],r2[1]))

    tempgffunsortname = gen_temp_gff(gffname, tmpdir)
    #sort the gff file
    tempgffname = os.path.join(tmpdir, "temp.remove.gff")
    tempbedname = os.path.join(tmpdir, "temp.keepregion.bed")
    command = "sort -k1,1 -k4,4n " + tempgffunsortname + " -o " + tempgffname
    try:
        p = subprocess.Popen(command.split())
        while True:
            retcode = p.poll()
            if retcode is None:
                write_formatted_string_withtime("Sorting GFF file.", 30, sys.stdout)
                time.sleep(5)
            else:
                write_formatted_string_withtime("GFF file sorting done.", 30, sys.stdout)
                break
    except Exception as e:
        sys.stderr.write("Error occurred when sorting the GFF file.\n")
        sys.stderr.write("Exception message: " + str(e) + "\n")
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
                        foutput.write(seqid+"\t"+str(region[1]-1)+"\t"+str(dict_len[seqid])+"\n")
                if int(sp[3]) >=minlen:
                    foutput.write(cur_seqid+"\t0"+"\t"+sp[3]+"\n")
                seqid = cur_seqid
                region = (int(sp[3]), int(sp[4]))
            else:
                ret = overlap(region, (int(sp[3]), int(sp[4])))
                if ret:
                    region = ret
                else:
                    if int(sp[3])-region[1] >=minlen:
                        foutput.write(cur_seqid+"\t"+str(region[1]-1)+"\t"+sp[3]+"\n")
                    region = (int(sp[3]), int(sp[4]))
                    seqid = cur_seqid
            #write the last one
        if dict_len[seqid] - region[1] >= minlen:
            foutput.write(seqid+"\t"+str(region[1]-1)+"\t"+str(dict_len[seqid])+"\n")
        foutput.close()
    return tempbedname


def gen_keep_regions_from_include_gff(gffname, tmpdir, minlen):
    '''
    Generate a BED format file which contains all the regions that are
    overlap with the regions of the features in the gff file.

    Note that gff file is 0 based, and the ending position is inclusive. BED
    file is 1 based, and the ending position is exclusive.
    '''
    tempgffunsortname = gen_temp_gff(gffname, tmpdir)
    #sort the gff file
    tempgffname = os.path.join(tmpdir, "temp.remove.gff")
    tempbedname = os.path.join(tmpdir, "temp.keepregion.bed")
    command = "sort -k1,1 -k4,4n " + tempgffunsortname + " -o " + tempgffname
    try:
        p = subprocess.Popen(command.split())
        while True:
            retcode = p.poll()
            if retcode is None:
                write_formatted_string_withtime("Sorting GFF file.", 30, sys.stdout)
                time.sleep(5)
            else:
                write_formatted_string_withtime("GFF file sorting done.", 30, sys.stdout)
                break
    except Exception as e:
        sys.stderr.write("Error occurred when sorting the GFF file.\n")
        sys.stderr.write("Exception message: " + str(e) + "\n")
        sys.exit(-1)
    foutput = open(tempbedname, 'w')
    with open(tempgffname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            sp = line.split()
            if int(sp[4])-int(sp[3]) >=minlen:
                foutput.write(sp[0]+"\t"+str(int(sp[3])-1)+"\t"+sp[4]+"\n")
        foutput.close()
    return tempbedname



def sam2bam(samfile, bamfile):
    command = "samtools view -bS -o " + bamfile + " " + samfile
    try:
        write_formatted_string_withtime("Command: "+command, 30, sys.stdout)
        subprocess.check_call(command.split())
        write_formatted_string_withtime("Done sam2bam.", 30, sys.stdout)
    except Exception as e:
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
        sys.stderr.write("Error occurred when filtering regions and sorting the BAM file\n")
        sys.exit(-1)
    return outbamprefix+".sort.bam"


def expand_bamfile(bamfile, maxdepth, outputbamfile, outputsamfile, tmpdir):
    #tempsamfile = tempfile.NamedTemporaryFile(mode='w',prefix = str(os.getpid())+"tempsam", suffix=".sam", delete=False)
    tempsamname = os.path.join(tmpdir, str(os.getpid())+".tempsam.sam")
    tempsamfile = open(tempsamname, 'w')
    command = "samtools view -h -o " + tempsamfile.name + " " + bamfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
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
            depth = int(sp[0].split("_")[-1].lstrip("x"))
            if depth > maxdepth:
                depth = maxdepth
            for i in xrange(depth):
                outf.write(line)
        outf.close()
    command = "samtools view -bS -o" + outputbamfile + " " + outputsamfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
        sys.stderr.write("Error occurred when converting SAM to BAM\n")
        sys.exit(-1)
    return outputbamfile, outputsamfile


def index_bam(bamfile, outindexfile):
    command = "samtools index " + bamfile + " "+outindexfile
    try:
        subprocess.check_call(command.split())
    except Exception as e:
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
        sys.stderr.write("Error occurred when filtering BAM file using flags.\n")
        sys.exit(-1)
    return outbamname


def prepare_data(dict_option, outtempfolder, logger):

    if logger:
        logger.info("Getting genomic sequence lengths.")
    write_formatted_string_withtime("Getting genomic sequence lengths.", 30, sys.stdout)
    #get the length of all the genomic sequences in the fasta/alignment files
    dict_len = get_length_from_sam(dict_option["ALIGNMENT_FILE"][0])

    #if the fasta file is not index, then create a index
    if not os.path.exists(dict_option["FASTA_FILE"]): #index the genome
        index_command = "samtools faidx " + dict_option["FASTA_FILE"]
        try:
            write_formatted_string_withtime("Indexing reference seqeunces.", 30, sys.stdout)
            write_formatted_string("Command: "+index_command, 30, sys.stdout)
            subprocess.check_call(index_command.split())
        except Exception as e:
            if logger:
                logger.error("Error occurred when indexing the genome file. "+
                             "Command: "+index_command)
            sys.stderr.write("Error occurred when indexing the genome file\n")
            sys.exit(-1)

    #convert the SAM files(generated by Bowtie) to BAM files
    if logger:
        logger.info("Converting SAM files to BAM files using samtools.")
    write_formatted_string_withtime("Converting SAM files to BAM files using SAMtools.", 30, sys.stdout)
    bamfiles = []
    for samname in dict_option["ALIGNMENT_FILE"]:
        if logger:
            logger.info("Converting "+samname)
            write_formatted_string_withtime("Converting "+samname, 30, sys.stdout)
        bamname = os.path.join(outtempfolder, os.path.basename(samname)+".bam")
        sam2bam(samname, bamname)
        bamfiles.append(bamname)

    combinedbamname = os.path.join(outtempfolder, "combined.bam")
    if len(bamfiles) > 1:
        #combine multiple BAM files from multiple sample together
        if logger:
            logger.info("Combining multiple BAM files from multiple samples together.")
        write_formatted_string_withtime("Combining multiple BAM files from multiple samples together.", 30, sys.stdout)
        combine_bamfiles(dict_option["ALIGNMENT_FILE"][0], combinedbamname, *bamfiles)
    else:
        shutil.copyfile(bamfiles[0], combinedbamname)

    #removing reads that are overlapped with features in the gff file, if provided.
    if os.path.exists(dict_option["GFF_FILE_EXCLUDE"]):
        #TODO: minlen should be user adjustable, not fix here.
        write_formatted_string_withtime("Removing reads that are overlapped with features in the GFF file.", 30, sys.stdout)
        if logger:
            logger.info("Removing reads that are overlapped with features in the GFF file.")

        tempkeepregion = gen_keep_regions_from_exclude_gff(dict_option["GFF_FILE_EXCLUDE"], outtempfolder, dict_len, 55)
        num_lines = 0
        with open(tempkeepregion) as f:
            for line in f:
                num_lines += 1
                if num_lines > 5:
                    break
        if num_lines==0:
            write_formatted_string_withtime("!!! No regions need to analyze after excluding regions in the GFF file, stop analyze!", 30, sys.stdout)
            if logger:
                logger.info("!!! No regions need to analyze after excluding regions in the GFF file, stop analyze!")
            exit(-1)

        combinedbamname = gen_keep_regions_sort_bam(combinedbamname, tempkeepregion, os.path.join(outtempfolder,"combined.filtered"))
    elif os.path.exists(dict_option["GFF_FILE_INCLUDE"]):
        write_formatted_string_withtime("GFF_FILE_INCLUDE specified, removing reads that are not overlap with features in the GFF file.", 30, sys.stdout)
        if logger:
            logger.info("GFF_FILE_INCLUDE specified, removing reads that are not overlap with features in the GFF file.")

        tempkeepregion = gen_keep_regions_from_include_gff(dict_option["GFF_FILE_INCLUDE"], outtempfolder,  55)
        num_lines = 0
        with open(tempkeepregion) as f:
            for line in f:
                num_lines += 1
                if num_lines > 5:
                    break
        if num_lines==0:
            write_formatted_string_withtime("!!! No regions in the GFF_FILE_INCLUDE file or all regions are shorter than 55, stop analyze!", 30, sys.stdout)
            if logger:
                logger.info("!!! No regions in the GFF_FILE_INCLUDE file or all regions are shorter than 55, stop analyze!")
            exit(-1)

        combinedbamname = gen_keep_regions_sort_bam(combinedbamname, tempkeepregion, os.path.join(outtempfolder,"combined.filtered"))

    else:
        combinedbamname = sort_index_bam(combinedbamname, os.path.join(outtempfolder,"combined.filtered"))

    expandedsamname = os.path.join(outtempfolder, "expanded.sam")
    expandedbamname = os.path.join(outtempfolder, "expanded.bam")
    expandedbam_plus = os.path.join(outtempfolder, "expanded.plus.bam")
    expandedbam_minus = os.path.join(outtempfolder, "expanded.minus.bam")
    if logger:
        logger.info("Generating expanded BAM and SAM files")
    expand_bamfile(combinedbamname, dict_option["READS_DEPTH_CUTOFF"], expandedbamname, expandedsamname, outtempfolder)
    if logger:
        logger.info("Generating "+expandedbam_plus)
    filter_bam_by_flag(expandedbamname, 16, expandedbam_plus, keep=False)
    if logger:
        logger.info("Generating "+expandedbam_minus)
    filter_bam_by_flag(expandedbamname, 16, expandedbam_minus, keep=True)
    return combinedbamname, expandedsamname, expandedbamname, expandedbam_plus, expandedbam_minus


def gen_contig_typeA(expandedbam_plus, expandedbam_minus, dict_option,
                     contig_minlen, logger, tmpdir):
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
            sys.stderr.write("Error occurred when generating contigs(typeA).\n")
            sys.exit(-1)

    depthfilename = os.path.join(tmpdir, "bam.depth.cut"+str(dict_option["READS_DEPTH_CUTOFF"]))
    f=open(depthfilename,'w')
    f.write(output.decode())
    f.close()
    if logger:
        logger.info("Generating contigs.")
    dict_contigs = {}
    cnt = 0
    for region in get_next_non_zero_region(depthfilename):

        if region[2] - region[1] < contig_minlen:
            continue
        cnt += 1
        if cnt % 20000 ==0:
            write_formatted_string_withtime(" .. generating peaks .." , 30, sys.stdout)
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


def dump_loci_seqs_and_alignment_multiprocess(dict_loci, piece_info_list,
                                              sortedbamname, fastaname,
                                              outputseqprefix, outputalnprefix,
                                              min_mature_depth,
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
                    if num_seq % 3000 == 0:
                        write_formatted_string_withtime(" .. generating candidate sequences, output file name: " + outputfastaname, 30, sys.stdout)
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
                        cur_strand = "+"
                        if "+" in strands:
                            seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo  # >ExtendRegion strand loci 0/L/R peaks
                            dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                        else:
                            cur_strand = "-"
                            seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo
                            dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                            seq = get_reverse_complement(seq)

                        alignments = samtools_view_region(sortedbamname, seqid,
                                                          extendregion[0][0],
                                                          extendregion[0][1])
                        dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                        matures = gen_possible_matures_loci(dict_loci_info, extendregion, cur_strand, min_mature_depth)
                        cPickle.dump([dump_key, extendregiontag.strip(), dict_loci_info, matures], foutdump, protocol=2)
                        for m in matures:
                            seqtag += " " + "M:" + str(m[0])+"-"+str(m[1]) + "/" + str(m[2]) + "/" + str(m[3])
                        seqtag += "\n"
                        fout.write(seqtag)
                        fout.write(seq+"\n")
                        num_seq += 1
                        continue
                    if len(strands) == 2:  # this region has peaks on both strands, generate precursors for each strand
                        alignments = samtools_view_region(sortedbamname, seqid,
                                                          extendregion[0][0],
                                                          extendregion[0][1])
                        dict_loci_info = gen_loci_alignment_info(alignments, seqid, extendregion)
                        matures_plus = gen_possible_matures_loci(dict_loci_info, extendregion, "+", min_mature_depth)
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "+"]
                        plus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures_plus])
                        matures_minus = gen_possible_matures_loci(dict_loci_info, extendregion, "-", min_mature_depth)
                        dump_key = [seqid, (extendregion[0][0],extendregion[0][1]), "-"]
                        minus_dumpinfo.append([dump_key, extendregiontag.strip(), dict_loci_info, matures_minus])

                        otherinfo_plus = ""
                        for peak in extendregion[1]:
                            if peak[2] == "+":
                                otherinfo_plus = otherinfo_plus + str(peak[0])+","+str(peak[1])+","+str(peak[2])+";"
                        otherinfo_plus = otherinfo_plus.rstrip(";")
                        seqtag = ">"+regionpos+" + "+ loci_pos + extendregiontag + otherinfo_plus
                        for m in matures_plus:
                            seqtag += " " + "M:" + str(m[0])+"-"+str(m[1]) + "/" + str(m[2]) + "/" + str(m[3])
                        seqtag += "\n"

                        plus_tag.append(seqtag)
                        plus_seq.append(seq)
                        num_seq += 1
                        otherinfo_minus = ""
                        for peak in extendregion[1]:
                            if peak[2] == "-":
                                otherinfo_minus = otherinfo_minus + str(peak[0])+","+str(peak[1])+","+str(peak[2])+";"
                        otherinfo_minus = otherinfo_minus.rstrip(";")
                        seqtag = ">"+regionpos+" - "+ loci_pos + extendregiontag + otherinfo_minus
                        for m in matures_minus:
                            seqtag += " " + "M:" + str(m[0])+"-"+str(m[1]) + "/" + str(m[2]) + "/" + str(m[3])
                        seqtag += "\n"
                        seq = get_reverse_complement(seq)
                        minus_tag.append(seqtag)
                        minus_seq.append(seq)
                        num_seq += 1

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
        write_formatted_string_withtime(" Done: generating candidate sequences, output file name: " + outputfastaname, 30, sys.stdout)
        queue.put(((outputfastaname, outputalnname), num_loci, num_seq))
        queue.put("done")
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
    num_joined = 0
    while True:
        try:
            if num_joined == len(jobs):
                for job in jobs:
                    job.join()
                break
            info = inforqueue.get_nowait()
            if info == "done":
                num_joined += 1
            else:
                total_loci += info[1]
                total_seq += info[2]
                finalresult.append(info[0])
        except Queue.Empty:
            time.sleep(2)
            continue
    # for job in jobs:
    #     job.join()
    #     info = inforqueue.get()
    #     total_loci += info[1]
    #     total_seq += info[2]
    #     finalresult.append(info[0])
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
            extend_both_side = (max_transcript - length) / 2
            left = region[0] - extend_both_side
            right = region[1] + extend_both_side
            if left < 0:
                left = 1
            if right > seqlen:
                right = seqlen + 1
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
                if total_loci%20000 == 0:
                    write_formatted_string_withtime(" .. combining peaks ..", 30, sys.stdout)

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
    # dict_pos_info = {"+": total count
    #                   "-":total cound
    #                   0: {dict_map}
    # }
    # dict_map = { pos1: {dict_strand}
    #             ....
    #             posn: {dict_strand}
    #}
    # dict_strand = { '+': [max_depth_list, [pos_info, ..., pos_info]]
    # }
    # max_depth_list = [readlen, this_depth, total_depth]
    # pos_info = [readlen, depth]
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
        #depth = int(sp[0].split("-")[-1])  # read names must be in "name-depth" format
        # read id could be in "name-N", "name_N","name-xN", or "name_xN" format
        depth = int(get_read_depth_fromID_as_string(sp[0]))
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


def gen_possible_matures_loci(dict_loci_info, loci, strand, min_mature_depth):
    def gen_matures_one_peak(max_num_mature, min_mature_depth, peak_start, peak_end):
        peak_start = peak_start - 20  # the most abundant read in the peak could start before the peak
        peak_end = peak_end
        matures = []
        higest = (0, 0, 0, 0)
        for pos in xrange(peak_start, peak_end):
            if pos in dict_loci_info[0]:
                if strand in dict_loci_info[0][pos]:
                    if dict_loci_info[0][pos][strand][0][1] > higest[3]:
                        higest = (pos, pos+dict_loci_info[0][pos][strand][0][0], strand, dict_loci_info[0][pos][strand][0][1])
                    if dict_loci_info[0][pos][strand][0][1] > min_mature_depth:
                        if len(matures) < max_num_mature:
                            matures.append((pos, pos+dict_loci_info[0][pos][strand][0][0], strand, dict_loci_info[0][pos][strand][0][1]))
                        elif dict_loci_info[0][pos][strand][0][1] > matures[-1][-1]:
                            matures[-1] = (pos, pos+dict_loci_info[0][pos][strand][0][0], strand, dict_loci_info[0][pos][strand][0][1])
                        matures.sort(key=lambda x: x[2], reverse=True)
        if not matures:
            matures.append(higest)
        return matures

    ret_matures = []
    for start, end, curstrand in loci[1]:
        if curstrand != strand:
            continue
        # mature_depth = 0
        # mature_length = 0
        # mature_pos = 0
        # for pos in xrange(start,end):
        #     if pos in dict_loci_info[0]:
        #         if strand in dict_loci_info[0][pos]:
        #             if dict_loci_info[0][pos][strand][0][1] > mature_depth:
        #                 mature_depth = dict_loci_info[0][pos][strand][0][1]
        #                 mature_length = dict_loci_info[0][pos][strand][0][0]
        #                 mature_pos = pos
        # mature = (mature_pos, mature_pos+mature_length, strand, mature_depth)
        #ret_matures.append(mature)
        matures = gen_matures_one_peak(2, min_mature_depth, start, end)
        ret_matures.extend(matures)
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

def get_structures_next_extendregion(rnalfoldoutname, minlen,
                                     minloop=3):
    '''
    Parse the RNALfold output file, return a tuple of (0/L/R, peak region, [
    norm_energy, fold_start, ss, sstype]) for an extend region each time (generator).

    Structure whose length is shorter than minlen are filtered.

    '''

    energypattern = re.compile(r'\(\s*(-?[0-9]+[\.]?[0-9]*)\s*\)')
    structures = []
    which = 0
    peak = ''
    first = True
    sp=[]
    with open(rnalfoldoutname) as f:
        for line in f:
            sp = line.strip().split()
            if line.startswith(">"):
                if not first:
                    yield (which, peak, structures)
                structures = []
                which = sp[3]  # 0, L, or R
                peak = sp[2]
                first = False
            else:
                if len(sp) >=3:  # structure line
                    if len(sp[0]) < minlen:
                        continue
                    norm_energy = energypattern.search(line).group(1)
                    norm_energy = float(norm_energy)/len(sp[0])
                    sstype = -1
                    if is_stem_loop(sp[0], minloop):
                        sstype = 0  # stemploop ()
                        structures.append((norm_energy,int(sp[-1]), sp[0], sstype))
                    else:
                        subss, totalout = filter_ss(sp[0])
                        if not subss:
                            continue
                        for cur_ss in subss:  #start, strcture
                            if is_stem_loop(cur_ss[1], minloop):
                                sstype = 0  # stemploop ()
                                structures.append((norm_energy,int(sp[-1])+cur_ss[0], cur_ss[1], sstype))
                                continue
                            if has_one_good_bifurcation(cur_ss[1]):
                                sstype = 1
                                #  TODO: use whole structure energy as sub ss energy now.
                                structures.append((norm_energy,int(sp[-1])+cur_ss[0], cur_ss[1], sstype))
                        # if has_one_good_bifurcation(sp[0]):
                        #     sstype = 1  # good bifurcation (()())
                        # else:
                        #     if two_parallel_stems(sp[0]):
                        #         sstype = 2  # ()()
                        #     else:
                        #         sstype = -1
                else:
                    continue
        yield (which, peak, structures)


def is_stem_loop(ss, minloopsize):
    last_open = ss.rfind("(")
    first_close = ss.find(")")
    if first_close - last_open-1 >= minloopsize:
        return True
    else:
        return False


def has_one_good_bifurcation(ss):
    # not an stem loop, but a structure with two bifurcation: (()()).
    bifurcation = 0
    stack = []
    last = 'Push'
    stack_pos = []
    first_bi_last = 0
    second_bi_first = 0
    last_pos = 0
    dict_pos = {}
    for idx, char in enumerate(ss):
        if char == '(':
            if idx !=0:
                if not stack:
                    if last == 'Pop': # ()() like structure
                        return False
                    last = "Push"
                    stack.append(char)
                    stack_pos.append(idx)
                    last_pos = idx
                else: # stack is not empty
                    if last == "Pop":
                        if bifurcation >= 1:
                            return False
                        else:
                            bifurcation = 1
                            first_bi_last = last_pos
                            second_bi_first = idx
                    stack.append(char)
                    last = "Push"
                    stack_pos.append(idx)
                    last_pos = idx
            else:
                stack.append(char)
                last = "Push"
                stack_pos.append(idx)
                last_pos = idx
        elif char == ")":
            last = "Pop"
            stack.pop()
            v = stack_pos.pop()
            dict_pos[idx] = v
            dict_pos[v] = idx
            last_pos = idx
    if float(dict_pos[second_bi_first]-dict_pos[first_bi_last])/len(ss) < 0.5:
        if float(dict_pos[first_bi_last])/len(ss) > 0.25:
            if float(dict_pos[second_bi_first])/len(ss) < 0.75:
                return True
    return False


def two_parallel_stems(ss):
    # not an stem loop, but a structure with two parallel stems: ()().
    stack = []
    last = ""
    stems = 0
    for idx, char in enumerate(ss):
        if char== ')':
            stack.pop()
            last = ')'
            if len(stack) == 0:
                stems += 1
        if char=='(':
            if last==')':
                if len(stack)!=0:
                    return False
            last = '('
            stack.append(char)
    if stems == 2:
        return True
    return False



def filter_ss(ss):
    #input is not stemloop, use is_stem_loop to filter out stemloops
    stack = []
    dict_pair = {}
    for idx, char in enumerate(ss):
        if char == '(':
            stack.append(idx)
            continue
        if char == ')':
            first = stack.pop()
            dict_pair[first] = idx
            dict_pair[idx] = first
    result = []

    prefirst = ss.find('(')
    prelast = dict_pair[prefirst]
    curfirst = prefirst
    curlast = prelast
    pregap = [0, prefirst]
    aftergap=[]
    totaloutstems = 1  # the number of outmost stems
    while curfirst != -1:
        curfirst = ss.find('(', prelast)
        if curfirst == -1:  # the previous is the last one
            aftergap = [prelast+1, len(ss)]
        else:
            aftergap = [prelast+1, curfirst]
            curlast = dict_pair[curfirst]
            totaloutstems += 1
        if prelast - prefirst > 0:
            result.append((pregap[0],aftergap[1]))
        pregap = aftergap
        prelast = curlast
        prefirst = curfirst
    ret = []
    for s, e in result:
        length = e - s
        if length > 55:
            ret.append((s, ss[s:e]))
    return (ret, totaloutstems)


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


def stat_duplex(mature, star):
    dict_bp = {}
    ss = mature + star
    openpos = ss.find("(")
    closepos = ss.find(")")
    openchar = "("
    closechar = ")"
    bulges = []
    loops = []
    if openpos > closepos:
        openchar = ")"
        closechar = "("
    stack = []
    for idx, char in enumerate(ss):
        if char == openchar:
            stack.append(idx)
        elif char == closechar:
            pos = stack.pop()
            dict_bp[pos] = idx
    keys =  sorted(dict_bp.keys())
    for i in range(len(keys)-1):
        if (keys[i+1]-keys[i]==1) and (dict_bp[keys[i]]-dict_bp[keys[i+1]]==1):
            continue
        else:
            maturebulgesize = keys[i+1]-keys[i]-1
            starbulgesize = dict_bp[keys[i]]-dict_bp[keys[i+1]]-1
            if maturebulgesize == starbulgesize:
                loops.append(maturebulgesize)
            else:
                bulges.append((maturebulgesize, starbulgesize))
    return loops, bulges


def pass_stat_duplex(loops, bulges):
    total = len(loops) + len(bulges)
    if total > 5:
        return ["FAIL_STRUCTURE_TOO_MANY_BULGE_OR_LOOP", False]
    num_loop_bigger_than_three = 0
    num_bulge_bigger_than_three = 0
    totalloopsize = 0
    max_bulge = 0
    for size in loops:
        totalloopsize += size
        if size >= 4:
            num_loop_bigger_than_three += 1

    for size1, size2 in bulges:
        max_bulge = max(max_bulge, size1, size2)
        if size1>=4 or size2>=4:
            num_bulge_bigger_than_three += 1
    #if num_loop_bigger_than_three >1 or num_bulge_bigger_than_three >0:
    #    return False
    if max_bulge >2:  # maxmium bulge size is bigger than 2
        return ["FAIL_STRUCTURE_MAX_BULGE_LARGE_THAN_2", False]
    if totalloopsize >5:  # total loop size bigger than 5. >=5???
        return ["FAIL_STRUCTURE_TOTAL_LOOP_SIZE_LARGER_THAN_5", False]
    if len(bulges) > 2:  #number of bulges is bigger than 2
        return ["FAIL_STRUCTURE_NUM_BULGE_MORE_THAN_2", False]
    return [None, True]


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
                return "FAIL_STRUCTURE_MATCHED_BASES"
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
        return "FAIL_STRUCTURE_MATURE_NOT_IN_FOLD_REGION"

    mature_ss = ss[mature_local_pos[0]:mature_local_pos[1]]
    #  Not in a arm of a stem, return None. This should not happen given
    #  that the mature is already in an arm of a stem.
    if mature_ss.find("(")!=-1 and mature_ss.find(")")!=-1:
        return "FAIL_STRUCTURE_MATURE_NOT_IN_ONE_ARM"
    if len(mature_ss) - mature_ss.count(".") <14:
        return "FAIL_STRUCTURE_MATURE_MATCH_SMALL_THAN_14"


    #  local positions
    star_start = 0
    star_end = 0
    star_ss = ""
    mature_sym = "("

    firstbp = ss.find("(", mature_local_pos[0], mature_local_pos[1])
    lastbp =  ss.rfind("(", mature_local_pos[0], mature_local_pos[1])
    prime5 = True
    if firstbp == -1:
        firstbp = ss.find(")",mature_local_pos[0], mature_local_pos[1])
        lastbp =  ss.rfind(")",mature_local_pos[0], mature_local_pos[1])
        prime5 = False
        mature_sym = ")"
    if firstbp == -1:  # all are dots
        return "FAIL_STRUCTURE_MATURE_MATCH_SMALL_THAN_14"

    start_unpaired_len = mature_local_pos[1]-1-lastbp
    end_unpaired_len = firstbp-mature_local_pos[0]
    star_start = dict_bp[lastbp] -start_unpaired_len + 2
    star_end = dict_bp[firstbp] + end_unpaired_len + 2 + 1  # +1 because the end pos is exclusive

    #mature and star could overlap...
    if mature_local_pos[0] <= star_start:
        if star_start - mature_local_pos[1] < 3:
            return "FAIL_STRUCTURE_MATURE_STAR_OVERLAP"
        if star_end > len(ss):
            return "FAIL_STRUCTURE_STAR_OUT_OF_FOLD_REGION"

    if star_start <= mature_local_pos[0]:
        if mature_local_pos[0] - star_end<3:
            return "FAIL_STRUCTURE_MATURE_STAR_OVERLAP"
        if star_start < 0:
            return "FAIL_STRUCTURE_STAR_OUT_OF_FOLD_REGION"

    mend = ss.rfind(mature_sym, mature_local_pos[0], mature_local_pos[1]-2)
    sstart = dict_bp[mend]
    send = dict_bp[firstbp]
    mature_duplex = ss[mature_local_pos[0]: mend+1]
    star_duplex = ss[sstart:send+1]

    total_dots = mature_duplex.count(".") + star_duplex.count(".")
    total_bps = len(mature_duplex) - mature_duplex.count(".")
    if total_bps<14:
        return "FAIL_STRUCTURE_MATURE_MATCH_SMALL_THAN_14"

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
        return "FAIL_STRUCTURE_STAR_NOT_IN_ONE_ARM"
    loops, bulges = stat_duplex(mature_duplex, star_duplex)
    pass_stat = pass_stat_duplex(loops, bulges)
    if not pass_stat[1]:
        return pass_stat[0]
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


def gen_mapinfo_each_sample(combinedsortedbamname, samplenames, chromID, start, end):
    # start and end are one based.
    dict_map_info = {}
    for sample in samplenames:
        dict_map_info[sample] = {"+" : {},
                                 "-" : {}
        }
    #SAM file and bam file (if accessed using samtools) are 1 based.
    alignments = samtools_view_region(combinedsortedbamname, chromID, start,
                                      end-1)
    for aln in alignments:
        sp = aln.split()
        startpos =0
        readlen = 0
        try:  # skip unmapped alignment
            startpos = int(sp[3])
            readlen = len(sp[9])
        except ValueError:
            continue
        if startpos < start or startpos+readlen > end:
            continue
        depth = int(get_read_depth_fromID_as_string(sp[0]))
        read_sample = "_".join(sp[0].split("_")[:-2]).strip(">")  # sample name could contain "_"
        aln_strand = "+"
        if int(sp[1]) & 16:  # on minus strand
            aln_strand = "-"
        if startpos not in dict_map_info[read_sample][aln_strand]:
            dict_map_info[read_sample][aln_strand][startpos] = []
        dict_map_info[read_sample][aln_strand][startpos].append((sp[9], depth))
    return dict_map_info


#def check_expression_new(dict_map_info, combinedsortedbamname, samplenames, chromID, start,
#                         end, maturepos_genome, mature_depth,
#                         starpos_genome, strand, allow_3nt_overhang):
def check_expression_new(dict_map_info, samplenames, start,
                         end, maturepos_genome, mature_depth,
                         starpos_genome, strand, allow_3nt_overhang):
    # position in dict_map_info are 1 based. All others are zero based.
    dict_info = {}
    mature_len = maturepos_genome[1] - maturepos_genome[0]
    star_len = starpos_genome[1] - starpos_genome[0]
    pre_len = end - start
    dict_info['samplenames'] = samplenames
    for sample in samplenames:
        dict_info[sample] = {}
        dict_info[sample]['reads_pre'] = 0  # number of reads mapped to the precursor region (on the same strand as the locus)
        dict_info[sample]['reads_mature'] = 0  # number of reads mapped exactly to the mature seq (on the same strand as the locus)
        dict_info[sample]['reads_star'] = 0  # number of reads mapped exactly to the star seq (on the same strand as the locus)
        dict_info[sample]['reads_antisense'] = 0  # number of reads mapped to the antisense of the locus.
        dict_info[sample]['reads_maps'] = {}  # reads starts from each position. key: position, value: a list of (readseq, count)
        dict_info[sample]['reads_mature_isoform'] = 0  # possible isoform
        dict_info[sample]['reads_mature_inside'] = 0
        dict_info[sample]['bases_with_reads_start'] = 0
        #dict_info[sample]['bases_with_read_coverage'] = [0] * precursor
        # the following entries are for imperfect star sequences. The three
        # values are for star sequence at [star_start, star_end+1],
        # [star_start+1, star_end+1], and [star_start+1, star_end]
        if allow_3nt_overhang:
            dict_info[sample]['reads_star_imperfect'] = [0, 0, 0]

    for aln_strand in ['+', '-']:
        for startpos in range(start, end):
            for read_sample in dict_map_info:
                if (startpos) in dict_map_info[read_sample][aln_strand]:
                    if strand == aln_strand:
                        dict_info[sample]['bases_with_reads_start'] += 1
                    for read, depth in dict_map_info[read_sample][aln_strand][startpos]:
                        if strand != aln_strand:
                            dict_info[read_sample]['reads_antisense'] = dict_info[read_sample]['reads_antisense'] + depth
                            continue
                        else:
                            dict_info[read_sample]['reads_pre'] = dict_info[read_sample]['reads_pre'] + depth

                        if startpos in dict_info[read_sample]['reads_maps']:
                            dict_info[read_sample]['reads_maps'][startpos].append((read,depth))
                        else:
                            dict_info[read_sample]['reads_maps'][startpos]= [(read,depth)]
                        readlen = len(read)
                        if (startpos==maturepos_genome[0]) and (readlen==(mature_len)):
                            dict_info[read_sample]['reads_mature'] = dict_info[read_sample]['reads_mature'] + depth
                        if (startpos==starpos_genome[0]) and (readlen==star_len):
                            dict_info[read_sample]['reads_star'] = dict_info[read_sample]['reads_star'] + depth
                        if (abs(startpos-maturepos_genome[0]) <=3) and (abs(readlen-mature_len) <=3):
                            dict_info[read_sample]['reads_mature_isoform'] = dict_info[read_sample]['reads_mature_isoform'] + depth
                        if (startpos>=maturepos_genome[0]) and (startpos+readlen <= maturepos_genome[1]):
                            dict_info[read_sample]['reads_mature_inside'] = dict_info[read_sample]['reads_mature_inside'] + depth

                        if allow_3nt_overhang:
                            if (startpos==starpos_genome[0]) and (readlen==star_len+1): # 2, 3
                                dict_info[read_sample]['reads_star_imperfect'][0] = dict_info[read_sample]['reads_star_imperfect'][0] + depth
                            elif (startpos-1==starpos_genome[0]) and (readlen==star_len-1): # 3, 2
                                dict_info[read_sample]['reads_star_imperfect'][1] = dict_info[read_sample]['reads_star_imperfect'][1] + depth
                            elif (startpos-1==starpos_genome[0]) and (readlen==star_len): # 3, 3
                                dict_info[read_sample]['reads_star_imperfect'][2] = dict_info[read_sample]['reads_star_imperfect'][2] + depth
    total_depth_just_this_strand = 0
    total_depth_anti = 0
    total_depth_mature = 0
    total_depth_isoform = 0
    total_depth_star = 0
    mature_depth_each_sample = []
    total_depth_imperfect_star = [0, 0, 0]
    for sample in samplenames:
        total_depth_just_this_strand += dict_info[sample]['reads_pre']
        total_depth_anti += dict_info[sample]['reads_antisense']
        total_depth_mature += dict_info[sample]['reads_mature']
        total_depth_isoform += dict_info[sample]['reads_mature_isoform']
        mature_depth_each_sample.append(dict_info[sample]['reads_mature'])
        total_depth_star += dict_info[sample]['reads_star']
        if allow_3nt_overhang:
            for i in range(3):
                total_depth_imperfect_star[i] +=dict_info[sample]['reads_star_imperfect'][i]
        dict_info[sample]['ratio_bases_with_reads_start'] = float(dict_info[sample]['bases_with_reads_start']) / pre_len

    dict_info['total_depth_just_this_strand'] = total_depth_just_this_strand
    dict_info['total_depth_anti'] = total_depth_anti
    dict_info['total_depth_mature'] = total_depth_mature
    dict_info['total_depth_isoform'] = total_depth_isoform
    dict_info['total_depth_star'] = total_depth_star
    dict_info['total_depth_imperfect_star'] = total_depth_imperfect_star
    dict_info['mature_depth_each_sample'] = mature_depth_each_sample
    if starpos_genome[0] > maturepos_genome[1]:
        dict_info['mature_star_distance'] = starpos_genome[0] - maturepos_genome[1]
    else:
        dict_info['mature_star_distance'] = maturepos_genome[0] - starpos_genome[1]

    def get_imperfect_star_pos():
        max_imperfect_star = max(total_depth_imperfect_star)
        if max_imperfect_star == 0:
            return (-1, 0, 0, 0)
        which = 0
        for i in range(len(total_depth_imperfect_star)):
            if max_imperfect_star == total_depth_imperfect_star[i]:
                which = i
                break
        if which == 0:
            return (which, max_imperfect_star, starpos_genome[0], starpos_genome[1]+1)
        if which == 1:
            return (which, max_imperfect_star, starpos_genome[0]+1, starpos_genome[1])
        if which == 2:
            return (which, max_imperfect_star, starpos_genome[0]+1, starpos_genome[1]+1)


    max_imperfect_star = 0
    if total_depth_star == 0 and allow_3nt_overhang:
        which, max_imperfect_star, star_start, star_end = get_imperfect_star_pos()
        if which == -1:
            dict_info['imperfect_star_start'] = star_start
            dict_info['max_imperfect_star'] = max_imperfect_star
        else:
            dict_info['imperfect_star_start'] = star_start
            dict_info['imperfect_star_end'] = star_end
            dict_info['max_imperfect_star'] = max_imperfect_star
            dict_info['imperfect_star_which'] = which
    mature_star_ratio_tatal =  float(total_depth_mature + max(total_depth_star, max_imperfect_star)) / total_depth_just_this_strand
    dict_info['mature_star_ratio_total'] = mature_star_ratio_tatal
    mature_star_ratio_tatal_both_strand =  float(total_depth_mature + max(total_depth_star, max_imperfect_star)) / (total_depth_just_this_strand + total_depth_anti)
    dict_info['mature_star_ratio_total_both_strand'] = mature_star_ratio_tatal_both_strand
    dict_info['mature_iso_star_ratio_total'] = float(total_depth_isoform + max(total_depth_star, max_imperfect_star)) / total_depth_just_this_strand
    if 'max_imperfect_star' in dict_info:
        dict_info['total_depth_star'] = max( total_depth_star, max_imperfect_star)
    return dict_info


def check_expression(ss, dict_extendregion_info, maturepos_genome, mature_depth, starpos_genome, strand):
    '''
    starpos_genome is the genome coordinate of the star sequence position.

    The format of starpos_genome is [start, end)
    The return value is the depth of the star sequence.
    '''
    starlen = starpos_genome[1] - starpos_genome[0]
    maturelen = maturepos_genome[1] - maturepos_genome[0]
    total_depth = 0
    star_depth = 0
    mature_isoform_depth = 0
    all_depth = []
    for pos in dict_extendregion_info[0]: # for each position
        for s in dict_extendregion_info[0][pos]: # for both strands
            total_depth += dict_extendregion_info[0][pos][s][0][-1]  # total depth on both strand?
            all_depth.append(dict_extendregion_info[0][pos][s][0][-1])
    if starpos_genome[0] not in dict_extendregion_info[0]:
        star_depth = 0
    elif strand not in dict_extendregion_info[0][starpos_genome[0]]:
        star_depth = 0

    #TODO  Maybe better to store length as key depth as value, easier for lookup.
    else:
        for length, depth in  dict_extendregion_info[0][starpos_genome[0]][strand][1]:
            if length == starlen:
                star_depth = depth
    for i in range (maturepos_genome[0]-3, maturepos_genome[0]+4):  # total isoform depth
        if i in dict_extendregion_info[0]:
            if strand not in dict_extendregion_info[0][i]:
                continue
            for length, depth in  dict_extendregion_info[0][i][strand][1]:
                if abs(length-maturelen) < 4:
                    mature_isoform_depth += depth
    ratio = (mature_depth+star_depth)/float(total_depth)
    ratio_matureisoform = (mature_isoform_depth)/float(total_depth)
    #print("T, M, S, R", total_depth, mature_depth, star_depth, ratio)
    return star_depth, ratio, total_depth, ratio_matureisoform, all_depth


def check_loci(structures, matures, region, dict_mapinfo_region, which,
               samplenames, allow_3nt_overhang, allow_no_star, min_mature_len,
               max_mature_len, peak_depth):
    # This function returns either the list 'miRNAs' or the
    # dict 'dict_why_not_miRNA_reasons'
    miRNAs = []
    # record the reason why a region is not predicted as a miRNA.
    # Only the first reason is recorded. For example, it's possible that
    # a region do not have stem-loop structure and not good expression pattern,
    # in this case only NO_STEM_LOOP is recorded because we check this first.
    dict_why_not_miRNA_reasons = {}
    #import pdb; pdb.set_trace();
    key = tuple(region)
    dict_why_not_miRNA_reasons[key] = {}
    dict_why_not_miRNA_reasons[key]['PEAK_PASS_DEPTH'] = "PASSED"
    # no stem loop structure found for this region
    if not structures[1]:
        dict_why_not_miRNA_reasons[key]['HAS_STEMLOOP_STRUCTURE'] = 'FAILED'
        return dict_why_not_miRNA_reasons
    else:
        dict_why_not_miRNA_reasons[key]['HAS_STEMLOOP_STRUCTURE'] = 'PASSED'

    num_matures = len(matures)
    num_matures_out_region = 0
    for m0, m1, strand, mdepth in matures:
        if m1-m0 < min_mature_len or m1-m0 > max_mature_len:
            num_matures_out_region += 1
    # sizes of all matures are out of bounds
    if num_matures_out_region == num_matures:
        dict_why_not_miRNA_reasons[key]['HAS_MATURE_SIZE_IN_RANGE'] = 'FAILED'
        return dict_why_not_miRNA_reasons
    else:
        dict_why_not_miRNA_reasons[key]['HAS_MATURE_SIZE_IN_RANGE'] = 'PASSED'


    for m0, m1, strand, mdepth in sorted(matures, key=lambda m: m[3], reverse=True):
        # TODO Move this if to previous stage.
        # if mature length is not in [18-24], then ignore this
        if m1-m0 < min_mature_len or m1-m0 > max_mature_len:
            continue
        lowest_energy = 0
        lowest_energy = 0
        outputinfo = []

        for energy, foldstart, ss, sstype in structures[1]:
            if energy > lowest_energy:
                continue
            else:
                #key = (tuple(region), m0, m1, strand, mdepth, ss, sstype)
                key1 = (m0, m1, strand, ss)
                dict_why_not_miRNA_reasons[key][key1] = {}
                structinfo = get_maturestar_info(ss, (m0,m1), foldstart,
                                                 foldstart+len(ss),
                                                 region[1][0],
                                                 region[1][1],
                                                 strand)
                if isinstance(structinfo, str): # failed
                    dict_why_not_miRNA_reasons[key][key1][structinfo] = "FAILED"
                else:
                    #exprinfo = check_expression(ss,dict_aln, (m0,m1),
                    #                            mdepth,
                    #                            (structinfo[0],structinfo[1]),
                    #                            strand)
                    exprinfo = check_expression_new(dict_mapinfo_region,
                                                    samplenames,
                                                    structinfo[2],
                                                    structinfo[3],
                                                    (m0, m1),
                                                    mdepth,
                                                    (structinfo[0], structinfo[1]),
                                                    strand,
                                                    allow_3nt_overhang
                                                 )
                    if exprinfo['mature_star_distance'] <= 4:
                           dict_why_not_miRNA_reasons[key][key1]["FAIL_EXPRESS_PATTERN_MATURE_STAR_TOO_CLOSE"] = "FAILED"
                           dict_why_not_miRNA_reasons[key][key1]['ss_info'] = [region[0], structinfo[2],
                                                                               structinfo[3], m0, m1, structinfo[0],
                                                                               structinfo[1], ss, strand, True, exprinfo]
                           continue
                    if exprinfo['total_depth_star'] > 0:  # has star expression
                        if exprinfo['mature_star_ratio_total'] < 0.2:
                            dict_why_not_miRNA_reasons[key][key1]["FAIL_EXPRESS_PATTERN_HAS_STAR_BUT_TOO_FEW_READS_MAPPED_TO_DUPLEX"] = "FAILED"
                            dict_why_not_miRNA_reasons[key][key1]['ss_info'] = [region[0], structinfo[2],
                                                                                structinfo[3], m0, m1, structinfo[0],
                                                                                structinfo[1], ss, strand, True, exprinfo]
                            continue
                        else:
                            #  The last 'True' means this is an confident miRNA
                            # regionID, fold start, fold end, mature start, mature end, star start,
                            # star end, ss, strand, has_star
                            outputinfo = [region[0], structinfo[2],
                                          structinfo[3], m0, m1, structinfo[0],
                                          structinfo[1], ss, strand, True, exprinfo]
                            if 'max_imperfect_star' in exprinfo:
                                outputinfo[5] = exprinfo['imperfect_star_start']
                                outputinfo[6] = exprinfo['imperfect_star_end']
                            lowest_energy = energy
                    else:  #  no star expression
                        if not allow_no_star:
                            dict_why_not_miRNA_reasons[key][key1]['FAIL_EXPRESS_PATTERN_NO_STAR_EXPRESSION_DISALLOW_NO_STAR'] = "FAILED"
                            dict_why_not_miRNA_reasons[key][key1]['ss_info'] = [region[0], structinfo[2],
                                                                                structinfo[3], m0, m1, structinfo[0],
                                                                                structinfo[1], ss, strand, False, exprinfo]
                            dict_why_not_miRNA_reasons[key][key1]['expression_info'] = exprinfo
                            continue
                        too_many_start = False
                        for sample in samplenames:
                            if exprinfo[sample]['ratio_bases_with_reads_start'] > 0.5:
                                too_many_start = True
                        if too_many_start:
                            dict_why_not_miRNA_reasons[key][key1]['FAIL_EXPRESS_PATTERN_NO_STAR_EXPRESSION_TOO_MANY_START'] = "FAILED"
                            dict_why_not_miRNA_reasons[key][key1]['ss_info'] = [region[0], structinfo[2],
                                                                                structinfo[3], m0, m1, structinfo[0],
                                                                                structinfo[1], ss, strand, False, exprinfo]
                            dict_why_not_miRNA_reasons[key][key1]['expression_info'] = exprinfo
                            continue
                        expressed_sample = len([x for x in exprinfo['mature_depth_each_sample'] if x > 0])
                        expressed_all = False
                        if expressed_sample == len(exprinfo['mature_depth_each_sample']):
                            expressed_all = True
                        if exprinfo['mature_iso_star_ratio_total'] >= 0.8 and (expressed_all or exprinfo['total_depth_mature'] >= 1000):
                            outputinfo = [region[0], structinfo[2],
                                          structinfo[3], m0, m1, structinfo[0],
                                          structinfo[1], ss, strand, False, exprinfo]
                            lowest_energy = energy
                        else:
                            if exprinfo['mature_iso_star_ratio_total'] < 0.8:
                                dict_why_not_miRNA_reasons[key][key1]['FAIL_EXPRESS_PATTERN_NO_STAR_MATURE_STAR_RATIO_TOO_SMALL'] = "FAILED"
                            if exprinfo['total_depth_mature'] <= 100:
                                dict_why_not_miRNA_reasons[key][key1]['FAIL_EXPRESS_PATTERN_NO_STAR_MATURE_DEPTH_TOO_SMALL'] = "FAILED"
                            if not expressed_all:
                                dict_why_not_miRNA_reasons[key][key1]['FAIL_EXPRESS_PATTERN_NO_STAR_MATURE_NOT_IN_ALL_SAMPLE'] = "FAILED"
                            dict_why_not_miRNA_reasons[key][key1]['ss_info'] = [region[0], structinfo[2],
                                                                                structinfo[3], m0, m1, structinfo[0],
                                                                                structinfo[1], ss, strand, True, exprinfo]
                            dict_why_not_miRNA_reasons[key][key1]['expression_info'] = exprinfo
        if outputinfo:  # the loci contains an miRNA
            miRNAs.append(outputinfo)
    if miRNAs:
        return miRNAs
    else:
        return dict_why_not_miRNA_reasons


def filter_next_loci(alndumpname, rnalfoldoutname, combinedsortedbamname,
                     samplenames, allow_3nt_overhang, allow_no_star,
                     output_details, min_mature_len, max_mature_len, peak_depth,
                     minlen=50):
    list_miRNA_loci = []
    alnf = open(alndumpname)
    # (which, peak region, [energy, fold_start, ss, sstype])
    ss_generator = get_structures_next_extendregion(rnalfoldoutname, minlen)
    while True:
        region = None
        which = None
        dict_aln = None
        matures = None
        try:
            region, which, dict_aln, matures = cPickle.load(alnf)
        except EOFError:
            raise StopIteration
        # position in dict_mapinfo_region, region are all one-based.
        dict_mapinfo_region = gen_mapinfo_each_sample(combinedsortedbamname,
                                                      samplenames, region[0],
                                                      region[1][0], region[1][1]
                                                  )
        structures = []
        if which == "0": #only one extend region for this loci
            structures = next(ss_generator)
            peak = structures[1]
            which_and_ss = (structures[0], structures[2])
            # all positions are one based.
            miRNAs = check_loci(which_and_ss, matures, region,
                                dict_mapinfo_region, which, samplenames,
                                allow_3nt_overhang,
                                allow_no_star,
                                min_mature_len,
                                max_mature_len,
                                peak_depth)
            if isinstance(miRNAs, dict): # not a miRNA
                if not output_details:
                    continue
                miRNAs['which'] = '0'
                miRNAs['peak'] = peak
            yield miRNAs
            continue

        else:
            structuresL = next(ss_generator)  # L
            structuresR = next(ss_generator)  # R
            peak = structuresL[1]
            which_and_ssL = (structuresL[0], structuresL[2])
            which_and_ssR = (structuresR[0], structuresR[2])
            region1, which1, dict_aln1, matures1 = cPickle.load(alnf)
            dict_mapinfo_region1 = gen_mapinfo_each_sample(combinedsortedbamname,
                                                           samplenames, region1[0],
                                                           region1[1][0], region1[1][1]
                                                      )
            miRNAsL = check_loci(which_and_ssL, matures, region,
                                 dict_mapinfo_region, which,
                                 samplenames, allow_3nt_overhang,
                                 allow_no_star,
                                 min_mature_len,
                                 max_mature_len,
                                 peak_depth)
            if isinstance(miRNAsL, dict):
                if not output_details:
                    pass
                miRNAsL['which'] = 'L'
                miRNAsL['peak'] = peak
            yield miRNAsL
            if isinstance(miRNAsL, dict):
                miRNAsR = check_loci(which_and_ssR, matures1, region1,
                                     dict_mapinfo_region1, which1,
                                     samplenames,allow_3nt_overhang,
                                     allow_no_star,
                                     min_mature_len,
                                     max_mature_len,
                                     peak_depth
                                 )
                if isinstance(miRNAsR, dict):
                    if not output_details:
                        continue
                    miRNAsR['which'] = 'R'
                    miRNAsR['peak'] = peak
                yield miRNAsR
                continue


def gen_miRNA_loci_nopredict(alndumpnames, rnalfoldoutnames, combinedsortedbamname, samplenames,
                             allow_3nt_overhang, allow_no_star, output_details,
                             min_mature_len, max_mature_len, minlen, peak_depth,
                             logger):

    def gen_miRNA_loci_local(queue, alndumpname, rnalfoldoutname, minlen,
                             samplenames, combinedsortedbamname,
                             allow_3nt_overhang, allow_no_star,
                             output_details, min_mature_len,
                             max_mature_len, peak_depth):
        mir_generator = filter_next_loci( alndumpname, rnalfoldoutname,
                                          combinedsortedbamname,
                                          samplenames, allow_3nt_overhang,
                                          allow_no_star,
                                          output_details,
                                          min_mature_len,
                                          max_mature_len,
                                          peak_depth,
                                          minlen=minlen)

        mirnas = []
        for mir in mir_generator:
            if isinstance(mir, list):
                mirnas.append((mir, True))
            else:
                mirnas.append((mir, False))
        queue.put(mirnas)
        queue.put("done")
        queue.close()
    miRNAqueue = multiprocessing.Queue()
    jobs = []
    finalresult = []
    failed_reasons = []
    for i in range(len(alndumpnames)):
        p = multiprocessing.Process(target = gen_miRNA_loci_local, args=(
            miRNAqueue, alndumpnames[i], rnalfoldoutnames[i],minlen,samplenames,
            combinedsortedbamname,
            allow_3nt_overhang,
            allow_no_star,
            output_details,
            min_mature_len,
            max_mature_len,
            peak_depth
        ))
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
                    if mirna[1]:
                        finalresult.append(mirna[0][0])
                    else:
                        failed_reasons.append(mirna[0])
        except Queue.Empty:
            time.sleep(2)
            continue


    return finalresult, failed_reasons


def convert_failure_reasons_list(list_failure_reasons):
    out_dict_reasons = {} # key is chromosome
    for dict_reason in list_failure_reasons:
        which = '0'
        peak = ''
        chrom = ''
        strand = ''
        region = ''
        for key in dict_reason: # key is a tuple: ('Chr1', (21784879, 21785192), '-')
            if key == 'which':
                which = dict_reason[key]
                continue
            if key == 'peak':
                peak = dict_reason[key]
                continue
            chrom = key[0]
            if chrom not in out_dict_reasons:
                out_dict_reasons[chrom] = {'+': {},
                                           '-': {}}
            strand = key[2]
            region = key[1]
            out_dict_reasons[chrom][strand][key[1]] = dict_reason[key]
        out_dict_reasons[chrom][strand][region]['which'] = which
        out_dict_reasons[chrom][strand][region]['peak'] = peak
    return out_dict_reasons


def write_dict_reasons(dict_reasons, fname, faileddumpname, fastaname, combinedsortedbamname, samplenames, outfolder):
    list_not_pass_expression = []
    dict_express_info_all = {}
    outf = open(fname, 'w')
    for chrom in dict_reasons:  # key is chromosome
        for strand in dict_reasons[chrom]:
            for region in dict_reasons[chrom][strand]:
                outf.write("===========================================================\n")
                regionkey1 = chrom + ":" + str(region[0]) + "-" + str(region[1])
                regionkey = regionkey1 + "\t" + strand + "\tpeak-region:" + dict_reasons[chrom][strand][region]['peak'] + "\t" + dict_reasons[chrom][strand][region]['which']
                outf.write( regionkey +"\n")
                dict_region = dict_reasons[chrom][strand][region]
                for key in dict_region:
                    if isinstance(key, str):
                        outf.write(key + ":" + dict_region[key] + "\n")
                for key in dict_region:
                    if isinstance(key, tuple):
                        outf.write("MATURE region: {0}-{1}, SS: {2}".format(key[0], key[1], key[3]) + "\n")
                        for v in dict_region[key]:
                            if isinstance(dict_region[key][v], str):
                                outf.write(v + "\t" + dict_region[key][v] + "\n")
                        if 'expression_info' in dict_region[key]:
                            dict_info = dict_region[key]['expression_info']
                            dict_express_info_all[regionkey1] = dict_info
                            #outf.write("total_depth: {0}, mature_depth: {1}, star_depth:{2}, mature_star_ratio:{3}, ".format(dict_info['total_depth_just_this_strand'], dict_info['total_depth_mature'], dict_info['total_depth_star'], dict_info['mature_star_ratio_total']) + "\n")
                            for k in dict_info:
                                if isinstance(dict_info[k], dict):
                                    continue
                                outf.write(k + "\t" + str(dict_info[k]) + "\n")
                        if 'ss_info' in dict_region[key]:
                            list_not_pass_expression.append(dict_region[key]['ss_info'])
                outf.write("\n")
    outdump = open(faileddumpname, 'w')
    cPickle.dump(dict_express_info_all, outdump)
    dict_not_mirna_info = gen_mirna_info(list_not_pass_expression, fastaname, combinedsortedbamname, samplenames)
    gen_map_result(dict_not_mirna_info, outfolder)


def get_mature_stemloop_star_seq(seqid, mstart, mend, start, end, sstart, send, fastaname):
    region1 = seqid+":"+str(mstart)+"-"+str(mend)
    region2 = seqid+":"+str(start)+"-"+str(end)
    region3 = seqid+":"+str(sstart)+"-"+str(send)
    mature_command = "samtools faidx " + fastaname + " " +region1
    stemloop_command = "samtools faidx " + fastaname + " " +region2
    star_command = "samtools faidx " + fastaname + " " +region3
    try:
        samtools_process = subprocess.Popen(mature_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
        mature_seq = str("".join(outmessage.decode().split("\n")[1:]))
        samtools_process = subprocess.Popen(stemloop_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
        stemloop_seq = str("".join(outmessage.decode().split("\n")[1:]))
        samtools_process = subprocess.Popen(star_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmessage, outerr = samtools_process.communicate()
        star_seq = str("".join(outmessage.decode().split("\n")[1:]))
    except Exception as e:
        sys.stderr.write("Error occurred when calling samtools in get_mature_stemloop_star_seq.\n")
        sys.exit(-1)
    return region1, mature_seq, region2, stemloop_seq, region3, star_seq


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
        sys.stderr.write("Error occurred when calling samtools in get_mature_stemloop_seq.\n")
        sys.exit(-1)
    return region1, mature_seq, region2, stemloop_seq


def adjust_mature_star(resultlist):
    for idx, mirna in enumerate(resultlist):
        if mirna[-1]['total_depth_mature'] < mirna[-1]['total_depth_star']:
            mirna[-1]['total_depth_mature'], mirna[-1]['total_depth_star'] = mirna[-1]['total_depth_star'], mirna[-1]['total_depth_mature']
            mirna[3], mirna[5] = mirna[5], mirna[3]
            mirna[4], mirna[6] = mirna[6], mirna[4]
            mirna[-1]['switch'] = True

def gen_gff_from_result(resultlist, gffname):

    f = open(gffname, 'w')
    resultlist.sort()
    for idx, mirna in enumerate(resultlist):
        maturename = "miRNA_"+ str(idx)
        mirname = "miRNA-precursor_"+ str(idx)
        has_star_expression = "y"
        overhangsize = "2:2"
        if mirna[-1]['total_depth_star'] == 0:
            has_star_expression = 'n'

        if 'max_imperfect_star' in mirna[-1] and mirna[-1]['max_imperfect_star'] != 0:
            if mirna[-1]['imperfect_star_which'] == 2:
                overhangsize = "3:3"
            else:
                overhangsize = "2:3"
        other = 'mature_expressed=y;star_expressed=' + has_star_expression + ";" + "overhangsize=" + overhangsize
        write_gff_line(mirna[0],mirna[1],mirna[2]-1,mirna[-3],mirname, mirname,
                       feature="miRNA-precursor",other=other, fout=f)
        write_gff_line(mirna[0],mirna[3],mirna[4]-1,mirna[-3],maturename,
                       maturename, feature="miRNA",fout=f)
    f.close()


def gen_mirna_info(resultlist, fastaname, combinedsortedbamname, samplenames):
    '''
    Each result:
    Seqid, locus_start, locus_end, mature_start, mature_end, star_start
    star_end, structure, strand, star_present
    '''

    dict_mirna_info = {'idlist':[],  # list of all miRNA ids
                       'info_dict':{},  # dict of miRNA infomations
                       'samples':samplenames
    }
    for idx, mirna in enumerate(resultlist):
        maturename = "miRNA_" + str(idx)
        mirname = "miRNA-precursor_" + str(idx)
        dict_mirna_info['idlist'].append(mirname)
        seqs = get_mature_stemloop_star_seq(mirna[0],mirna[3], mirna[4]-1, mirna[1], mirna[2]-1, mirna[5], mirna[6]-1, fastaname)
        matureseq = seqs[1].upper().replace("T", "U")
        stemloopseq = seqs[3].upper().replace("T", "U")
        starseq = seqs[5].upper().replace("T", "U")
        structure = mirna[-4]  # structure is already reversed if on minus
        strand = mirna[8]
        maturestart_inseq = mirna[3] - mirna[1]
        matureend_inseq = mirna[4] - mirna[1]
        starstart_inseq = mirna[5] - mirna[1]
        starend_inseq = mirna[6] - mirna[1]

        #if mirna[-2] == "-":
        #    matureseq = get_reverse_complement(matureseq)
        #    stemloopseq = get_reverse_complement(stemloopseq)
        #    starseq = get_reverse_complement(starseq)
        dict_mirna_info['info_dict'][mirname] = {}
        dict_mirna_info['info_dict'][mirname]['precursor'] = stemloopseq
        dict_mirna_info['info_dict'][mirname]['matureseq'] = matureseq
        dict_mirna_info['info_dict'][mirname]['starseq'] = starseq
        dict_mirna_info['info_dict'][mirname]['starseq'] = starseq
        dict_mirna_info['info_dict'][mirname]['strand'] = mirna[8]
        dict_mirna_info['info_dict'][mirname]['precursor_id'] = mirname
        dict_mirna_info['info_dict'][mirname]['mature_id'] = maturename
        dict_mirna_info['info_dict'][mirname]['ss'] = structure
        dict_mirna_info['info_dict'][mirname]['chr'] = mirna[0]
        dict_mirna_info['info_dict'][mirname]['locus_start'] = mirna[1]
        dict_mirna_info['info_dict'][mirname]['locus_end'] = mirna[2]  # [start, end)
        dict_mirna_info['info_dict'][mirname]['mature_start'] = mirna[3]  # [start, end)
        dict_mirna_info['info_dict'][mirname]['mature_end'] = mirna[4]  # [start, end)
        dict_mirna_info['info_dict'][mirname]['star_start'] = mirna[5]  # [start, end)
        dict_mirna_info['info_dict'][mirname]['star_end'] = mirna[6]  # [start, end)
        for sample in samplenames:
            dict_mirna_info['info_dict'][mirname][sample] = {}
            dict_mirna_info['info_dict'][mirname][sample]['reads_pre'] = 0  # number of reads mapped to the precursor region (on the same strand as the locus)
            dict_mirna_info['info_dict'][mirname][sample]['reads_mature'] = 0  # number of reads mapped exactly to the mature seq (on the same strand as the locus)
            dict_mirna_info['info_dict'][mirname][sample]['reads_star'] = 0  # number of reads mapped exactly to the star seq (on the same strand as the locus)
            dict_mirna_info['info_dict'][mirname][sample]['reads_antisense'] = 0  # number of reads mapped to the antisense of the locus.
            dict_mirna_info['info_dict'][mirname][sample]['reads_maps'] = {}  # reads starts from each position. key: position, value: a list of (readseq, count)
        alignments = samtools_view_region(combinedsortedbamname, mirna[0],
                                          mirna[1], mirna[2]-1)
        # import pdb; pdb.set_trace();
        for aln in alignments:
            sp = aln.split()
            startpos =0
            readlen = 0
            try:  # skip unmapped alignment
                startpos = int(sp[3])
                readlen = len(sp[9])
            except ValueError:
                continue
            if startpos < mirna[1] or startpos+readlen > mirna[2]:
                continue
            depth = int(get_read_depth_fromID_as_string(sp[0]))
            read_sample = "_".join(sp[0].split("_")[:-2]).strip(">")  # sample name could contain "_"
            aln_strand = "+"
            if int(sp[1]) & 16:  # on minus strand
                aln_strand = "-"
            if aln_strand != strand:  # antisense strand
                dict_mirna_info['info_dict'][mirname][read_sample]['reads_antisense'] = dict_mirna_info['info_dict'][mirname][read_sample]['reads_antisense'] + depth
                continue
            dict_mirna_info['info_dict'][mirname][read_sample]['reads_pre'] = dict_mirna_info['info_dict'][mirname][read_sample]['reads_pre'] + depth
            if startpos in dict_mirna_info['info_dict'][mirname][read_sample]['reads_maps']:
                dict_mirna_info['info_dict'][mirname][read_sample]['reads_maps'][startpos].append((sp[9],depth))
            else:
                dict_mirna_info['info_dict'][mirname][read_sample]['reads_maps'][startpos]= [(sp[9],depth)]
            if (startpos==mirna[3]) and (readlen==len(matureseq)):
                dict_mirna_info['info_dict'][mirname][read_sample]['reads_mature'] = dict_mirna_info['info_dict'][mirname][read_sample]['reads_mature'] + depth
            if (startpos==mirna[5]) and (readlen==len(starseq)):
                dict_mirna_info['info_dict'][mirname][read_sample]['reads_star'] = dict_mirna_info['info_dict'][mirname][read_sample]['reads_star'] + depth
    return dict_mirna_info



def gen_miRNA_stat(dict_mirna_info):
    dict_length_distrubution = {}
    first_char_distribution = {}
    total_prediction = len(dict_mirna_info['idlist'])
    for mirna in dict_mirna_info['info_dict']:
        mirnalen = len(dict_mirna_info['info_dict'][mirna]['matureseq'])
        dict_length_distrubution[mirnalen] = dict_length_distrubution.get(mirnalen, 0) + 1
        first_char_distribution[dict_mirna_info['info_dict'][mirna]['matureseq'][0]] \
            = first_char_distribution.get(dict_mirna_info['info_dict'][mirna]['matureseq'][0], 0) +1
    return total_prediction, dict_length_distrubution, first_char_distribution


def gen_csv_table(dict_mirna_info, csvname):
    f = open(csvname, 'w')
    header1 = "miRNAID, Seqid(chromosome), start position, end position, strand, precursor sequence, secondary structure, mature sequence, star sequence, "
    header2 = "miRNAID, Seqid(chromosome), start position, end position, strand, precursor sequence, secondary structure, mature sequence, star sequence, "
    for sample in dict_mirna_info['samples']:
        header1 = header1 + sample +"," + sample +","+ sample +","+ sample +","
        header2 = header2 + "reads mapped to precursor, reads mapped to mature, reads mapped to star, reads mapped to antisense region,"
    f.write(header1+"\n")
    f.write(header2+"\n")
    for mirna in dict_mirna_info['idlist']:
        preseq = dict_mirna_info['info_dict'][mirna]['precursor']
        matureseq = dict_mirna_info['info_dict'][mirna]['matureseq']
        starseq = dict_mirna_info['info_dict'][mirna]['starseq']
        if dict_mirna_info['info_dict'][mirna]['strand']=="-":
            preseq = get_reverse_complement(preseq)
            matureseq = get_reverse_complement(matureseq)
            starseq = get_reverse_complement(starseq)
        outstr = [dict_mirna_info['info_dict'][mirna]['precursor_id'],
                  dict_mirna_info['info_dict'][mirna]['chr'],
                  str(dict_mirna_info['info_dict'][mirna]['locus_start']),
                  str(dict_mirna_info['info_dict'][mirna]['locus_end']),
                  dict_mirna_info['info_dict'][mirna]['strand'],
                  preseq,
                  dict_mirna_info['info_dict'][mirna]['ss'],
                  matureseq,
                  starseq]
        for sample in dict_mirna_info['samples']:
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_pre']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_mature']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_star']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_antisense']))
        f.write(", ".join(outstr)+"\n")


def gen_search_miRBase_str(matureseq, taxon):
    formstr = '''<form name="myform" action="http://www.mirbase.org/cgi-bin/blast.pl" method="POST"> \
    <div align="center"> \
    <input type="hidden" size="25" name="sequence" value=" ''' +matureseq + '''"> \
    <input type="hidden" size="25" name="seqfile" value=""> \
    <input type="hidden" size="25" name="type" value="mature"> \
    <input type="hidden" size="25" name="search_method" value="blastn"> \
    <input type="hidden" size="25" name="evalue" value="10"> \
    <input type="hidden" size="25" name="maxalign" value="100"> \
    <input type="hidden" size="25" name="taxon" value=" '''+ taxon + '''"> \
    <input type="submit" value="Search mature on miRBase ( ''' + taxon + ''' )" onclick="this.form.target='_blank';return true;"> \
    </div> \
    </form>'''
    return formstr

def gen_html_table_file(dict_mirna_info, htmlname):
    num_of_samples = len(dict_mirna_info['samples'])
    colors = ['#A9E2F3', '#ACFA58',  '#F5A9BC']
    f = open(htmlname, 'w')

    total_count, dict_len, dict_first = gen_miRNA_stat(dict_mirna_info)
    f.write('<h1 > microRNAs predicted by miR-PREFeR </h1>\n')
    f.write('<div>\n')
    f.write('<h2 > Total number of prediction:' + str(len(dict_mirna_info['idlist']))+ '  </h2>\n')

    f.write('<h3>Distribution of the lengths of the mature sequences</h3>\n')
    f.write('<table border="1">\n')
    f.write("\t<thead>\n")
    f.write("\t\t<tr>\n")
    f.write('\t\t<th>Length</th>\n')
    f.write('\t\t<th>Count</th>\n')
    f.write("\t\t</tr>\n")
    f.write("\t</thead>\n")

    f.write("\t<tbody>\n")
    for k in sorted(dict_len.keys()):
        f.write('\t\t<tr>\n')
        f.write('\t\t\t<td>' + str(k) + ' </td>\n')
        f.write('\t\t\t<td>' + str(dict_len[k]) + ' </td>\n')
        f.write('\t\t</tr>\n')
    f.write("\t</tbody>\n")
    f.write("</table>\n")
    f.write('</div>\n')
    f.write('<h3>Distribution of the nucleotide of the first base of the mature sequences</h3>\n')
    f.write('<table border="1">\n')
    f.write("\t<thead>\n")
    f.write("\t\t<tr>\n")
    f.write('\t\t<th>Nucleotide</th>\n')
    f.write('\t\t<th>Count</th>\n')
    f.write("\t\t</tr>\n")
    f.write("\t</thead>\n")

    f.write("\t<tbody>\n")
    for k in sorted(dict_first.keys()):
        f.write('\t\t<tr>\n')
        f.write('\t\t\t<td>' + str(k) + ' </td>\n')
        f.write('\t\t\t<td>' + str(dict_first[k]) + ' </td>\n')
        f.write('\t\t</tr>\n')
    f.write("\t</tbody>\n")
    f.write("</table>\n")
    f.write('</div>')

    f.write('<div>')
    f.write('<h3>Detailed infomation </h3>\n')
    f.write('<table border="1">\n')

    f.write('<colgroup>\n')
    f.write('\t<col span=8 style="background-color:#CECEF6">\n')
    for i in range(num_of_samples):
        f.write('\t<col span="4" style="background-color:'+colors[i%len(colors)]+'">')
    f.write('</colgroup>\n')
    f.write("\t<thead>\n")
    f.write("\t\t<tr>\n")
    col1 = ['miRNA precursor ID', 'Chromosome', 'start position', 'end position', 'strand', 'precursor sequence and secondary structure','mature sequence', 'star sequence']
    col2 = []
    for c in col1:
        f.write('\t\t\t<th rowspan="2">' + c +  '</th>\n')

    for sample in dict_mirna_info['samples']:
        f.write('\t\t\t<th colspan="4">'+sample+'</th>\n')
    f.write('\t\t\t<th colspan="2">read mappings</th>\n')

    f.write("\t\t</tr>\n")
    f.write("\t\t<tr>\n")
    for c in col2:
        f.write('\t\t\t<th>' + c +  '</th>\n')
    for sample in dict_mirna_info['samples']:
        f.write('\t\t\t<th> reads mapped to precursor </th>\n')
        f.write('\t\t\t<th> reads mapped to mature </th>\n')
        f.write('\t\t\t<th> reads mapped to star </th>\n')
        f.write('\t\t\t<th> reads mapped to antisense region </th>\n')
    f.write("\t\t</tr>\n")
    f.write("\t</thead>\n")
    f.write("\t<tbody>\n")

    for mirna in dict_mirna_info['idlist']:
        f.write("\t\t<tr>\n")
        preseq = dict_mirna_info['info_dict'][mirna]['precursor']
        matureseq = dict_mirna_info['info_dict'][mirna]['matureseq']
        starseq = dict_mirna_info['info_dict'][mirna]['starseq']
        if dict_mirna_info['info_dict'][mirna]['strand']=="-":
            preseq = get_reverse_complement(preseq)
            matureseq = get_reverse_complement(matureseq)
            starseq = get_reverse_complement(starseq)
        outstr = [dict_mirna_info['info_dict'][mirna]['precursor_id'],
                  dict_mirna_info['info_dict'][mirna]['chr'],
                  str(dict_mirna_info['info_dict'][mirna]['locus_start']),
                  str(dict_mirna_info['info_dict'][mirna]['locus_end']),
                  dict_mirna_info['info_dict'][mirna]['strand']]
        for s in outstr:
            f.write('\t\t\t<td nowrap>'+s+'</td>\n')
        f.write('\t\t\t<td nowrap> <code>'+preseq+'<BR>'+dict_mirna_info['info_dict'][mirna]['ss']+ ' </code></td>')
        outstr = [ matureseq + gen_search_miRBase_str(matureseq, "Viridiplantae") + gen_search_miRBase_str(matureseq, "ALL"),
                  starseq]
        for sample in dict_mirna_info['samples']:
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_pre']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_mature']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_star']))
            outstr.append(str(dict_mirna_info['info_dict'][mirna][sample]['reads_antisense']))
        for s in outstr:
            f.write('\t\t\t<td nowrap>'+s+'</td>\n')
        f.write('\t\t\t<td><a href="readmapping/'+mirna+'.map.txt"'+ ' target="_blank">Click to see detailed mapping.</a></td>')
        f.write("\t\t</tr>\n")

    f.write("\t</tbody>\n")
    f.write("</table>\n")
    f.write('</div>\n')


def gen_map_result(dict_mirna_info, foldername):
    for mirna in dict_mirna_info['idlist']:
        #htmlfile = open(mirna+".map.html",'w')
        txtfile = open(os.path.join(foldername,mirna+".map.txt"),'w')
        outlines = []
        idline = ">"+dict_mirna_info['info_dict'][mirna]['precursor_id'] + " " \
                 + dict_mirna_info['info_dict'][mirna]['chr'] + ":" +          \
                  str(dict_mirna_info['info_dict'][mirna]['locus_start']) +    \
                 "-" + str(dict_mirna_info['info_dict'][mirna]['locus_end']) + \
                 " " + dict_mirna_info['info_dict'][mirna]['strand']
        outlines.append(idline)
        preseq = dict_mirna_info['info_dict'][mirna]['precursor']
        matureseq = dict_mirna_info['info_dict'][mirna]['matureseq']
        starseq = dict_mirna_info['info_dict'][mirna]['starseq']
        ss = dict_mirna_info['info_dict'][mirna]['ss']
        if dict_mirna_info['info_dict'][mirna]['strand']=="-":
            preseq = get_reverse_complement(preseq)
            matureseq = get_reverse_complement(matureseq)
            starseq = get_reverse_complement(starseq)

        for sample in dict_mirna_info['samples']:
            outlines.append(">> Read mappings for sample: "+sample)
            outlines.append("5'->3'")
            outlines.append(preseq+"\ttotal_mapped_reads="+str(dict_mirna_info['info_dict'][mirna][sample]['reads_pre']))
            outlines.append(ss)
            sortedkeys = sorted(dict_mirna_info['info_dict'][mirna][sample]['reads_maps'].keys())
            if dict_mirna_info['info_dict'][mirna]['strand'] == "-":
                sortedkeys = sorted(dict_mirna_info['info_dict'][mirna][sample]['reads_maps'].keys(), reverse=True)
            for startpos in sortedkeys:
                for r in sorted(dict_mirna_info['info_dict'][mirna][sample]['reads_maps'][startpos], key=lambda tup: len(tup[0])):
                    outstr = ""
                    padchar = '.'
                    if (startpos==dict_mirna_info['info_dict'][mirna]['mature_start']) and (len(r[0])==len(matureseq)):
                        padchar = 'm'
                    elif (startpos==dict_mirna_info['info_dict'][mirna]['star_start']) and (len(r[0])==len(starseq)):
                        padchar = 's'
                    prepadlen = startpos - dict_mirna_info['info_dict'][mirna]['locus_start']
                    for i in range(prepadlen):
                        outstr = outstr + padchar
                    outstr = outstr + r[0]
                    postpadlen = len(preseq) - len(outstr)
                    for i in range(postpadlen):
                        outstr = outstr + padchar
                    if dict_mirna_info['info_dict'][mirna]['strand'] == '-':
                        outstr = get_reverse_complement(outstr).replace("U", "T")
                    outstr = outstr + "\tdepth=" + str(r[1]) + ", length="+str(len(r[0]))
                    if padchar=='m':
                        outstr = outstr + " [mature]"
                    if padchar=='s':
                        outstr = outstr + " [star]"
                    outlines.append(outstr)
        for line in outlines:
            txtfile.write(line+"\n")



def gen_mirna_fasta_ss_from_result(resultlist, maturename, stemloopname,
                                   fastaname, ssname):
    '''
    Each result:
    Seqid, locus_start, locus_end, mature_start, mature_end, star_start
    star_end, structure, strand, star_present, mapinfo
    '''
    f1 = open(maturename, 'w')
    f2 = open(stemloopname, 'w')
    f3 = open(ssname, 'w')
    for idx, mirna in enumerate(resultlist):
        maturename = "miRNA_" + str(idx)
        mirname = "miRNA-precursor_" + str(idx)
        seqs = get_mature_stemloop_star_seq(mirna[0],mirna[3], mirna[4]-1, mirna[1], mirna[2]-1, mirna[5], mirna[6]-1, fastaname)
        matureid = ">" + seqs[0] + " " + mirna[-3]+ " " + mirname
        stemloopid = ">" + seqs[2] + " " + mirna[-3]+ " " + mirname
        matureseq = seqs[1].upper().replace("T", "U")
        stemloopseq = seqs[3].upper().replace("T", "U")
        structure = mirna[-4]
        maturestart = mirna[3] - mirna[1]
        matureend = mirna[4] - mirna[1]
        starstart = mirna[5] - mirna[1]
        starend = mirna[6] - mirna[1]
        mature_M = "M" * (matureend - maturestart)
        star_S = "S" * (starend - starstart)
        seq_dot = ""
        if maturestart < starstart:
            seq_dot = seq_dot + "." * maturestart
            seq_dot = seq_dot + mature_M
            seq_dot = seq_dot + "." * (starstart-matureend)
            seq_dot = seq_dot + star_S
            seq_dot = seq_dot + "." * (len(stemloopseq)-starend)
        else:
            seq_dot = seq_dot + "." * starstart
            seq_dot = seq_dot + star_S
            seq_dot = seq_dot + "." * (maturestart-starend)
            seq_dot = seq_dot + mature_M
            seq_dot = seq_dot + "." * (len(stemloopseq)-matureend)
        if mirna[-3] == "-":
            matureseq = get_reverse_complement(matureseq)
            stemloopseq = get_reverse_complement(stemloopseq)
            #structure = structure[::-1]
            #structure = structure.replace('(','<')
            #structure = structure.replace(')','(')
            #structure = structure.replace('<',')')
            seq_dot = seq_dot[::-1]
        f1.write(matureid+"\n")
        f1.write(matureseq+"\n")
        f2.write(stemloopid+"\n")
        f2.write(stemloopseq+"\n")
        f3.write(stemloopid+"\n")
        f3.write(stemloopseq+"\n")
        f3.write(structure+"\n")
        f3.write(seq_dot+"\n")
    f1.close()
    f2.close()
    f3.close()


def gen_next_chunk(inputname, outputpath, start, num_of_lines):
    f = open(inputname, 'r')
    f.seek(start)
    outname = inputname +".chunk"+str(start)
    outchunkfile = open(outname,'w')
    count = 0
    while count < num_of_lines:
        line = f.readline()
        if not line:
            break
        outchunkfile.write(line)
        count += 1
        if count == num_of_lines:
            outchunkfile.close()
            return  f.tell(), outname, count
        else:
            continue
    if count >0 :
        outchunkfile.close()
        return  f.tell(), outname, count
    else:
        os.unlink(outname)
        return None


def fold_use_RNALfold(inputfastalist, tempfolder, dict_option, maxspan, chunksize):
    '''
    Fold the sequences using RNALfold and write results to outputname file.
    '''

    def fold(inputname, outputname, tempfolder, recovername, chunksize):
        command = "RNALfold -L " + str(maxspan)
        try:
            recoverfile = open(recovername, 'r')
            dict_recover_infold = cPickle.load(recoverfile)
            recoverfile.close()
            next_chunk = gen_next_chunk(inputname, tempfolder,  dict_recover_infold["nextpos"], chunksize)
            while next_chunk:
                rnaoutname = next_chunk[1] + ".rnafoldout"
                outputfile = open(rnaoutname, 'w')
                inputfile = open(next_chunk[1])
                write_formatted_string("[In fold]: Folding " + inputname, 30, sys.stdout)
                subprocess.check_call(command.split(), stdin=inputfile, stdout=outputfile)
                outputfile.close()
                recoverfile = open(recovername, 'r')
                dict_recover_infold = cPickle.load(recoverfile)
                recoverfile.close()
                dict_recover_infold["nextpos"] = next_chunk[0]  # the position for next round.
                dict_recover_infold["finishedchunks"] = dict_recover_infold["finishedchunks"] + 1
                dict_recover_infold["finishedseq"] = dict_recover_infold["finishedseq"] + next_chunk[2]/2
                dict_recover_infold["rnafoldoutfiles"].append(rnaoutname)
                recovername_temp = recovername + ".temp"
                recoverfile_temp = open(recovername_temp, 'w')
                cPickle.dump(dict_recover_infold, recoverfile_temp)
                recoverfile_temp.flush()
                os.fsync(recoverfile_temp.fileno())
                recoverfile_temp.close()
                os.rename(recovername_temp, recovername)  # This is atomic on Unix, not work on Windows
                write_formatted_string("[In fold]: " + inputname + ": "+str(dict_recover_infold["finishedseq"]) +" sequences folded. Checkpoint done.", 30, sys.stdout)
                # remove current chunk
                os.unlink(next_chunk[1])
                next_chunk = gen_next_chunk(inputname, tempfolder,  next_chunk[0], chunksize)
            #all chunk done, combine file
            outputname_temp = outputname+".temp"
            fout = open(outputname_temp, 'w')
            recoverfile = open(recovername, 'r')
            dict_recover_infold = cPickle.load(recoverfile)
            recoverfile.close()
            # cat the rnafold fold result together
            for name in dict_recover_infold["rnafoldoutfiles"]:
                with open(name) as f:
                    for line in f:
                        fout.write(line)
            fout.flush()
            os.fsync(fout.fileno())
            fout.close()
            os.rename(outputname_temp, outputname)  # This is atomic on Unix, not work on Windows
            #remove the rnafold result file
            #for name in dict_recover_infold["rnafoldoutfiles"]:
            #   os.unlink(name)

        except subprocess.CalledProcessError:
            sys.stderr.write("Error occurred when folding sequences (when calling RNALfold).\n")
            sys.stderr.write("Input file: "+inputname +", current chunk:")
            sys.exit(-1)
    inputfastalist.sort()
    outputnames = []
    for i in range(len(inputfastalist)):
        outname = os.path.join(tempfolder, dict_option["NAME_PREFIX"]+"_rnalfoldoutput_"+str(i))
        outputnames.append(outname)
    jobs = []
    for i in range(len(inputfastalist)):
        p = multiprocessing.Process(target = fold, args=(inputfastalist[i], outputnames[i], tempfolder, inputfastalist[i]+".recover", chunksize))
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

def run_check_sam_format(dict_option):
    write_formatted_string_withtime("Checking SAM file format", 30, sys.stdout)
    sam_not_right = False
    for name in dict_option["ALIGNMENT_FILE"]:
        write_formatted_string("Checking " + name, 30, sys.stdout)
        ret = check_sam_format(name)
        if ret[0] and ret[1] and ret[2]:
            write_formatted_string("*** SAMfile OK.", 30, sys.stdout)
            continue
        else:
            sam_not_right = True
        if not ret[0]:
            write_formatted_string("!!! SAMfile does not has header.", 30, sys.stdout)
        if not ret[1]:
            write_formatted_string("!!! SAMfile read IDs are not in the right format, please refer to the README.", 30, sys.stdout)
        if not ret[2]:
            write_formatted_string("!!! SAMfile flags are not in the right format, please refer to the README.", 30, sys.stdout)
    if sam_not_right:
        write_formatted_string("!!! Please make sure the SAM format is in the right format as described in the README.", 30, sys.stdout)
        return False
    return True

def run_check_fasta_format(dict_option):
    write_formatted_string_withtime("Checking reference sequence format.", 30, sys.stdout)
    ret = check_reference(dict_option["FASTA_FILE"])
    if ret[0]:
        write_formatted_string("*** Reference file OK.", 30, sys.stdout)
        return True
    else:
        if ret[1]==0:
            write_formatted_string("!!! The format of the input fasta file in not correct. Because the pipeline uses samtools to manipulate fasta file, there are some requirement for the fasta file:", 30, sys.stdout)
            write_formatted_string("1. The ID line of each sequence must start with '>'.", 30, sys.stdout)
            write_formatted_string("2. The sequence ID starts from the fisrt non-whitespace character until the first whitespace character. All characters after the first whitespace character are ignored. The ID of each sequence must only contain 0-9, a-z, A-Z, underscore(_), and dash(-).", 30, sys.stdout)
            write_formatted_string("3. The ID of each sequence must be unique.", 30, sys.stdout)
            write_formatted_string("4. Each sequence line of one sequence must be the same length, except the last line.", 30, sys.stdout)
            write_formatted_string("5. There should be no blank lines.", 30, sys.stdout)
        else:
            write_formatted_string("!!! " + ret[1], 30, sys.stdout)
        return False

def run_check_gff(gffname):
    write_formatted_string_withtime("Checking the format of the GFF file.", 30, sys.stdout)
    gffsize = os.path.getsize(gffname)/1024/1024
    if gffsize > 50:
        write_formatted_string_withtime("!!! Warning: large GFF file size: " + str(gffsize) + "MB. Make sure the GFF file only contains regions that need to be excluded from miRNA prediction.", 30, sys.stdout)
    ret = check_gff(gffname)
    if ret[0]:
        write_formatted_string("*** GFF file OK.\n", 30, sys.stdout)
        return True
    else:
        write_formatted_string("!!! " + ret[1], 30, sys.stdout)
        return False

def run_check(dict_option, outtempfolder, recovername):
    write_formatted_string_withtime("Checking RNALfold and samtools", 30, sys.stdout)
    if not check_RNALfold():
        write_formatted_string("!!! RNALfold is required but not installed or not in the PATH.", 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file for how to install and configuration RNALfold.", 30, sys.stdout)
    else:
        RNALfoldversion = get_RNALfold_version()

        if is_bug_RNALfold(RNALfoldversion.strip()):
            write_formatted_string("!!! The version of RNALfold (2.0.4) has a bug. ", 30, sys.stdout)
            write_formatted_string("!!! Please use the latest version of RNALfold", 30, sys.stdout)
            write_formatted_string("!!! Please refer to the README file for how to install and configurate RNALfold.", 30, sys.stdout)
        else:
            write_formatted_string("*** RNALfold is ready.", 30, sys.stdout)
        write_formatted_string("!!! RNALfold version: " + RNALfoldversion, 30, sys.stdout)

    samtoolcheck =  check_samtools()
    if not samtoolcheck[0]:
        write_formatted_string("!!! "+samtoolcheck[1], 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file for how to install and configurate samtools.", 30, sys.stdout)
        exit(-1)
    else:
        write_formatted_string("*** SAMtools is ready.", 30, sys.stdout)
        write_formatted_string("*** SAMtools version: " + samtoolcheck[1] + "\n", 30, sys.stdout)

    write_formatted_string_withtime("Checking SAM file format", 30, sys.stdout)
    sam_not_right = False
    for name in dict_option["ALIGNMENT_FILE"]:
        write_formatted_string("Checking " + name, 30, sys.stdout)
        ret = check_sam_format(name)
        if ret[0] and ret[1] and ret[2]:
            write_formatted_string("!!! SAMfile OK.", 30, sys.stdout)
            continue
        else:
            sam_not_right = True
        if not ret[0]:
            write_formatted_string("!!! SAMfile does not has header.", 30, sys.stdout)
        if not ret[1]:
            write_formatted_string("!!! SAMfile read IDs are not in the right format, please refer to the README.", 30, sys.stdout)
        if not ret[2]:
            write_formatted_string("!!! SAMfile flags are not in the right format, please refer to the README.", 30, sys.stdout)
    if sam_not_right:
        write_formatted_string_withtime("!!! Please make sure the SAM format is in the right format as descripted in the README.", 30, sys.stdout)

    run_check_fasta_format(dict_option)

    if os.path.exists(dict_option["GFF_FILE_EXCLUDE"]):
        run_check_gff(dict_option["GFF_FILE_EXCLUDE"])

    if os.path.exists(dict_option["GFF_FILE_INCLUDE"]):
        run_check_gff(dict_option["GFF_FILE_INCLUDE"])

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
        write_formatted_string("!!! Please refer to the README file for how to install and configurate RNALfold.", 30, sys.stdout)
        allgood = False
    else:
        RNALfoldversion = get_RNALfold_version()
        if is_bug_RNALfold(RNALfoldversion):
            write_formatted_string("!!! The version of RNALfold (2.0.4) has a bug. ", 30, sys.stdout)
            write_formatted_string("!!! Please use the latest version of RNALfold", 30, sys.stdout)
            write_formatted_string("!!! Please refer to the README file for how to install and configurate RNALfold.", 30, sys.stdout)
            allgood = False
        else:
            write_formatted_string("*** RNALfold is ready.", 30, sys.stdout)
        write_formatted_string("*** RNALfold version: " + RNALfoldversion, 30, sys.stdout)
    samtoolscheck =  check_samtools()
    if samtoolscheck[0]:
        write_formatted_string("*** SAMtools is ready.", 30, sys.stdout)
        write_formatted_string("*** SAMtools version: "+samtoolscheck[1] + "\n", 30, sys.stdout)
    else:
        #write_formatted_string("!!! samtools is required but not installed/not in the PATH/wrong version.", 30, sys.stdout)
        #write_formatted_string("!!! Please refer to the README file for how to install and configurate samtools.", 30, sys.stdout)
        write_formatted_string("!!! "+samtoolscheck[1], 30, sys.stdout)
        write_formatted_string("!!! Please refer to the README file for how to install and configurate samtools.", 30, sys.stdout)
        allgood = False

    return allgood


def get_samplename_from_sam(samname):
    samplename = ""
    with open(samname) as f:
        for line in f:
            if line.startswith("@"):
                continue
            else:
                samplename = "_".join(line.split()[0].split("_")[0:-2])
                return samplename


def write_failed_reason(dict_reasons, fname):
    outf = open(fname, 'w')
    for key in dict_reasons:
        outf.write(str(key))
        outf.write("\n")
    outf.close()

def run_prepare(dict_option, outtempfolder, recovername):
    logger = None
    if dict_option["LOG"]:
        logger = dict_option["LOGGER"]
        logger.info("Starting preparing data for the 'candidate' stage.")
    write_formatted_string_withtime("Starting preparing data for the 'candidate' stage.", 30, sys.stdout)
    combinedbamname, expandedsamname, expandedbamname, expandedbam_plus, expandedbam_minus = prepare_data(dict_option, outtempfolder, logger)
    samplenames = []
    for name in dict_option["ALIGNMENT_FILE"]:
        samplenames.append(get_samplename_from_sam(name))
    samplenamefilename = os.path.join(outtempfolder, "samplenames.dump")
    samplenamefile = open(samplenamefilename, 'w')
    cPickle.dump(samplenames, samplenamefile)
    samplenamefile.close()

    if logger:
        logger.info("Start writing recover infomation to the recovery file")
    dict_recover["last_stage"] = "prepare"
    dict_recover["finished_stages"][ "prepare"]={}
    dict_recover["finished_stages"]["prepare"]["combinedbamname"] = combinedbamname
    dict_recover["finished_stages"]["prepare"]["expandedbamname"] = expandedbamname
    dict_recover["finished_stages"]["prepare"]["expandedsamname"] = expandedsamname
    dict_recover["finished_stages"]["prepare"]["expandedbam_plus"] = expandedbam_plus
    dict_recover["finished_stages"]["prepare"]["expandedbam_minus"] = expandedbam_minus
    dict_recover["finished_stages"]["prepare"]["samplenamefilename"] = samplenamefilename
    dict_recover["files"]["prepare"] = []
    dict_recover["files"]["prepare"].append(combinedbamname)
    dict_recover["files"]["prepare"].append(expandedsamname)
    dict_recover["files"]["prepare"].append(expandedbamname)
    dict_recover["files"]["prepare"].append(expandedbam_plus)
    dict_recover["files"]["prepare"].append(expandedbam_minus)
    dict_recover["files"]["prepare"].append(samplenamefilename)
    recoverfile = open(recovername,'w')
    cPickle.dump(dict_recover,recoverfile)
    if logger:
        logger.info("Recovery file successfully updated.")
    write_formatted_string_withtime("Done (prepare stage)\n", 30, sys.stdout)
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

    if not dict_len:
        write_formatted_string_withtime("Can not get the sequence length from the input SAM files. Make sure SAM files have headers.", 30, sys.stdout)
        exit(-1)

    combinedbamname = dict_recover["finished_stages"]["prepare"]["combinedbamname"]
    #expandedbamname = dict_recover["finished_stages"]["prepare"]["expandedbamname"]
    #expandedsamname = dict_recover["finished_stages"]["prepare"]["expandedsamname"]
    expandedbam_plus = dict_recover["finished_stages"]["prepare"]["expandedbam_plus"]
    expandedbam_minus = dict_recover["finished_stages"]["prepare"]["expandedbam_minus"]
    #if the length of a contig is smaller than min_contig_len, then ignore it
    min_contig_len = 19
    write_formatted_string_withtime("Generating peaks.", 30, sys.stdout)
    depthfilename, dict_contigs = gen_contig_typeA(expandedbam_plus,expandedbam_minus, dict_option, min_contig_len, logger, outtempfolder)
    write_formatted_string_withtime( "All peaks generated.", 30, sys.stdout)
    write_formatted_string_withtime( "Combining peaks to form candidate regions.", 30, sys.stdout)
    dict_loci, piece_info = gen_candidate_region_typeA(dict_contigs,dict_len,dict_option,outtempfolder, logger)
    #  a list of (rnafoldinfasta, lociinfodump)  pairs
    write_formatted_string_withtime( "All candidate regions generated.", 30, sys.stdout)
    write_formatted_string_withtime( "Getting sequences of candidate regions.", 30, sys.stdout)
    num_loci, num_fasta, retnames = dump_loci_seqs_and_alignment_multiprocess(dict_loci,
                                                                              piece_info,
                                                                              combinedbamname,
                                                                              dict_option["FASTA_FILE"],
                                                                              rnalfoldinputprefix,
                                                                              outputdumpprefix,
                                                                              dict_option['READS_DEPTH_CUTOFF'] * 0.5,
                                                                              logger)
    write_formatted_string_withtime( "All sequences generated.", 30, sys.stdout)
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
    write_formatted_string_withtime("Done (candidate stage)\n", 30, sys.stdout)
    if logger:
        logger.info(outstr)
    if logger:
        logger.info("=========================Done (candidate stage)=======================\n\n")


def run_fold(dict_option, outtempfolder, recovername, tryrecover):
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

    for name in sorted(dict_recover["finished_stages"]["candidate"]["fasta"]):
        dict_recover_fold = load_recover_file(name+".recover")
        if dict_recover_fold and tryrecover:
            write_formatted_string_withtime("File "+name +" has  "+str(dict_recover_fold["finishedseq"]) +" sequences folded when stopped last time, start from the next sequence.", 30, sys.stdout)
        else:
            dict_recover_fold = {"nextpos":0,
                                 "finishedchunks":0,
                                 "finishedseq":0,
                                 "rnafoldoutfiles":[]
            }
            foldrecovername_temp = name+".recover.temp"
            foldrecovername = name+".recover"
            foldrecoverfile = open(foldrecovername_temp,'w')
            cPickle.dump(dict_recover_fold, foldrecoverfile)
            foldrecoverfile.flush()
            os.fsync(foldrecoverfile.fileno())
            foldrecoverfile.close()
            os.rename(foldrecovername_temp, foldrecovername)

    foldnames = fold_use_RNALfold(dict_recover["finished_stages"]["candidate"]["fasta"],
                                  outtempfolder,
                                  dict_option,
                                  dict_option["PRECURSOR_LEN"],
                                  2*dict_option["CHECKPOINT_SIZE"])  # Do checkpointing every CHECKPOINT_SIZE sequences.

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
    write_formatted_string_withtime("Done (fold stage)\n", 30, sys.stdout)
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
    samplenamefilename = dict_recover["finished_stages"]["prepare"]["samplenamefilename"]
    samplenames = cPickle.load(open(samplenamefilename))
    combinedsortedbamname = os.path.join(outtempfolder,"combined.filtered.sort.bam")

    allow_3nt_overhang = dict_option['ALLOW_3NT_OVERHANG']
    allow_no_star = dict_option['ALLOW_NO_STAR_EXPRESSION']
    output_details = dict_option['OUTPUT_DETAILS_FOR_DEBUG']
    min_mature_len = dict_option['MIN_MATURE_LEN']
    max_mature_len = dict_option['MAX_MATURE_LEN']

    foldnames = dict_recover["finished_stages"]["fold"]["foldnames"]
    result, failed_reasons = gen_miRNA_loci_nopredict(dict_recover["finished_stages"]["candidate"]["infodump"],
                                                      foldnames, combinedsortedbamname,
                                                      samplenames,  allow_3nt_overhang,
                                                      allow_no_star, output_details,
                                                      min_mature_len,
                                                      max_mature_len,
                                                      55,
                                                      dict_option['READS_DEPTH_CUTOFF'],
                                                      logger)
    failedname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_reason_why_not_miRNA.txt")
    faileddumpname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_failed.dict.dump")
    if failed_reasons:
        failedreadmappingfolder = os.path.join(dict_option['OUTFOLDER'],"failed_readmapping")
        if not os.path.exists(failedreadmappingfolder):
            os.mkdir(failedreadmappingfolder)

        formatted_fail_reasons = convert_failure_reasons_list(failed_reasons)
        write_dict_reasons(formatted_fail_reasons, failedname, faileddumpname, dict_option["FASTA_FILE"],
                           os.path.join(outtempfolder,
                                        "combined.filtered.sort.bam"),
                           samplenames, failedreadmappingfolder)

    if len(result)==0:
        write_formatted_string_withtime("0 miRNA identified. No result files generated.", 30, sys.stdout)
        if logger:
            logger.info("0 miRNA identified. No result files generated.")
        sys.exit(0)

    gffname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.gff3")
    maturename = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.mature.fa")
    stemloopname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.precursor.fa")
    ssname = os.path.join(dict_option["OUTFOLDER"],dict_option["NAME_PREFIX"]+"_miRNA.precursor.ss")


    write_formatted_string_withtime("The output files are in " + dict_option["OUTFOLDER"], 30, sys.stdout)
    if logger:
        logger.info("The output files are in " + dict_option["OUTFOLDER"])
        logger.info("Generating a gff file contains all predictions. GFF file name: " + gffname)
    adjust_mature_star(result)
    gen_gff_from_result(result,gffname)
    if logger:
        logger.info("Generating two fasta files contains predicted mature sequences and precursor sequences. Fasta file names: [mature]: " +
                    maturename + ", [precurosor]: " + stemloopname)

    readmappingfolder = os.path.join(dict_option['OUTFOLDER'],"readmapping")
    if not os.path.exists(readmappingfolder):
        os.mkdir(readmappingfolder)
    statname =  os.path.join(dict_option['OUTFOLDER'],"miRNA.stat.txt")
    htmlfile =  os.path.join(dict_option['OUTFOLDER'],dict_option['NAME_PREFIX']+"_miRNA.detail.html")
    csvfile =  os.path.join(dict_option['OUTFOLDER'],dict_option['NAME_PREFIX']+"_miRNA.detail.csv")
    if logger:
        logger.info("Generating precursor region read mapping files in folder: "+readmappingfolder)
    dict_mirna_info = gen_mirna_info(result, dict_option["FASTA_FILE"],
                                     os.path.join(outtempfolder,
                                                  "combined.filtered.sort.bam"),
                                     samplenames)
    gen_map_result(dict_mirna_info, readmappingfolder)
    if logger:
        logger.info("Generating html and csv format table files which contain the details of the predicted miRNAs. ")
        logger.info("Generating miRNA statistics file: "+statname)
    gen_csv_table(dict_mirna_info, csvfile)
    gen_html_table_file(dict_mirna_info, htmlfile)
    total_prediction, dict_len, dict_first = gen_miRNA_stat(dict_mirna_info)
    statfile = open(statname, 'w')
    statfile.write("Total number of predicted miRNAs: "+ str(total_prediction)+'\n')
    statfile.write("Distribution of the length of the mature miRNAs:\n")
    for k in sorted(dict_len.keys()):
        statfile.write(str(k) + ": " +str(dict_len[k]) + '\n')
    statfile.write("Distribution of the nucleotide of the first base of the mature miRNAs:\n")
    for k in sorted(dict_first.keys()):
        statfile.write(str(k) + ": " +str(dict_first[k]) + '\n')

    resultlistname = os.path.join(outtempfolder,dict_option["NAME_PREFIX"]+"_miRNA.info.dump")
    resultfile = open(resultlistname,'w')
    cPickle.dump(result,resultfile)

    dict_info_name = os.path.join(outtempfolder,dict_option["NAME_PREFIX"]+"_miRNA.dict.info.dump")
    resultfile = open(dict_info_name,'w')
    cPickle.dump(dict_mirna_info,resultfile)


    gen_mirna_fasta_ss_from_result(result, maturename, stemloopname, dict_option["FASTA_FILE"], ssname)
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
    outstr = "{0} miRNAs identified.".format(len(result))
    write_formatted_string(outstr, 30, sys.stdout)
    write_formatted_string_withtime("Done (predict stage)\n", 30, sys.stdout)
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
        sys.stderr.write("Error occurred when removing the tmp folder.\n")
        sys.stderr.write("You can remove the folder by yourself.\n")
        return 1


if __name__ == '__main__':
    dict_option = parse_option_optparse()
    display_dict_option(dict_option)
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
    if dict_option['TMPFOLDER']:  # user specified a tmp folder
        outtempfolder = dict_option['TMPFOLDER']
    recovername = os.path.join(outtempfolder, dict_option["NAME_PREFIX"]+"_recover")


    if dict_option['ACTION'] == 'check':
        run_check(dict_option, outtempfolder, recovername)
        exit(0)

    if not os.path.exists(outfolder):
        os.mkdir(dict_option["OUTFOLDER"])
    if not os.path.exists(outtempfolder):
        os.mkdir(outtempfolder)
    allgood = check_dependency()
    samgood = run_check_sam_format(dict_option)
    if not allgood or not samgood:
        exit(-1)
    sys.stdout.write("\n")

    fasta_good = run_check_fasta_format(dict_option)
    if not fasta_good:
        exit(-1)

    sys.stdout.write("\n")
    if os.path.exists(dict_option["GFF_FILE_EXCLUDE"]):
        gff_good = run_check_gff(dict_option["GFF_FILE_EXCLUDE"])
        if not gff_good:
            exit(-1)
    elif os.path.exists(dict_option["GFF_FILE_INCLUDE"]):
        gff_good = run_check_gff(dict_option["GFF_FILE_INCLUDE"])
        if not gff_good:
            exit(-1)
    else:
        write_formatted_string_withtime("No GFF file.", 30, sys.stdout)

    sys.stdout.write("\n")
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
        run_fold(dict_option, outtempfolder, recovername, False)
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
        run_fold(dict_option, outtempfolder, recovername, False)
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
            run_fold(dict_option, outtempfolder, recovername, False)
            run_predict(dict_option, outtempfolder, recovername)
            if not dict_option["KEEPTMP"]:
                run_removetmp(outtempfolder)
                if logger:
                    logger.info("Temporary folder removed.")
        elif last_stage == "candidate":
            write_formatted_string("*** Starting the pipeline from stage 'fold'.", 30, sys.stdout)
            run_fold(dict_option, outtempfolder, recovername, True)
            run_predict(dict_option, outtempfolder, recovername)
            if not dict_option["KEEPTMP"]:
                run_removetmp(outtempfolder)
                if logger:
                    logger.info("Temporary folder removed.")
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
