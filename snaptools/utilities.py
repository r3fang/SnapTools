import os
import sys
import gzip 
import bz2
import pysam
import argparse
import pybedtools
import re
import collections 

magic_dict = {
    "\x1f\x8b\x08": "gz",
    "\x42\x5a\x68": "bz2",
    "\x50\x4b\x03\x04": "zip"
    }

max_len = max(len(x) for x in magic_dict)

def file_type(filename):
    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "txt"
    
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def checkFileFormat(fname):
    """Check if a file format based on its extension.
    
    Args:
        fname: file name

    Return:
        file format bed, bam or snap
    """
    
    if fname.endswith('.bed') or fname.endswith('.bed.gz') or fname.endswith('.bed.bz2'):
        if is_bed(fname):
            return "bed";
    elif fname.endswith('.bam') or fname.endswith('.sam'):
        if is_bam(fname):
            return "bam";
    elif fname.endswith('.snap'):
        return "snap";
    elif fname.endswith('.gtf.gz') or fname.endswith('.gtf') or fname.endswith('.gff') or fname.endswith('.gff.gz'):
        if is_gff(fname):
            return "gtf";
    else:
        return "NA"
    
def is_bed(fname):
    """Check if a given file is a bed file
    
    Args:
        fname: a bed file name

    Returns:
        True if file is a bed file, otherwise False
    """
    try:
        # open fname using pysam
        fileType = pybedtools.BedTool(fname).file_type;
        if fileType == "bed":
            return True
    except IndexError as e:
        # handle the errors
        return False
    return False

def is_gff(fname):
    """Check if a given file is a gff file
    
    Args:
        fname: a gff file name

    Returns:
        True if file is a gff or gtf file, otherwise False
    """
    try:
        # open fname using pysam
        fileType = pybedtools.BedTool(fname).file_type;
        if fileType == "gff":
            return True
    except IndexError as e:
        # handle the errors
        return False
    return False

def is_bam(fname):
    """Check if a given file is a bam file
    
    Args:
        fname: file name

    Returns:
        True if file is a bam file, otherwise False
    """
    try:
        fileType = pybedtools.BedTool(fname).file_type;
        if fileType == "bam":
            return True
    except ValueError as e:
        # handle the errors
        return False
    return True

def isQuerynameSorted(fname, file_type):    
    """Check if a given file file is sorted based on the queryname.
    
    Args:
        fname: a input file name

    Returns:
        True if fname is a bam file sorted by read name, otherwise False
    """
    if file_type == "bam":
        return(isBamQuerynameSorted(fname))
    elif file_type == "bed":
        return(isBedQuerynameSorted(fname))
    else:
        print ('error: %s is not a supported file format!' % file_type)
        sys.exit(1)
        
def isBamQuerynameSorted(fname):
    """Check if a given bam file is sorted based on the queryname.
    
    Args:
        fname: bam file name

    Returns:
        True if fname is a bam file sorted by read name, otherwise False
    """
    if not is_bam(fname):
        print "Error @is_sorted_queryname: " + fname + " is not a bam/sam file!"
        return False;
    samfile = pysam.AlignmentFile(fname, "rb");
    flag = False
    header = samfile.header;
    if("HD" in header):
        if("SO" in header["HD"]):
            if(header["HD"]["SO"] == "queryname"):
                flag = True
    samfile.close();
    return flag

def isBedQuerynameSorted(fname):
    """Check if a given bed file is sorted based on the queryname (the 4th column).
    
    Args:
    -----
        fname: a bed file name

    Returns:
    -----
        True if fname is a bed file sorted by read name, otherwise False 
    """
    
    if not is_bed(fname):
        print "Error @is_bed_queryname_sorted: " + fname + " is not a bed file!"
        return False;
    
    # read the first 1 million reads and see if the read name is sorted
    fin = pybedtools.BedTool(fname);
    flag = True;
    num_instance = 0
    for item in fin:
        cur_barcode = str(item).split()[3];
        if num_instance == 0:
            prev_barcode = cur_barcode;
        num_instance += 1;
        if num_instance >= 1000000: break;
        if cur_barcode < prev_barcode: 
            flag = False;
            break;
    del fin
    return flag

def read_genome_size_from_bam(fname):
    """Read genome information from bam file header.
    
    Args:
        fname: bam file

    Returns:
        A dictionary contains SQ as key and SL as value, otherwise None
    """
    # first check if fname is a bam file
    if not is_bam(fname): 
        print "Error @read_genome_size_from_bam: " + fname + " is not a bam file!"
        return None;
    samfile = pysam.AlignmentFile(fname, "rb");
    header = samfile.header;
    res = None;
    try:
        header["SQ"];
        res = dict([(item["SN"], item["LN"]) for item in header["SQ"]]);
    except KeyError as e:
        res = None
    samfile.close()
    return res

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def readGenomeSizeFromTxt(fname):
    """
    Read genome information.
    
    Args:
    -----
        fname: a txt file contains genome information

    Returns:
    -----
        A dictionary contains SQ as key and SL as value, otherwise None
    """
    # first check if fname is a bam file
    res = dict();
    with open(fname) as fin:
        for line in fin:
            elems = line.split();
            chr_name = elems[0];
            chr_len = int(elems[1]);
            res[chr_name] = chr_len;
    return res

def readAlignmentInfo(fname, file_type):
    """
    Read the alingment information from a bam or bed file.
    
    Args:
    -----
        fname: a bam or bed file.
        
        file_type: file format.
    
    Returns:
    -----
        A dictionary contains alingment details. 
    """
    
    if file_type == "bam":
        res = readAlignmentInfoFromBam(fname)
    elif file_type == "bed":
        res = dict()
        res["PN"] = "NA"
        res["ID"] = "NA"
        res["VN"] = "NA"
        res["CL"] = "NA"
    else:
        print ('error: %s is not a supported file format!' % file_type)
        sys.exit(1)
    if "PN" not in res: res["PN"] = "NA";
    if "ID" not in res: res["ID"] = "NA";
    if "VN" not in res: res["VN"] = "NA";
    if "CL" not in res: res["CL"] = "NA";
    return res

def readAlignmentInfoFromBam(fname):
    """Read alingment info from bam file header.
    
    Args:
    -----
        fname: a bam file

    Returns:
    -----
        A dictionary contains alingment info, otherwise None
    """
    if not is_bam(fname): 
        print "Error @read_alignment_info: " + fname + " is not a bam file!"
        return None;

    samfile = pysam.AlignmentFile(fname, "rb");
    header = samfile.header;
    res = None;
    try:
        header["PG"];
        res = header["PG"][0]
    except KeyError as e:
        print "error: failed to read alingment info from bam file"
        sys.exit(1)
    samfile.close()
    return res
    
def minEditDistanceString(s1, S):
    """ Find string in list S that has the smallest edit distance to s1
    
    Args:
        s1: inqury string
        S:  A list of strings

    Returns:
        strings in S that has smallest edit distance to s1
    """
    d_min = sys.maxint
    # first find the min edit distance
    for s2 in S:
        d = editDistance(s1, s2)
        if d < d_min:
            d_min = d;

    # second find all strings in S that have d_min
    s_min = [];
    for s2 in S:
        d = editDistance(s1, s2)
        if d == d_min:
            s_min.append(s2)
    return s_min

def getBinsFromGenomeSize(genome_dict, bin_size):
    """Create a dictionary contains all bins of the same size across the genome
    
    Attributes:
        binSize: bin size (i.e. 5000)
    
        genomeDict: a dictionary contains chromosome sizes
    
    Return:
        A dictionary contains all bins and its index (start from 1)
    """
    bin_dict = collections.OrderedDict();
    i = 1;
    for _chrom in genome_dict:
    	for _start in range(1, genome_dict[_chrom], bin_size):
            _end = min(_start + bin_size - 1, genome_dict[_chrom]);
            _binId = (_chrom , _start, _end);
            bin_dict[_binId] = i;
            i = i +1;
    return bin_dict;

def editDistance(s1, s2):
    """ Calculate edit distance between 2 strings
    
    Args:
        s1: string 1
        s2: string 2

    Returns:
        Edit distance between s1 and s2
    """
    # ignore the case
    s1 = s1.upper()
    s2 = s2.upper()
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]
