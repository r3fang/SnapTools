# -*- coding: utf-8 -*-
""" 

The MIT License

Copyright (c) 2018 Rongxin Fang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""

import sys, os 
import collections
import gzip
import operator
import argparse
import datetime
import urllib2
import time
import subprocess
import shlex
import tempfile

try:
    import snaptools.utilities
    import snaptools.global_var
    import snaptools.gtf
    from snaptools.utilities import file_type
except Exception:
    print "Package snaptools not installed!"
    sys.exit(1)

try:
	import numpy as np
except Exception:
    print "Package numpy not installed!"
    sys.exit(1)

try:
	import h5py
except Exception:
    print "Package numpy not installed!"
    sys.exit(1)

try:
    import pysam
except Exception:
    print "Package pysam not installed!"
    sys.exit(1)

try:
    import pybedtools
except Exception:
    print "Package pybedtools not installed!"
    sys.exit(1)


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

def getBarcodesFromSnap(fname):
    """Read barcodes from a snap file
    
    Attributes:
        fname - a snap-format file

    Return:
        a dictionary contains all barcode in the snap file
    """
    try:
        f = h5py.File(fname, 'r');
    except IOError:
        print("error: unable to open fname, check if it is a snap file")
        sys.exit(1)

    barcode_dict = collections.OrderedDict();
    i = 1;
    for item in f["BD/name"]:
        barcode_dict[item] = qc();         
        barcode_dict[item].id = i; 
        barcode_dict[item].total = f["BD"]["TN"][i-1]
        barcode_dict[item].mapped = f["BD"]["UM"][i-1]
        barcode_dict[item].single = f["BD"]["SE"][i-1]
        barcode_dict[item].secondary = f["BD"]["SA"][i-1]
        barcode_dict[item].paired = f["BD"]["PE"][i-1]
        barcode_dict[item].proper_paired = f["BD"]["PP"][i-1]
        barcode_dict[item].proper_flen = f["BD"]["PL"][i-1]
        barcode_dict[item].usable = f["BD"]["US"][i-1]
        barcode_dict[item].uniq = f["BD"]["UQ"][i-1]
        barcode_dict[item].chrM = f["BD"]["CM"][i-1]
        barcode_dict[item].final = f["BD"]["UQ"][i-1] - f["BD"]["CM"][i-1]        
        i = i + 1;
    f.close()
    return barcode_dict

def getBarcodesFromInput(fname, ftype):
    """Read barcodes from a given input file (bam or bed)
    
    Attributes:
        fname: 
            A bam or bed file
        ftype:
            File format, bed or bam
    Returns:
        a dictionary contains all barcode present in the input file
    """    
    # automatically determine the file format
    if ftype == "bam":
        return getBarcodesFromBam(fname);
    if ftype == "bed":
        return getBarcodesFromBed(fname);
    else:
        print "error: unrecognized file format for %s, only support bam or bed file" % ftype
        sys.exit(1)

def getBarcodesFromTxt(fname):
    """Read barcodes from a given txt file
    
    Attributes:
        fname: A txt file that contains pre-defined barcodes
                
    Returns:
        a dictionary contains barcode in the txt file
    """   
    barcode_list = [];     
    with open(fname) as fin:
        for line in fin:
            if line.startswith("#"): continue;
            if(len(line.split()) > 1):
                barcode = line.split()[0].upper()
                barcode_list.append(barcode)            
    barocde_num = len(set(barcode_list));
    if barocde_num < len(barcode_list):
        print "warning: duplicate barcodes identified, remove it first"
        sys.exit(1)
    barcode_dict = collections.OrderedDict(zip(sorted(set(barcode_list)), range(1,(barocde_num+1))));
    return barcode_dict

def getBarcodesFromBam(input_bam):
    """Identify unique barcodes from the bam file
    
    Args:
        input_bam: a bam file

    Returns:
        A dictionary contains all barcodes, otherwise None
    """
    
    barcode_dict = collections.OrderedDict();
    samfile = pysam.AlignmentFile(input_bam, "rb");
    i = 1
    for _read in samfile:
        barcode = _read.qname.split(":")[0].upper();
        if barcode not in barcode_dict:
            barcode_dict[barcode] = i
            i = i + 1
    samfile.close();
    return barcode_dict;

def getBarcodesFromBed(input_bed):
    """Identify unique barcodes from a bed file
    
    Args:
        input_bed: a bed file.
    
    Returns:
        A dictionary contains all barcodes in the bed file.
    
    """
    
    # check the file type
    barcode_dict = collections.OrderedDict();   
    if snaptools.utilities.file_type(input_bed) == "gz":
        fin = gzip.open(input_bed, 'rb');
    elif snaptools.utilities.file_type(input_bed) == "bz2":
        fin = bz2.BZ2File(input_bed, 'r');
    elif snaptools.utilities.file_type(input_bed) == "txt":
        fin = open(input_bed, 'r');
    else:
	    print "error: unrecoginized bed file format, only supports .gz, .bz2, .bed"
	    sys.exit(1)    
    i = 1

    for _read in fin:
        barcode = _read.split()[3].split(":")[0].upper();
        if barcode not in barcode_dict:
            barcode_dict[barcode] = i
            i = i + 1
    fin.close();
    return barcode_dict;

def getBarcodeCov(barcode_list, input_file, file_format):
    """
    Get barcode coverage from a given input file
    
    Args:
        barcode_dict: 
            a list of pre-defined barcodes

        input_file: 
            a bam or bed file

        file_format: 
            file format bed or bam
    
    Returns:
        a dictionary contains barcode coverage
    """    
    
    if file_format == "bam":
        return get_barcode_cov_from_bam(barcode_list, input_file)
    elif file_format == "bed":
        return get_barcode_cov_from_bed(barcode_list, input_file)
    else:
        print "error: unrecognized file format %s " % file_format
        sys.exit(1)
    
def get_barcode_cov_from_bam(barcode_list, input_bam):
    """
    Get barcode coverage from bam file
    
    Args:
        barcode_dict: a list of pre-defined barcodes

        input_bam: a bam file
    
    Returns:
        a dictionary contains barcode coverage
    """    
    
    if not snaptools.utilities.is_bam(input_bam):
        print "error: @get_barcode_cov_from_bam: " + input_bam + " is not a bam file!";   
        sys.exit(1)
    
    if len(barcode_list) == 0:
        print "error: @get_barcode_cov_from_bam: barcode_list is empty!";   
        sys.exit(1)
            
    barcode_dict = collections.defaultdict(lambda : 0);    
    samfile = pysam.AlignmentFile(input_bam, "rb");
    for _read in samfile:
        barcode = _read.qname.split(":")[0];
        # approximate counting, a read is half fragment
        barcode_dict[barcode] += 0.5;
    return barcode_dict;

def get_barcode_cov_from_bed(barcode_list, input_bed):
    """
    Get barcode coverage from bed file
    
    Args:
    -----
    barcode_dict: 
        a list of pre-defined barcodes

    input_bed: 
        a bed file
    
    Returns:
    ------
    a dictionary contains barcode coverage
    """    
        
    if len(barcode_list) == 0:
        print "error: @get_barcode_cov_from_bam: barcode_list is empty!";   
        sys.exit(1)
    
    if file_type(input_bed) == "gz":
        fin = gzip.open(input_bed, 'rb');
    elif file_type(input_bed) == "bz2":
        fin = bz2.BZ2File(input_bed, 'r');
    elif file_type(input_bed) == "txt":
        fin = open(input_bed, 'r');
    else:
	    print "error: unrecoginized bed file format, only supports .gz, .bz2, .fastq"
	    sys.exit(1)
    
    barcode_dict = collections.defaultdict(lambda : 0);    
    for _read in fin:
        barcode = _read.split()[3].split(":")[0].upper();
        # approximate counting, a read is half fragment
        barcode_dict[barcode] += 1;
    fin.close()
    return barcode_dict;

def group_reads_by_barcode(input_file, file_format):
    if file_format == "bam":
        return group_reads_by_barcode_bam(input_file);
    elif file_format == "bed":
        return group_reads_by_barcode_bed(input_file);
    else:
        print "error: unrecoganized file format for %s " % file_format;
        sys.exit(1);
    
def group_reads_by_barcode_bam(input_bam):
    """ Group reads based on the barcodes
    
    Args:
        input_bam: a bam file

    Returns:
        Generator that contains reads sharing the same barcode
    """
    if not os.path.exists(input_bam): 
        print "Error @group_reads_by_barcode_bam: " + input_bam + " does not exist!"
 
    if not snaptools.utilities.is_bam(input_bam): 
        print "Error @group_reads_by_barcode_bam: " + input_bam + " is not a bam file!"    
 
    read_group_list = []; 
    pre_barcode = "";
    samfile = pysam.AlignmentFile(input_bam, "rb");
    for cur_read in samfile:
        cur_barcode = cur_read.qname.split(":")[0];
        if cur_barcode == pre_barcode:
            read_group_list.append(cur_read)
        else:
            if pre_barcode != "":
                # return read group
                yield (x for x in read_group_list)
            read_group_list = [cur_read] # add the first read
            pre_barcode = cur_barcode
    # reads from the last barcode
    yield (x for x in read_group_list)

def group_reads_by_barcode_bed(input_bed):
    """ Group fargments based on the barcodes
    
    Args:
        input_bed: a bed file

    Returns:
        Generator that contains reads sharing the same barcode
    """
    if not os.path.exists(input_bed): 
        print "Error @group_reads_by_barcode_bam: " + input_bed + " does not exist!"
 
    read_group_list = []; 
    pre_barcode = "";
    
    if file_type(input_bed) == "gz":
        fin = gzip.open(input_bed, 'rb');
    elif file_type(input_bed) == "bz2":
        fin = bz2.BZ2File(input_bed, 'r');
    elif file_type(input_bed) == "txt":
        fin = open(input_bed, 'r');
    else:
	    print "error: unrecoginized fastq file format, only supports .gz, .bz2, .fastq"
	    sys.exit(1)
    
    for cur_read in fin:
        cur_barcode = cur_read.split()[3].split(":")[0].upper();
        if cur_barcode == pre_barcode:
            read_group_list.append(cur_read)
        else:
            if pre_barcode != "":
                # return read group
                yield (x for x in read_group_list)
            read_group_list = [cur_read] # add the first read
            pre_barcode = cur_barcode
    # reads from the last barcode
    yield (x for x in read_group_list)
    fin.close()

def pairReadsByName(read_list):
    """ Pair reads based on read names
    
    Args:
        read_list: a list of reads that share the same barcode

    Returns:
        Generator contains read pairs from the same fragment
        and a boolen variable indicates whether it is supplementary alignment
    """
    # pair up 
    for read1 in read_list:
        # read until read1 is not a supplementary alignment
        while(read1.is_supplementary):
            yield (read1, None, False, True)
            try:
                #print "Warning: skipping ", read1.qname;
                read1 = read_list.next();
            except:
            	break        
        try:
            read2 = read_list.next();
        except:
        	break
        while(read2.is_supplementary):
            yield (read2, None, False, True)
            try:
                #print "Warning: skipping ", read2.qname;
                read2 = read_list.next();
            except:
            	break
        if(read1.qname != read2.qname):
            while (read1.qname != read2.qname):
                yield(read1, None, False, False);
                read1 = read2;
                try:
                    read2 = read_list.next();
                    while(read2.is_supplementary):
                        try:
                            #print "Warning: skipping ", read2.qname;
                            read2 = read_list.next();
                        except:
                        	break
                except:
                    break;
        yield (read1, read2, True, False)
 
def readPairToFragment(read1, read2, is_secondary):
    """ convert read pairs to fragments
    
    Args:
        read1: R1 read

        read2: R2 read
    Returns:
        Generator that contains a fragment object
    """
    try:
        read1.qname;
        read2.qname;
        if read1.qname != read2.qname:
            sys.exit('read_pair_to_fragment: read1 and read2 name does not match!');  
    except ValueError as e:
        sys.exit('read_pair_to_fragment: can not get read1 or read2 name!');   
    barcode = read1.qname.split(":")[0];
    mapq = min(read1.mapq, read2.mapq);
    try:
        chrom1 = read1.reference_name;
        start1 = read1.reference_start;
        strand1 = "-" if read1.is_reverse else "+";    
        chrom2 = read2.reference_name;
        start2 = read2.reference_start;
        strand2 = "-" if read1.is_reverse else "+";    
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        flen2 = read2.reference_length if read2.reference_length != None else 0
        end1   = start1 + flen1;
        end2   = start2 + flen2;
        start = min(start1, end1, start2, end2);
        end = max(start1, end1, start2, end2);
        return fragment(read1.qname, chrom1, start, abs(start - end), mapq, False, is_secondary, read1.is_proper_pair);
    except ValueError as e:
        return fragment(read1.qname, None, None, None, mapq, False, is_secondary, False);

def readToFragment(read1, is_secondary):
    """ convert a single read to fragment
    
    Args:
        read1: a single read

    Returns:
        Generator contains a fragment object
    """
    try:
        read1.qname;
    except ValueError as e:
        sys.exit('readto_fragment: can not get read1 name!');   

    barcode = read1.qname.split(":")[0];
    mapq = read1.mapq;
    try:
        chrom1 = read1.reference_name;
        start1 = read1.reference_start;
        strand1 = "-" if read1.is_reverse else "+";    
        # it is possible that flen1 is None  
        flen1 = read1.reference_length if read1.reference_length != None else 0
        end1   = start1 + flen1;
        start = min(start1, end1);
        end = max(start1, end1);
        #(qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        return fragment(read1.qname, chrom1, start, abs(start - end), mapq, True, is_secondary, False);
    except ValueError as e:
        return fragment(read1.qname, None, None, None, mapq, True, is_secondary, False);

class qc(object):
    """A quality control object that has the following attributes:

    Attributes:
        total: total number of sequenced fragments.

        mapped: number of mappable fragments.

        chrM: number of fragments mapped to chrM.

        paired: number of fragments are paired.

        single: number of fragments are from single read.

        proper_paired: number of paired reads are properly paired.

        usable: number of usable fragments.

        uniq: number of unique fragments.

        isize: average insert size distribution.
    """
    def __init__(self):
        """Return a qc object"""
        self.id = 0
        self.total = 0
        self.mapped = 0
        self.single = 0
        self.secondary = 0
        self.paired = 0
        self.proper_paired = 0
        self.proper_flen = 0
        self.usable = 0
        self.uniq = 0
        self.chrM = 0
        self.final = 0

class fragment(object):
    """A fragment object that has the following attributes:

    Attributes:
        chrom: chromsome name

        start: start position

        end: end position

        mapq: mapping quality

        is_proper_pair: whether properly paired

        is_single: whether it is a single read

        is_secondary: whether it is a secondary alignment read

        flen: fragment length
    """
    def __init__(self, qname, chrom, pos, flen, mapq, is_single, is_secondary, is_proper_pair):
        """Return a qc object"""
        self.qname = qname
        self.chrom = chrom
        self.pos = pos
        self.flen = flen
        self.mapq = mapq
        self.is_single = is_single
        self.is_secondary = is_secondary
        self.is_proper_pair = is_proper_pair


def run_snap_view_head(input_snap):
    """
    View header session of a snap file

    Args:
    --------
        input_snap: a snap.
    """
    # 1. whether this file exists
    # 2. whether this is a snap-format file
    # 3. check the version of the snap-format file
    # 4. read the header info
    # 5. read the alignment info
    # 5. read the reference genome info
    # 4. print the header session on the screen 
    
    if not os.path.exists(input_snap):
        print 'error: ' + input_snap + ' does not exist!'
        sys.exit(1);

    #if not snaptools.utilities.is_snap(input_snap):
    #    print 'error: ' + input_snap + ' is not a snap format file!'
    #    sys.exit(1);
    
    f = h5py.File(input_snap, 'r');
    
    if "HD" not in f:
        print input_snap + " HD session does not exist!"
        sys.exit(1)
    
    # read the header info
    header = collections.defaultdict();
    header["MG"] = f["HD"]["MG"].value
    header["DT"] = f["HD"]["DT"].value
    header["VN"] = f["HD"]["VN"].value
    header["CL"] = f["HD"]["CL"].value
    header["CW"] = f["HD"]["CW"].value

    # read the alignment info
    align_dict = collections.defaultdict();
    align_dict["PN"] = f["HD"]["AL"]["PN"].value
    align_dict["ID"] = f["HD"]["AL"]["ID"].value
    align_dict["VN"] = f["HD"]["AL"]["VN"].value
    align_dict["CL"] = f["HD"]["AL"]["CL"].value
        
    # read the genome info
    genome_name = f["HD"]["SQ"]["ID"].value
    genome_dict = dict(zip(f["HD"]["SQ"]["SN"], f["HD"]["SQ"]["SL"]))
    
    print "@HD\tMG:%s"%header["MG"]
    print "@HD\tDT:%s"%header["DT"]
    print "@HD\tVN:%s"%header["VN"]
    print "@HD\tCL:%s"%header["CL"]
    print "@HD\tCW:%s"%header["CW"] 
    print "@AL\tPN:%s"%align_dict["PN"]
    print "@AL\tID:%s"%align_dict["ID"]
    print "@AL\tVN:%s"%align_dict["VN"]
    print "@AL\tVN:%s"%align_dict["CL"]
    print "@SQ\tID:%s"%(genome_name);
    for chr_name in genome_dict:
        print "@SQ\tSN:%s\tSL:%d"%(chr_name, genome_dict[chr_name]);
    f.close();

def getCellPeakVector(f, peak_file, peak_dict, barcode_dict, barcode, tmp_folder):
    """
    """
    idx_list = [];
    idy_list = [];
    count_list = [];
    barcode_id = barcode_dict[barcode];
    barcode = f["BD/name"][barcode_id];
    _chroms = f["FM"]["fragChrom"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
    _start = f["FM"]["fragStart"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
    _end = f["FM"]["fragLen"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)];
    frag_list_uniq = zip(_chroms, _start, _start + _end);
    
    # write the fragments as bed file in a temporary file
    fout_frag = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder);
    for item in frag_list_uniq:
        fout_frag.write("\t".join(map(str, item))+"\n");
    fout_frag.close();
    
    # find overlaps
    frag_bt = pybedtools.BedTool(fout_frag.name); 
    peak_bt = pybedtools.BedTool(peak_file);      
    
    # count for frequency            
    peak_cov = collections.defaultdict(lambda : 0);            
    for item in frag_bt.intersect(peak_bt, wa=True, wb=True):
        elems = str(item).split();
        key = (elems[3], int(elems[4]), int(elems[5]));
        peak_cov[key] += 1;
    
    # save the count matrix for current barcode            
    for key in peak_cov:
        if key in peak_dict and barcode in barcode_dict:
            idx_list.append(barcode_dict[barcode]);
            idy_list.append(peak_dict[key]);
            count_list.append(peak_cov[key]);
    
    # remove the tmporary file
    subprocess.check_call(["rm", fout_frag.name]);
    return((idx_list, idy_list, count_list))


def getFragFromBarcode(fname, barcode_list):
    """Extract fragments of given barcodes

    Attributes:
        fname: a snap file

        barcode_list: a list of selected barcodes
    
    Return:
        a list of fragments
    """

    # extract the fragments of given barcodes;
    barcode_dict = getBarcodesFromSnap(fname);
    # get the fragments of selected barcodes;
    f = h5py.File(fname, "r", libver='earliest');
    frag_list = []
    for barcode in barcode_list:
        if barcode not in barcode_dict: continue;
        barcode_id = barcode_dict[barcode].id;
        barcodePos = f["FM"]["barcodePos"][barcode_id-1] - 1
        barcodeLen = f["FM"]["barcodeLen"][barcode_id-1]
        _chroms = f["FM"]["fragChrom"][barcodePos:(barcodePos + barcodeLen)];
        _start = f["FM"]["fragStart"][barcodePos:(barcodePos + barcodeLen)];
        _len = f["FM"]["fragLen"][barcodePos:(barcodePos + barcodeLen)];        
        frag_list = frag_list + zip(_chroms, _start, _start + _len, [barcode] * len(_chroms));
    f.close()
    return frag_list
  
def dump_read(snap_file,
              output_file,
              buffer_size,
              barcode_file,
              tmp_folder,
              overwrite):

    """
    Dump reads from snap file into a bed file (of selected barcodes)

    Required Args:
    --------
    snap_file: 
        a snap format file.

    output_file: 
        a bed format file contains reads of selected barcodes.
    
    Optional Args:
    --------
    buffer_size:
        max number of barcodes whose reads are stored before dumping.
              
    barcode_file:
        a txt file contains selected barcodes.
    
    tmp_folder:
        a tmp folder that stores potential temoprial files.
              
    overwrite:
        a boolen variable indicates whether to overwrite a file if it already exists

    """
    # check if input snap file exists
    if not os.path.exists(snap_file):
        print 'error: %s does not exist!' % snap_file;
        sys.exit(1);

    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(snap_file);
    if file_format != "snap":
        print "error: input file %s is not a snap file!" % snap_file;

    # check the output bed file exists
    if os.path.exists(output_file):
        if overwrite == True:
            subprocess.check_call(["rm", output_file]);
        else:
            print 'error: %s already exists, change --overwrite or remove it first!' % output_file
            sys.exit(1);    

    # identify the barcodes
    if barcode_file is None:
        barcode_dict = getBarcodesFromSnap(snap_file);
    else:
        barcode_dict = getBarcodesFromTxt(barcode_file);
    
    # check if tmp_folder is given
    if(tmp_folder!=None):
        if(not os.path.isdir(tmp_folder)):
            print('error: %s does not exist!' % tmp_folder)
            sys.exit(1);
        
    # check if FM session exists
    f = h5py.File(snap_file, "r", libver='earliest');
    if "FM" not in f:
        print "error: FM session does not exit in the snap file!";
        sys.exit(1);
    f.close();
    
    # cut barcodes into small chunks
    barcode_list = snaptools.utilities.chunks(barcode_dict.keys(), buffer_size);
    barcode_list = [list(barcode_chunk) for barcode_chunk in barcode_list];
    nChunk = len(barcode_list);

    # write fragments down
    if output_file.endswith(".gz"):
        fout = gzip.open(output_file, "wb")
    elif output_file.endswith(".bz2"):
        fout = gzip.open(output_file, "wb")
    else:
        fout = open(output_file, "w")

    # cut the barcode into chunks and write down seperately
    for i in range(nChunk):
        # extract the fragment list of given barcodes
        frag_list = getFragFromBarcode(snap_file, barcode_list[i]);
        for item in frag_list:
            fout.write("\t".join(map(str, item)) + "\n")
        del frag_list
    fout.close()
    return 0
