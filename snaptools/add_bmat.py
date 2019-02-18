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
    from snaptools.snap import *
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


def snap_bmat(snap_file,
              bin_size_list,
              tmp_folder,
              verbose):

    """
    Pre-processing to create a snap file from a bam that contains alignments or a bed file that contains fragments.

    Args:
    --------
    snap_file: 
        a snap format file.
    
    Optional
    --------
    bin_size_list: 
            a list object contains all bin sizes [5000]
    
    verbose: 
            a boolen variable indicates whether to output the progress [True];
    """

    if not os.path.exists(snap_file):
        print 'error: ' + snap_file + ' does not exist!'
        sys.exit(1);
    
    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(snap_file);
    if file_format != "snap":
        print "error: input file %s is not a snap file!" % snap_file;
    
    # extract the barcodes
    barcode_dict = snaptools.snap.getBarcodesFromSnap(snap_file);

    # create the bin list
    f = h5py.File(snap_file, "a", libver='earliest');

    try:
        genome_dict = dict(zip(f["HD"]["SQ"]["SN"][:], f["HD"]["SQ"]["SL"][:]))
    except KeyError:
        print "error: unable to read genome information"
        sys.exit(1)
    
    if "AM" in f:
        print "error: AM - cell x bin accessibility matrix already exists, delete it first using --snap-del "
        sys.exit(1)

    bin_dict_list = collections.defaultdict(dict);

    for bin_size in bin_size_list:
        bin_dict = snaptools.utilities.getBinsFromGenomeSize(genome_dict, bin_size);
        bin_dict_list[bin_size] = bin_dict;

    num_barcode = len(barcode_dict);
    if verbose:
        print "===== reading the barcodes and bins ======"    
        print "@AM\tnBinSize:%d"%len(bin_dict_list.keys())
        print "@AM\tbinSizeList:", bin_dict_list.keys();    
        for bin_size in bin_dict_list.keys():
            print "@AM\tbinSize:%d\tnBin:%d"%(bin_size, len(bin_dict_list[bin_size]));
        
    idxList   = collections.defaultdict(list);        # barcode index list
    idyList   = collections.defaultdict(list);        # bin index list
    countList = collections.defaultdict(list);        # number of count
            
    barcode_id = 0
    for barcode in f["BD"]["name"]:        
        _chroms = f["FM"]["fragChrom"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
        _start = f["FM"]["fragStart"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
        _len = f["FM"]["fragLen"][(f["FM"]["barcodePos"][barcode_id] - 1):(f["FM"]["barcodePos"][barcode_id] + f["FM"]["barcodeLen"][barcode_id] - 1)]
        frag_list_uniq = zip(_chroms, _start, _start + _len);    

        for bin_size in bin_dict_list:
            bin_dict = bin_dict_list[bin_size];
            bins = collections.defaultdict(lambda : 0);
            for item in frag_list_uniq:
                bin_chr = item[0];
                for bin_pos in set([item[1]/bin_size * bin_size + 1, item[2]/bin_size * bin_size + 1]):
                    bins[(bin_chr, bin_pos, bin_pos + bin_size - 1)] += 1;
        
            for key in bins:
                if key in bin_dict and barcode in barcode_dict:
                    idyList[bin_size].append(bin_dict[key]);
                    countList[bin_size].append(bins[key]);
                    idxList[bin_size].append(barcode_dict[barcode].id);
        
        barcode_id += 1;
        del bin_dict, bins, frag_list_uniq;
    
    dt = h5py.special_dtype(vlen=bytes)    
    f.create_dataset("AM/nBinSize", data=len(bin_dict_list),  dtype="uint32");
    f.create_dataset("AM/binSizeList", data=bin_dict_list.keys(),  dtype="uint32");
    
    for bin_size in bin_dict_list:
        f.create_dataset("AM/"+str(bin_size)+"/binChrom",data=[key[0] for key in bin_dict_list[bin_size]], dtype=h5py.special_dtype(vlen=bytes), compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/binStart",data=[key[1] for key in bin_dict_list[bin_size]], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/idx", data=idxList[bin_size], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/idy", data=idyList[bin_size], dtype="uint32", compression="gzip", compression_opts=9);
        f.create_dataset("AM/"+str(bin_size)+"/count", data=countList[bin_size], dtype="uint8", compression="gzip", compression_opts=9);    
    f.close()  
    return 0
