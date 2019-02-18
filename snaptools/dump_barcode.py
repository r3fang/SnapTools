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


def dump_barcode(snap_file,
                 output_file,
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
        a txt format file contains attributes of selected barcodes.
    
    Optional Args:
    --------
                               
    barcode_file:
        a txt file contains selected barcodes.
                  
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
    
    # check if BD session exists
    f = h5py.File(snap_file, "r", libver='earliest');
    if "BD" not in f:
        print "error: BD session does not exit in the snap file!";
        sys.exit(1);
    f.close();
    
    barcode_dict = getBarcodesFromSnap(snap_file);
    
    # identify the barcodes
    if barcode_file is None:
        barcode_list = barcode_dict.keys();
    else:
        barcode_list = getBarcodesFromTxt(barcode_file).keys();
    
    # write fragments down
    fout = open(output_file, "w");
    res = ["Barcode", "TN", "UM", "SE", "SA", "PE", "PP", "PL", "US", "UQ", "CM", "\n"];
    fout.write("\t".join(map(str, res)))    
    for barcode in barcode_list:
        if barcode in barcode_dict:
            res = [barcode, barcode_dict[barcode].total,  \
                     barcode_dict[barcode].mapped, \
                     barcode_dict[barcode].single, \
                     barcode_dict[barcode].secondary, \
                     barcode_dict[barcode].proper_paired, \
                     barcode_dict[barcode].proper_flen, \
                     barcode_dict[barcode].usable, \
                     barcode_dict[barcode].uniq, \
                     barcode_dict[barcode].chrM, \
                     "\n"]
        else:
            res = [barcode] +  ["0"] * 9
        fout.write("\t".join(map(str, res)))          
    return 0
    
