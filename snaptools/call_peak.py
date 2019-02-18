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

def call_peak(snap_file,
              output_prefix,
              barcode_file,
              gsize,
              path_to_macs,
              buffer_size,
              macs_options,
              tmp_folder):
    
    """
    Call peaks using selected barcodes.

    Required:
    --------
    snap_file: 
        a snap format file;

    barcode_file: 
        a txt file contains selected barcodes as the first column;

    output_prefix: 
        experiment name, which will be used to generate output;

    path_to_macs: 
        a path to the folder contains excutable file macs2;

    Optional:
    --------
    buffer_size:
        max number of barcodes to be stored in the memory
    
    macs_options:
        a list of strings indicating options you'd like passed to aligner.
        (default: "--nomodel --qval 1e-2 -B --SPMR --call-summits --keep-dup all");

    tmp_folder:
        folder to store intermedia files;
    """
    
    # if the path_to_macs path given, need to check the existance of MACS1
    if path_to_macs != None:
        path_to_macs+="/"
        if not os.path.isdir(path_to_macs):
            print('Error: ' + path_to_macs + ' is not a folder')
            sys.exit(1);
        if not os.path.exists(path_to_macs+"macs2"):
            print('Error: macs2 does not exist')
            sys.exit(1);
    else:
        try:
            # pipe output to /dev/null for silence
            null = open("/dev/null", "w")
            subprocess.Popen(macs2, stdout=null, stderr=null)
            null.close()
        except OSError as e:
            print('Error: macs2 does not exist!');
            sys.exit(1);
        path_to_macs=""

    # check temp folder
    if(tmp_folder!=None):
        if not os.path.isdir(tmp_folder):
            print("Error: 'tmp_folder' is not a folder or does not exist")
            sys.exit(1);
    
    # check wheather snap file exists
    if not os.path.exists(snap_file):
        print 'error: ' + snap_file + ' does not exist!'
        sys.exit(1);

    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(snap_file);
    if file_format != "snap":
        print "Error: input file %s is not a snap file!" % snap_file;

    if not os.path.exists(barcode_file):
        print 'error: ' + barcode_file + ' does not exist!'
        sys.exit(1);
    
    # default aligner option
    if macs_options is None:
        macs_options = ["--nomodel", "--qval 1e-2", "-B", "--SPMR", "--call-summits", "--keep-dup all"];
    
    # extract the barcodes
    barcode_dict = snaptools.snap.getBarcodesFromSnap(snap_file);

    # read barcode from barcode file 
    barcode_sel = snaptools.snap.getBarcodesFromTxt(barcode_file);
        
    # compare selected barcodes and total barcodes in snap file
    if len(list(set(barcode_sel.keys()) & set(barcode_dict.keys()))) == 0:
        print 'Error: selected barcodes does not exist in the snap file!'
        sys.exit(1);
    
    # first cut the fragments into small piecies, write them down
    fout_frag = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder);
    snaptools.snap.dump_read(snap_file, fout_frag.name, buffer_size, barcode_file, tmp_folder, True);

    # call peaks using macs
    args = [path_to_macs+"macs2"];
    args.append("callpeak");
    args.append("-f AUTO");
    args.append("-t " + fout_frag.name);
    args.append("-n " + output_prefix);
    args.append("-g " + gsize);
    args.extend(macs_options);
    ftmp = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)

    try:
        subprocess.check_call(" ".join(args), stdout=ftmp, shell=True, executable='/bin/bash');
    except subprocess.CalledProcessError as e:
        sys.exit('error: fail to run macs2!');    
    ftmp.close();    
    
    # remove the temporary files
    subprocess.check_call(["rm", fout_frag.name]);
    subprocess.check_call(["rm", ftmp.name]);
    return 0

