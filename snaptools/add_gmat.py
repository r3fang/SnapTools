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

def snap_gmat(snap_file,
              gene_file,
              buffer_size,
              tmp_folder,
              verbose):
              
    """
    Create a cell x peak matrix from snap file.

    Required:
    --------
    snap_file: 
        a snap format file;

    gene_file: 
        a bed format file contains gene locations;
    
    Optional:
    --------
    tmp_folder:
        folder to store temporarily created files;

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
    
    if not os.path.exists(gene_file):
        print 'error: ' + gene_file + ' does not exist!'
        sys.exit(1);
    
    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(gene_file);
    if file_format != "bed":
        print "error: input file %s is not a bed file!" % snap_file;

    # check if gene file is a bed file with the last column as gene
    gene_list = set([str(item.name) for item in pybedtools.BedTool(gene_file)]);    
    gene_dict = collections.OrderedDict(zip(gene_list, range(1, len(gene_list) + 1)));
    
    
    # extract the barcodes
    barcode_dict = getBarcodesFromSnap(snap_file);

    # check if GM session already exists
    f = h5py.File(snap_file, "r", libver='earliest');
    if "GM" in f:
        print "error: cell x peak matrix already exists";
        sys.exit(1);
    f.close();
    
        
    # first cut the fragments into small piecies, write them down
    fout_frag = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder);
    dump_read(snap_file, fout_frag.name, buffer_size, None, tmp_folder, True);

    # in parallel find the overlap cell and peaks
    frag_bt = pybedtools.BedTool(fout_frag.name); 
    peak_bt = pybedtools.BedTool(gene_file);      
    
    # count for frequency            
    cell_peak_arr = collections.defaultdict(list);            
    for item in frag_bt.intersect(peak_bt, wa=True, wb=True):
        key = str(item.fields[7]);
        if key in gene_dict:
            idy = gene_dict[key];
            barcode = item.name.split(":")[0];
            idx = barcode_dict[barcode].id;
            cell_peak_arr[idx].append(idy);        
    
    IDX_LIST = [];
    IDY_LIST = [];
    VAL_LIST = [];
    for barcode_id in cell_peak_arr:
        d = collections.Counter(cell_peak_arr[barcode_id]);
        IDX_LIST += [barcode_id] * len(d);
        for peak_id in d:
            IDY_LIST.append(peak_id);
            VAL_LIST.append(d[peak_id]);
        
    f = h5py.File(snap_file, "a", libver='earliest');
    dt = h5py.special_dtype(vlen=bytes);
    f.create_dataset("GM/name", data=gene_dict.keys(), dtype=h5py.special_dtype(vlen=bytes), compression="gzip", compression_opts=9);    
    f.create_dataset("GM/idx",       data=IDX_LIST,      dtype="uint32", compression="gzip", compression_opts=9);
    f.create_dataset("GM/idy",       data=IDY_LIST,      dtype="uint32", compression="gzip", compression_opts=9);
    f.create_dataset("GM/count",     data=VAL_LIST,    dtype="uint8", compression="gzip", compression_opts=9);     
    f.close()  
    # remove the temporary files
    subprocess.check_call(["rm", fout_frag.name]);
    return 0

