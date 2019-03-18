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
import time
import subprocess
import shlex
import tempfile
import numpy as np
import warnings

try:
    import snaptools.utilities
    import snaptools.global_var
    import snaptools.gtf
    from snaptools.snap import *
    from snaptools.utilities import file_type
except Exception:
    print("Package snaptools not installed!")
    sys.exit(1)

try:
    import numpy as np
except Exception:
    print("Package numpy not installed!")
    sys.exit(1)

try:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        import h5py
    import h5py
except Exception:
    print("Package numpy not installed!")
    sys.exit(1)

try:
    import pysam
except Exception:
    print("Package pysam not installed!")
    sys.exit(1)

try:
    import pybedtools
except Exception:
    print("Package pybedtools not installed!")
    sys.exit(1)


def snap_del(snap_file,
             session_name
             ):

    """
    Delete an existing attribute in a snap file.

    Args:
    --------
    snap_file: 
        a snap format file.

    session_name: 
        attribute to delete ["AM", "GM", "PM", "FM"].
    
    """
    
    if not os.path.exists(snap_file):
        print(('error: ' + snap_file + ' does not exist!'));
        sys.exit(1);
    
    # check if snap_file is a snap-format file
    file_format = snaptools.utilities.checkFileFormat(snap_file);
    if file_format != "snap":
        print(("error: input file %s is not a snap file!" % snap_file));
        sys.exit(1);
    
    fin = h5py.File(snap_file, "r", libver='earliest');
    
    if session_name not in list(fin.keys()):
        print(("error: --session-name %s does not exist in %s" % (session_name, snap_file)));
        sys.exit(1);
    
    fout_name = tempfile.NamedTemporaryFile(delete=False, dir=None);
    fout = h5py.File(fout_name.name, "a", libver='earliest');
    
    session_name_list = list(fin.keys());
    session_name_list.remove(session_name);
    
    for group_name in session_name_list:
        fout.copy(fin[group_name],group_name,shallow=False)
        
    fin.close()
    fout.close()
    subprocess.check_call("\t".join(["mv", fout_name.name, snap_file]), shell=True);


