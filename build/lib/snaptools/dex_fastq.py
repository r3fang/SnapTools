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

import os.path
import sys 
import gzip 
import collections
from snaptools.utilities import file_type
import bz2

def dex_fastq(input_fastq,
              output_fastq,
              index_fastq_list
              ):
    
    """
    De-multiplex fastq files by adding barcode to the beginning of each read name.

    Required:
    --------
    input_fastq: 
        a fastq format file that contains the sequencing reads;

    output_fastq: 
        a fastq file contains output fastq file;

    index_fastq_list: 
        a list of fastq files that contains the barcode

    """

    # check wheather snap file exists
    if not os.path.exists(input_fastq):
        print(('error: ' + input_fastq + ' does not exist!'));
        sys.exit(1);

    if os.path.exists(output_fastq):
        print(('error: ' + output_fastq + ' already exists, remove it first!'));
        sys.exit(1);
    
    for index_fastq in index_fastq_list:
        if not os.path.exists(index_fastq):
            print(('error: ' + index_fastq + ' does not exist!'));
            sys.exit(1);
    
    if file_type(input_fastq) == "gz":
        fr1 = gzip.open(input_fastq, 'rb')
    elif file_type(input_fastq) == "bz2":
        fr1 = bz2.BZ2File(input_fastq, 'r')
    elif file_type(input_fastq) == "txt":
        fr1 = open(input_fastq, 'r')
        
    index_files = []    
    for index_fastq in index_fastq_list:
        if file_type(index_fastq) == "gz":
            fix = gzip.open(index_fastq, 'rb')
        elif file_type(index_fastq) == "bz2":
            fix = bz2.BZ2File(index_fastq, 'r')
        elif file_type(index_fastq) == "txt":
            fix = open(index_fastq, 'r')
        index_files.append(fix)
    

    if output_fastq.endswith("gz"):
        fout = gzip.open(output_fastq, 'wb')
    elif output_fastq.endswith("bz2"):
        fout = bz2.BZ2File(output_fastq, 'w')
    else:
        fout = open(output_fastq, 'w')            
            
    while True:
        cur_r1_name = fr1.readline().strip()[1:]
        if cur_r1_name == "": break        
        cur_r1_read = fr1.readline().strip()
        cur_r1_plus = fr1.readline().strip()
        cur_r1_qual = fr1.readline().strip()
        
        cur_idex_list = []
        for fix in index_files:
            cur_name = fix.readline().strip()[1:]
            cur_read = fix.readline().strip()
            cur_plus = fix.readline().strip()
            cur_qual = fix.readline().strip()
            cur_idex_list.append(cur_read)
        cur_barcode = "".join(cur_idex_list)
            
        if not (cur_name.split()[0] == cur_r1_name.split()[0]): sys.exit("read name does not match")        
        fout.write('@' + cur_barcode + ':' + cur_r1_name+"\n")
        fout.write(cur_r1_read+"\n")
        fout.write("+\n")
        fout.write(cur_r1_qual+"\n")                   

    fout.close()
    fr1.close()
    for fix in index_files:
        fix.close()

