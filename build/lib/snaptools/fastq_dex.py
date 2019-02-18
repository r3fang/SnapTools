import sys
import os
import subprocess
import collections 
import gzip 
import bz2 
from snaptools.utilities import file_type

def dex_fastq(input_fastq,
              output_fastq,
              index1_fastq,
              index2_fastq,
              index_list):
    """
    Decomplex fastq file by adding barcode to the beginning of the read name.

    Args:
        input_fastq: fastq file contains sequencing reads (demo.R1.fastq.gz), support .fastq, .gz, .bz2 
              
        index1_fastq: fastq file contains r7, i7 barcode (demo.I1.fastq.gz), support .fastq, .gz, .bz2 

        index2_fastq: fastq file contains r5, i5 barcode (demo.I2.fastq.gz), support .fastq, .gz, .bz2 

        index_list: txt file contains pre-designed r7, i7, r5, i5 barcodes (barcodes.txt)
    
    """
    
    # check if those files exist
    if not os.path.isfile(input_fastq): exit("error: \'%s\' not exist" % input_fastq);
    if not os.path.isfile(index1_fastq): exit("error: \'%s\' not exist" % index1_fastq);
    if not os.path.isfile(index2_fastq): exit("error: \'%s\' not exist" % index2_fastq);
    if not os.path.isfile(index_list): exit("error: \'%s\' not exist" % index_list);
    
    if os.path.isfile(output_fastq): exit("error: \'%s\' already exists, remove it first" % output_fastq);
    # check if they are fastq file
    #if not is_fastq(input_fastq): exit("error: \'%s\' is not a fastq file" % input_fastq)
    #if not is_fastq(input_fastq): exit("error: \'%s\' is not a fastq file" % index1_fastq)
    #if not is_fastq(input_fastq): exit("error: \'%s\' is not a fastq file" % index2_fastq)

    # check barcodes
    r7_dict = collections.defaultdict(int)
    i7_dict = collections.defaultdict(int)
    r5_dict = collections.defaultdict(int)
    i5_dict = collections.defaultdict(int)    
    with open(index_list) as fin:
        for line in fin:
            elems = line.split()
            if(len(elems) != 2): continue;
            if elems[1] == "r7":
                r7_dict[elems[0].upper()] = 0
            elif elems[1] == "i7":
                i7_dict[elems[0].upper()] = 0
            elif elems[1] == "r5":
                r5_dict[elems[0].upper()] = 0
            elif elems[1] == "i5":
                i5_dict[elems[0].upper()] = 0
            else:
                exit("error: unorganized index  \'%s\', only support r7, i7, r5, i5" % elems[1])
    
    # check if index is the same length
    if(len(set(map(len, r7_dict.keys()))) != 1):
        exit("error: r7 index has different length");
    
    if(len(set(map(len, i7_dict.keys()))) != 1):
        exit("error: i7 index has different length");
    
    if(len(set(map(len, r5_dict.keys()))) != 1):
        exit("error: r5 index has different length");
    
    if(len(set(map(len, i5_dict.keys()))) != 1):
        exit("error: i5 index has different length");

    r7_len = len(r7_dict.keys()[0]);
    i7_len = len(i7_dict.keys()[0]);
    r5_len = len(r5_dict.keys()[0]);
    i5_len = len(i5_dict.keys()[0]);
    
    if file_type(index1_fastq) == "gz":
        fi1 = gzip.open(index1_fastq, 'rb');
    elif file_type(index1_fastq) == "bz2":
        fi1 = bz2.BZ2File(index1_fastq, 'r');
    elif file_type(index1_fastq) == "txt":
        fi1 = open(index1_fastq, 'r');

    if file_type(index2_fastq) == "gz":
        fi2 = gzip.open(index2_fastq, 'rb');
    elif file_type(index2_fastq) == "bz2":
        fi2 = bz2.BZ2File(index2_fastq, 'r');
    elif file_type(index2_fastq) == "txt":
        fi2 = open(index2_fastq, 'r');

    if file_type(input_fastq) == "gz":
        fr1 = gzip.open(input_fastq, 'rb');
    elif file_type(input_fastq) == "bz2":
        fr1 = bz2.BZ2File(input_fastq, 'r');
    elif file_type(input_fastq) == "txt":
        fr1 = open(input_fastq, 'r');
    
    if output_fastq.endswith("gz"):
        fout = gzip.open(output_fastq, 'wb')
    elif output_fastq.endswith("bz2"):
        fout = bz2.BZ2File(output_fastq, 'w')
    else:
        fout = open(output_fastq, 'w')            
    
    
    TOTAL_READS = 0 # number of totally sequenced reads
    QUALI_READS = 0 # number of usable reads
    
    while True:
        cur_i1_name = fi1.readline().strip()[1:]
        cur_i1_read = fi1.readline().strip()
        cur_i1_plus = fi1.readline().strip()
        cur_i1_qual = fi1.readline().strip()
    
        cur_i2_name = fi2.readline().strip()[1:]
        cur_i2_read = fi2.readline().strip()
        cur_i2_plus = fi2.readline().strip()
        cur_i2_qual = fi2.readline().strip()
    
        cur_r1_name = fr1.readline().strip()[1:]
        cur_r1_read = fr1.readline().strip()
        cur_r1_plus = fr1.readline().strip()
        cur_r1_qual = fr1.readline().strip()
            
        if cur_i1_name == "" or cur_i2_name == "" or cur_r1_name == "": break                
        if not (cur_i1_name.split()[0] == cur_i2_name.split()[0] == cur_r1_name.split()[0]): sys.exit("error: read name does not match")        
        TOTAL_READS += 1
        cur_r7 = cur_i1_read[:r7_len].upper()
        cur_i7 = cur_i1_read[-i7_len:].upper()
        cur_i5 = cur_i2_read[:r5_len].upper()
        cur_r5 = cur_i2_read[-i5_len:].upper()        
        if (cur_i5 in i5_dict) and (cur_i7 in i7_dict) and (cur_r5 in r5_dict) and (cur_r7 in r7_dict):
            QUALI_READS += 1
            cur_barcode_cur = cur_r7 + cur_i7 + cur_i5 + cur_r5
            fout.write('@' + cur_barcode_cur + ':' + cur_r1_name+"\n")
            fout.write(cur_r1_read+"\n")
            fout.write("+\n")
            fout.write(cur_r1_qual+"\n")                   
            r7_dict[cur_r7] += 1
            r5_dict[cur_r5] += 1
            i7_dict[cur_i7] += 1
            i5_dict[cur_i5] += 1                        
    fi1.close()
    fi2.close()
    fr1.close()
    fout.close()
    #### generate a report
    print "Total number of sequencing reads: ", TOTAL_READS
    print "Total number of usable reads: ", QUALI_READS
    print "=========================================="
    for key in r7_dict:
        print '%s\t%s\t%.2f%%' % (key, "r7", round(float(r7_dict[key])/QUALI_READS * 100, 2));
    for key in i7_dict:
        print '%s\t%s\t%.2f%%' % (key, "i7", round(float(i7_dict[key])/QUALI_READS * 100, 2));
    for key in r5_dict:
        print '%s\t%s\t%.2f%%' % (key, "r5", round(float(r5_dict[key])/QUALI_READS * 100, 2));
    for key in i5_dict:
        print '%s\t%s\t%.2f%%' % (key, "i5", round(float(i5_dict[key])/QUALI_READS * 100, 2));
    return 0

