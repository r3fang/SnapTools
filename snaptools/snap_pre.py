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

def snap_pre(input_file,
             output_snap,
             genome_name,
             genome_size,
             min_mapq,
             min_flen,
             max_flen,
             min_cov,
             barcode_file,
             keep_chrm,
             keep_single,
             keep_secondary,
             keep_discordant,
             tmp_folder,
             overwrite,
             qc_file,
             verbose):

    """
    Pre-processing to create a snap file from a bam or bed that contains alignments or a bed file that contains fragments.

    Args:
    --------
    input_file: 
        a bed or bam file contains fragments 
        
    output_snap: 
        output snap file.
    
    Optional
    --------
    min_mapq: 
            min mappability [10]. fragments with mappability less than 10 will be filtered
               
    min_flen: 
            min fragment size [0]. fragments of length shorter than min_flen will be filtered

    max_flen: 
            max fragment size [1000]. fragments of length bigger than min_flen will be filtered

    min_cov:
            min coverage per barcode. barcodes with sequencing fragments less than min_cov will be filtered before processed  
    
    barcode_file: 
            a txt file that contains selected barcodes to create count matrix [None]

    keep_chrm:
            a boolen variable indicates whether to keep reads mapped to chrM [True]

    keep_single:
            a boolen variable indicates whether to keep single-end reads [False]

    keep_secondary:
            a boolen variable indicates whether to keep secondary alingments [False]

    keep_discordant:
            a boolen variable indicates whether to keep discordant alingments [False]
            
    tmp_folder: 
            where to store the temporary files [None];
    
    qc_file: 
            a boolen variable indicates whether to create a master qc file [True];

    verbose: 
            a boolen variable indicates whether to output the progress [True];
    """
    
    if not os.path.exists(input_file):
        print 'error: ' + input_file + ' does not exist!';
        sys.exit(1);

    if min_flen < snaptools.global_var.MIN_FRAGMENT_LEN:
        print('error: --min-flen can not be smaller than %s ' % str(snaptools.global_var.MIN_FRAGMENT_LEN));
        sys.exit(1); 
        
    if max_flen > snaptools.global_var.MAX_FRAGMENT_LEN:
        print('error: --max-flen can not be smaller than %s ' % str(snaptools.global_var.MAX_FRAGMENT_LEN));
        sys.exit(1); 

    if barcode_file != None:
        if not os.path.exists(barcode_file):
            print('error: --barcode-file %s does not exist!' % barcode_file);
            sys.exit(1);
    
    if os.path.exists(output_snap):
        if overwrite == True:
            subprocess.check_call(["rm", output_snap]);
        else:
            print ('error: %s already exists, change --overwrite!' % output_snap);
            sys.exit(1);    
    
    if(tmp_folder!=None):
        if(not os.path.isdir(tmp_folder)):
            print('error: %s does not exist!' % tmp_folder);
            sys.exit(1);

    # automatically determine the file format
    file_format = snaptools.utilities.checkFileFormat(input_file);
    
    if file_format != "bed" and file_format != "bam":
        print ('error: unrecognized input format %s, only supports bed and bam file!' % input_file)
        sys.exit(1);    
    
    if not os.path.exists(genome_size):
        print 'error: ' + genome_size + ' does not exist!'
        sys.exit(1);
    
    if genome_name not in snaptools.global_var.GENOMELIST:
        print('error: --genome-name unrecoginized genome identifier %s' % str(genome_name));
        print("Do you mean: %s ?" % ",".join(snaptools.utilities.minEditDistanceString(genome_name, snaptools.global_var.GENOMELIST)))
        sys.exit(1)
            
    # read the reference genome info from given genome size txt file
    genome_dict = snaptools.utilities.readGenomeSizeFromTxt(genome_size);
        
    # enforce the following parameters if it is a bed file, assuming the bed file has been filtered
    if file_format == "bed":
        keep_single=False;
        keep_secondary=False;
        keep_discordant=False;
        min_mapq=0;
    
    # check if a file is sorted by read name or queryname
    if not snaptools.utilities.isQuerynameSorted(input_file, file_format):
        print ('error: %s is not a sorted by read name!' % input_file)
        sys.exit(1);    
    
    # create the header session
    header = collections.defaultdict();
    header["MG"] = snaptools.global_var.MAGIC_STRING;
    header["DT"] = str(datetime.datetime.now());
    header["VN"] = snaptools.__version__;
    header["CL"] = "\t".join(map(str, sys.argv));
    header["CW"] = os.getcwd();

    # bed file does not contain any alignment information, force it to "NA"
    align_dict = snaptools.utilities.readAlignmentInfo(input_file, file_format);

    if barcode_file != None:
        barcode_dict = snaptools.snap.getBarcodesFromTxt(barcode_file);
    else:
        barcode_dict = snaptools.snap.getBarcodesFromInput(input_file, file_format);

    if min_cov > 0:
        barcode_cov = snaptools.snap.getBarcodeCov(barcode_dict.keys(), input_file, file_format);
        barcodes = [key for key in barcode_cov if barcode_cov[key] > min_cov];
        barcode_dict = collections.OrderedDict(zip(sorted(barcodes), range(1,(len(barcodes)+1))));        
    
    f = h5py.File(output_snap, "w", libver='earliest');
    f.create_group("HD")
    f.create_dataset("HD/MG", data=header["MG"]);
    f.create_dataset("HD/VN", data=header["VN"]);
    f.create_dataset("HD/DT", data=header["DT"]);
    f.create_dataset("HD/CL", data=header["CL"]);
    f.create_dataset("HD/CW", data=header["CW"]);
    
    # header/alignment
    f.create_dataset("HD/AL/PN", data=align_dict["PN"]);
    f.create_dataset("HD/AL/ID", data=align_dict["ID"]);
    f.create_dataset("HD/AL/VN", data=align_dict["VN"]);
    f.create_dataset("HD/AL/CL", data=align_dict["CL"]);
    
    fragment_chrom_ds = f.create_dataset("FM/fragChrom", (0,), 
                     dtype=h5py.special_dtype(vlen=bytes), 
                     maxshape=(None,), chunks=(10**4,), 
                     compression="gzip", compression_opts=9);        

    fragment_start_ds = f.create_dataset("FM/fragStart", (0,), 
                     dtype="uint32", 
                     maxshape=(None,), chunks=(10**4,),
                     compression="gzip", compression_opts=9);        

    fragment_len_ds = f.create_dataset("FM/fragLen", (0,), 
                     dtype="uint16", 
                     maxshape=(None,), chunks=(10**4,),
                     compression="gzip", compression_opts=9);        

    fragment_chrom_ds[:] = [];
    fragment_start_ds[:] = [];
    fragment_len_ds[:] = [];
    
    check_point = 0;
    start_time = time.time();
    qc_dict = collections.defaultdict(qc);   
    num_barcode = len(barcode_dict);
    
    for read_group in snaptools.snap.group_reads_by_barcode(input_file, file_format):
        frag_list = [];    
        if file_format == "bam":
            for (read1, read2, is_paired, is_secondary) in pairReadsByName(read_group):   
                if is_paired:
                    frag = snaptools.snap.readPairToFragment(read1, read2, is_secondary);
                else: # supplementary alignments or unmated reads;
                    frag = snaptools.snap.readToFragment(read1, is_secondary);            
                # extract the barcode
                barcode = frag.qname.split(":")[0].upper();
                # only for printing the progress
                check_point += 1
                if verbose and check_point%100000 == 0:
                    print "%d\ttags, %s seconds " % (check_point, time.time() - start_time);  
                # skip this barcode if it is not in the barcode list 
                if barcode not in barcode_dict: break
                # total number of sequencing fragments (exclude supplementary alignments)
                if frag.is_secondary == False:
                    qc_dict[barcode].total += 1;
                ## 1. Filter non-uniquely mapped fragments
                if frag.mapq < min_mapq: continue;            
                qc_dict[barcode].mapped += 1;
                # 2. if it is single read, keep it only if keep_unpaired is true
                #    if it is paired, keep it only if it is approperly paired
                if frag.is_single == True: 
                    if frag.is_secondary:
                        qc_dict[barcode].secondary += 1
                        if not keep_secondary: continue;
                    else:
                        qc_dict[barcode].single += 1;
                        if not keep_single: continue;                    
                else: # paired-end reads
                    qc_dict[barcode].paired += 1;
                    if frag.is_proper_pair:
                        qc_dict[barcode].proper_paired += 1;
                    else:
                        if not keep_discordant: continue;    
                # 3. check fragment size
                if frag.flen > min_flen and frag.flen < max_flen:
                    qc_dict[barcode].proper_flen += 1;
                else:
                    continue
                # 4. combine single and paired as fragments
                frag_list.append((frag.chrom, frag.pos, frag.pos+frag.flen, barcode));        
        else:
            for _read in read_group:
                elems = _read.split();
                frag = fragment("NA", elems[0], int(elems[1]), int(elems[2]) - int(elems[1]), 60, False, False, True)
                # extract the barcode
                barcode = elems[3].split(":")[0].upper();
                # only for printing the progress
                check_point += 1;
                if verbose and check_point%100000 == 0:
                    print "%d\ttags, %s seconds " % (check_point, time.time() - start_time);  
                # skip this barcode if it is not in the barcode list 
                if barcode not in barcode_dict: break;
                # total number of sequencing fragments (exclude supplementary alignments)
                if frag.is_secondary == False:
                    qc_dict[barcode].total += 1;
                ## 1. Filter non-uniquely mapped fragments
                if frag.mapq < min_mapq: continue;            
                qc_dict[barcode].mapped += 1;
                # 2. if it is single read, keep it only if keep_unpaired is true
                #    if it is paired, keep it only if it is approperly paired
                if frag.is_single == True: 
                    if frag.is_secondary:
                        qcDict[barcode].secondary += 1
                        if not keep_secondary: continue;
                    else:
                        qc_dict[barcode].single += 1;
                        if not keep_single: continue;                    
                else: # paired-end reads
                    qc_dict[barcode].paired += 1;
                    if frag.is_proper_pair:
                        qc_dict[barcode].proper_paired += 1;
                    else:
                        if not keep_discordant: continue;    

                # 3. check fragment size
                if frag.flen > min_flen and frag.flen < max_flen:
                    qc_dict[barcode].proper_flen += 1;
                else:
                    continue
                # 4. combine single and paired as fragments
                frag_list.append((frag.chrom, frag.pos, frag.pos+frag.flen, barcode));        
        
        if barcode not in barcode_dict: continue  
        # 5. remove duplicate fragments
        qc_dict[barcode].usable = len(frag_list);
        frag_list_uniq = set(frag_list); # remove duplicated fragments
        qc_dict[barcode].uniq = len(frag_list_uniq);
                
        ## 6. keep mt reads if keep_chrm is true
        if(len(frag_list_uniq) == 0): 
            qc_dict[barcode].chrM = 0
            continue

        for item in frag_list_uniq:
            if item[0] == "chrM":
                qc_dict[barcode].chrM += 1                   
        if not keep_chrm:
            frag_list_uniq = [item for item in frag_list_uniq if item[0] != "chrM"];            

        # just in case the only fragments are chrM
        qc_dict[barcode].final = len(frag_list_uniq);
        if len(frag_list_uniq) == 0: continue;
        
        fragment_chrom_ds.resize(fragment_chrom_ds.shape[0]+len(frag_list_uniq), axis=0);
        fragment_start_ds.resize(fragment_start_ds.shape[0]+len(frag_list_uniq), axis=0);
        fragment_len_ds.resize(fragment_len_ds.shape[0]+len(frag_list_uniq), axis=0);
        
        fragment_chrom_ds[-len(frag_list_uniq):] = [item[0] for item in frag_list_uniq];
        fragment_start_ds[-len(frag_list_uniq):] = [item[1] for item in frag_list_uniq];
        fragment_len_ds[-len(frag_list_uniq):] = [item[2] - item[1] for item in frag_list_uniq];

        del frag_list, frag_list_uniq;
    
    # genome information
    f.create_dataset("HD/SQ/ID", data = genome_name);
    f.create_dataset("HD/SQ/SN", data = genome_dict.keys());
    f.create_dataset("HD/SQ/SL", data = genome_dict.values());

    # barcode
    f.create_dataset("BD/name",data=(barcode_dict.keys()), compression="gzip", compression_opts=9);
    f.create_dataset("BD/TN", data=[qc_dict[key].total for key in barcode_dict],  dtype="uint32");
    f.create_dataset("BD/UM", data=[qc_dict[key].mapped for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/SE", data=[qc_dict[key].single for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/SA", data=[qc_dict[key].secondary for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/PE", data=[qc_dict[key].paired for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/PP", data=[qc_dict[key].proper_paired for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/PL", data=[qc_dict[key].proper_flen for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/US", data=[qc_dict[key].usable for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/UQ", data=[qc_dict[key].uniq for key in barcode_dict], dtype="uint32");
    f.create_dataset("BD/CM", data=[qc_dict[key].chrM for key in barcode_dict], dtype="uint32");            
    barcode_frag_num = [qc_dict[barcode].final for barcode in barcode_dict];
    barcode_frag_start = np.cumsum(barcode_frag_num) - barcode_frag_num;
    f.create_dataset("FM/barcodePos", data=barcode_frag_start + 1, dtype="uint32");        
    f.create_dataset("FM/barcodeLen", data=barcode_frag_num, dtype="uint32");        
    f.close()  

    if qc_file:
        with open(output_snap+".qc", "w") as fout:
            fout.write("Total number of unique barcodes:             %d\n"%num_barcode)
            fout.write("TN - Total number of fragments:              %d\n"%sum([qc_dict[key].total for key in qc_dict]))
            fout.write("UM - Total number of uniquely mapped:        %d\n"%sum([qc_dict[key].mapped for key in qc_dict]))
            fout.write("SE - Total number of single ends:            %d\n"%sum([qc_dict[key].single for key in qc_dict]))
            fout.write("SA - Total number of secondary alignments:   %d\n"%sum([qc_dict[key].secondary for key in qc_dict]))
            fout.write("PE - Total number of paired ends:            %d\n"%sum([qc_dict[key].paired for key in qc_dict]))
            fout.write("PP - Total number of proper paired:          %d\n"%sum([qc_dict[key].proper_paired for key in qc_dict]))
            fout.write("PL - Total number of proper frag len:        %d\n"%sum([qc_dict[key].proper_flen for key in qc_dict]))
            fout.write("US - Total number of usable fragments:       %d\n"%sum([qc_dict[key].usable for key in qc_dict]))
            fout.write("UQ - Total number of unique fragments:       %d\n"%sum([qc_dict[key].uniq for key in qc_dict]))
            fout.write("CM - Total number of chrM fragments:         %d\n"%sum([qc_dict[key].chrM for key in qc_dict]))
    return 0

