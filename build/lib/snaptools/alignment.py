import sys
import os
import subprocess
import shlex
import tempfile
import pysam
from snaptools.utilities import file_type
import gzip
import bz2
import collections

def count_barcode_cov_from_fastq(fname):
    """
    Count barcode coverage from fastq file

    args:
    -----
	fname: 
        a fastq file, support .gz, .txt, .bz2 file

    output:
    -----
	a dictionary contains barode and its coverage
    """ 
    if file_type(fname) == "gz":
        fin = gzip.open(fname, 'rb');
    elif file_type(fname) == "bz2":
        fin = bz2.BZ2File(fname, 'r');
    elif file_type(fname) == "txt":
        fin = open(fname, 'r');
    else:
	    print "error: unrecoginized fastq file format, only supports .gz, .bz2, .fastq"
	    sys.exit(1)
    
    barcode_cov = collections.defaultdict(lambda : 0)
    while True:
        cur_name = fin.readline().strip()[1:]
        cur_read = fin.readline().strip()
        cur_plus = fin.readline().strip()
        cur_qual = fin.readline().strip()
        if cur_name == "": break                        
        cur_barcode = cur_name.split(":")[0]
        barcode_cov[cur_barcode] += 1
    
    fin.close()
    return(barcode_cov);

def filter_fastq(fname, barcodes, tmp_folder):
    """
    Filter reads belonging to unselected barcodes

    args:
    ------
	fname: 
        a fastq file, support .gz, .txt, .bz2 file

    barcodes:
	    a list contains selected barcodes
    
    tmp_folder:
        folder to store temp file
    
    output:
    ------
    temporary file name
    """

    if file_type(fname) == "gz":
        fin = gzip.open(fname, 'rb');
    elif file_type(fname) == "bz2":
        fin = bz2.BZ2File(fname, 'r');
    elif file_type(fname) == "txt":
        fin = open(fname, 'r');
    else:
        print "error: unrecoginized fastq " + fname +  " file format, only supports .gz, .bz2, .fastq"
        sys.exit(1)    
    
    if len(barcodes) == 0:
	print "error: no barcode is selected"
	sys.exit(1)
    else:
	barcodes = set(barcodes)
    
    fout = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)
    fout_name = fout.name;
    while True:
        cur_name = fin.readline()
        cur_read = fin.readline()
        cur_plus = fin.readline()
        cur_qual = fin.readline()
        if cur_name == "": break                        
        cur_barcode = cur_name.split(":")[0][1:]
        if cur_barcode in barcodes:
            fout.write(cur_name)
            fout.write(cur_read)
            fout.write(cur_plus)
            fout.write(cur_qual)
    fin.close()
    fout.close()
    return(fout_name)
    
def index_ref(input_fasta,
              output_prefix,
              path_to_aligner,
              aligner,
              num_threads):
    """
    Index sequences in the FASTA format 

    Required
    --------
    input_fasta: a fasta file containing the reference genome (mm10.fa)

    Optional
    --------
                
    output_prefix: prefix of output index (optional)
    
    path_to_aligner: directory path access to the aligner

    aligner: aligner name "bwa", "bowtie", "bowtie2" or "minimap2"

    num_threads: number of indexing threads [3];
    """
    
    # if the aligner path given, need to check the existance of the aligner
    if path_to_aligner != None:
        path_to_aligner+="/"
        if not os.path.isdir(path_to_aligner):
            print('Error: path_to_aligner is not a folder')
            sys.exit(1);
        if not os.path.exists(path_to_aligner+aligner):
            print('Error: aligner does not exist')
            sys.exit(1);
    else:
        try:
            # pipe output to /dev/null for silence
            null = open("/dev/null", "w")
            subprocess.Popen(aligner, stdout=null, stderr=null)
            null.close()
        except OSError as e:
            print('Error: ' + aligner + ' does not exist!');
            sys.exit(1);
        path_to_aligner=""
            
    if not os.path.exists(input_fasta):
        print('Error: ' + input_fasta + ' does not exist!');
        sys.exit(1);    
    
    if output_prefix == None:
        output_prefix = os.path.splitext(input_fasta)[0]
    
    aligner = aligner.lower()
    if aligner not in ["bwa", "bowtie", "bowtie2", "minimap2"]:
        print('Error: only support bwa, bowtie, bowtie2, minimap2')
        sys.exit(1);
    # minimap2
    if aligner.lower() == "minimap2":        
        subprocess.check_call([path_to_aligner+"minimap2",
                               "-t",str(num_threads),
                               "-d",output_prefix+".mmi",
                               input_fasta])        
    # bowtie2
    if aligner.lower() == "bowtie2":
        subprocess.check_call([path_to_aligner + "bowtie2-build", 
                               "-f", input_fasta, 
                               output_prefix])
    # bowtie
    if aligner.lower() == "bowtie":
        subprocess.check_call([path_to_aligner + "bowtie-build", 
                               "-f", input_fasta, 
                               output_prefix])
    # bwa
    if aligner.lower() == "bwa":
        subprocess.check_call([path_to_aligner+"bwa", 
                               "index", input_fasta])
    return 0

def run_align_pe(input_reference,
                 input_fastq1,
                 input_fastq2,
                 output_bam,
                 aligner,
                 path_to_aligner,
                 read_fastq_command,
                 num_threads,
		         min_cov,
                 aligner_options,
                 if_sort,
                 tmp_folder,
                 overwrite,
                 verbose):

    """
    Map paired-end single-cell ATAC-seq reads

    args:
    -----
    input_reference: 
        reference genome file generated by index_reference

    input_fastq1: 
        a fastq file contains R1 reads, supports .fq, .fastq, .gz, .bz2

    input_fastq2: 
        a fastq file contains R2 reads, supports .fq, .fastq, .gz, .bz2

    output_bam: 
        a master bam file contains alignments
    
    path_to_aligner: 
        directory path access to the aligner

    aligner: 
        aligner name "bwa", "bowtie", "bowtie2" or "minimap2"
    
    aligner_options:
        a list of strings indicating options you'd like passed to aligner.
        (default for bowtie2)
        (default for bwa: "mem")
        (default for minimap2: "-ax sr")
        (default for bowtie: "-S -k 2 --best --strata --chunkmbs 64 -n 1")
                 
    num_threads: 
        number of mapping threads [1];

    if_sort: 
        if sort the alignment based on read name [True];

    tmp_folder: 
        where to store the temporary files [None];
    
    min_cov: 
        barcodes of fragments fewer than min_cov will be filtered before alingment; 

    read_fastq_command: 
        command to uncompress a compressed fastq file i.e. 'zcat', 'bzcat' [None];

    overwrite:
        whether to overwrite the output file if it already exists [False];
                 
    verbose:
        a boolen variable indicates whether to output the progress [True];         
    """
    # if the aligner path given, need to check the existance of the aligner
    if path_to_aligner != None:
        path_to_aligner+="/"
        if not os.path.isdir(path_to_aligner):
            print('Error: path_to_aligner is not a folder')
            sys.exit(1);
        if not os.path.exists(path_to_aligner+aligner):
            print('Error: aligner does not exist')
            sys.exit(1);
    else:
        try:
            # pipe output to /dev/null for silence
            null = open("/dev/null", "w")
            subprocess.Popen(aligner, stdout=null, stderr=null)
            null.close()
        except OSError as e:
            print('Error: ' + aligner + ' does not exist!');
            sys.exit(1);
        path_to_aligner=""
        
    if(tmp_folder!=None):
        if not os.path.isdir(tmp_folder):
            print('Error: tmp_folder is not a folder or does not exist')
            sys.exit(1);
    
    # check the existance of input and output files
    if not os.path.exists(input_fastq1):
        sys.exit('Error: ' + input_fastq1 + ' does not exist!');
    if not os.path.exists(input_fastq2):
        sys.exit('Error: ' + input_fastq2 + ' does not exist!');
    
    if os.path.isfile(output_bam): 
        if overwrite:
            subprocess.check_call(["rm", output_bam]);
        else:
            sys.exit("error: \'%s\' already exists, remove it first" % output_bam);

    if input_fastq1 == input_fastq2:
        sys.exit("error: --input_fastq1 and --input_fastq2 are same file");
    
    # check if can create the output_bam file
    try:
        with open(output_bam, "w") as outfile:
            outfile.write('Hello World')
        subprocess.check_call(["rm", output_bam]);
    except IOError:
        print "error: could not create %s, check if the folder exists." % output_bam;
        sys.exit(1)
        
    if min_cov > 0:
        barcode_dict = count_barcode_cov_from_fastq(input_fastq1);
        barcode_sel = set([key for key in barcode_dict if barcode_dict[key] > min_cov]);
        if len(barcode_sel) == 0:
            print "error: no barcode contains fragments more than --min-cov, lower --min-cov and try it again!"
            sys.exit(1)
        input_fastq1 = filter_fastq(input_fastq1, barcode_sel, tmp_folder);
        input_fastq2 = filter_fastq(input_fastq2, barcode_sel, tmp_folder);
        read_fastq_command = "cat"
        
    # check validity of aligner
    aligner = aligner.lower()
    if aligner not in ["bwa", "bowtie", "bowtie2", "minimap2"]:
        sys.exit('Error: only support bwa, bowtie, bowtie2, minimap2');
    
    # default aligner option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-X 1000", "-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]
            aligner_options.append("--phred33-quals")
        elif aligner.lower() == "bowtie2": # bowtie2
            aligner_options = []
            aligner_options.append("--phred33-quals")
        elif aligner.lower() == "bwa": # bowtie2
            aligner_options = ["mem"]
    options = aligner_options
    
    # update num_threads if it is given in the aligner_options
    if aligner in ["bowtie", "bowtie2"]:
        if " ".join(options).find(" -p ") == -1:
            options.append("-p "+str(num_threads))
    elif aligner in ["minimap2", "bwa"]:
        if " ".join(options).find(" -t ") == -1:
            options.append("-t "+str(num_threads))
    else:
        sys.exit('Error: only support bwa, bowtie, bowtie2, minimap2');
        
    # if cat_cmd is not given, automatically detect file type and choose cat_cmd
    if read_fastq_command == None:
        if file_type(input_fastq1) == "gz":
            read_fastq_command = "zcat"
        elif file_type(input_fastq1) == "bz2":
            read_fastq_command = "bzcat"
        elif file_type(input_fastq1) == "txt": # .fq or fastq file
            read_fastq_command = "cat"
        else:
            sys.exit('Error: unrecoganized fastq file, supports .fq, .fastq, .gz, .bz2 file');

    # mapping and write the alignments into a temporary file        
    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append(input_reference)
        args.append("<(" + read_fastq_command + " " + input_fastq1  + ")")
        args.append("<(" + read_fastq_command + " " + input_fastq2  + ")")
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append(input_reference)
        args.append("-1 " + "<(" + read_fastq_command + " " + input_fastq1  + ")")
        args.append("-2 " + "<(" + read_fastq_command + " " + input_fastq2  + ")")
    elif aligner.lower() == "bowtie2": # bowtie2
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("-x " + input_reference)        
        args.append("-1 " + "<(" + read_fastq_command + " " + input_fastq1  + ")")
        args.append("-2 " + "<(" + read_fastq_command + " " + input_fastq2  + ")")
    else:
        args = [path_to_aligner+"bwa"]
        args.extend(options)
        args.append(input_reference)
        args.append("<(" + read_fastq_command + " " + input_fastq1  + ")")
        args.append("<(" + read_fastq_command + " " + input_fastq2  + ")")        
    
    ftmp = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)
    try:
        subprocess.check_call(" ".join(args), stdout=ftmp, shell=True, executable='/bin/bash');
    except subprocess.CalledProcessError as e:
        sys.exit('error: fail to run alignment, check if aligner and reference genome is correct!');    
    ftmp.close();    

    if(if_sort):    
        pysam.sort("-n", "-@", str(num_threads), "-o", output_bam, ftmp.name);
    else:
        samfile = pysam.AlignmentFile(ftmp.name, "r")
        fout = pysam.AlignmentFile(output_bam, "wb", template=samfile)
        for read in samfile.fetch():
            fout.write(read)
        fout.close()
        samfile.close()        
    subprocess.check_call(["rm", ftmp.name]);
    
    # remove tmp fastq file after alignment
    if min_cov > 0:
        subprocess.check_call(["rm", input_fastq1]);
        subprocess.check_call(["rm", input_fastq2]);
    return 0

def run_align_se(input_reference,
                 input_fastq1,
                 output_bam,
                 aligner="bwa",
                 path_to_aligner=None,
                 read_fastq_command=None,
                 num_threads=3,
                 aligner_options=None,
                 if_sort=True,
                 tmp_folder=None,
                 overwrite=False):

    """
    Map single-cell ATAC-seq reads in single-end mode

    Required
    --------
    input_reference: reference genome file generated by index_reference

    input_fastq1: a fastq file contains R1 reads, supports .fq, .fastq, .gz, .bz2

    output_bam: a bam file contains alignments
    
    Optional
    --------
    path_to_aligner: directory path access to the aligner

    aligner: aligner name "bwa", "bowtie", "bowtie2" or "minimap2"

    aligner_options is a list of strings indicating options you'd like passed to aligner.
        (default for bowtie2: "-X 1000 -k 2 --no-mixed --no-discordant")
        (default for bowtie: "-X 1000 -S -k 1 -m 1 --best --strata --chunkmbs 64 -n 1")
        (default for bwa: "mem")
        (default for minimap2: "-ax sr --secondary=no")
                 
    num_threads: number of mapping threads [3];

    if_sort: if sort the alignment based on read name [True];

    tmp_folder: where to store the temporary files [None];
    
    read_fastq_command: command to uncompress a compressed fastq file i.e. 'zcat', 'bzcat' [None];

    overwrite: whether to overwrite the output file if it already exists [False];
    """
    # if the aligner path given, need to check the existance of the aligner
    if path_to_aligner != None:
        path_to_aligner+="/"
        if not os.path.isdir(path_to_aligner):
            print('Error: path_to_aligner is not a folder')
            sys.exit(1);
        if not os.path.exists(path_to_aligner+aligner):
            print('Error: aligner does not exist')
            sys.exit(1);
    else:
        try:
            # pipe output to /dev/null for silence
            null = open("/dev/null", "w")
            subprocess.Popen(aligner, stdout=null, stderr=null)
            null.close()
        except OSError as e:
            print('Error: ' + aligner + ' does not exist!');
            sys.exit(1);
        path_to_aligner=""
        
    if(tmp_folder!=None):
        if not os.path.isdir(tmp_folder):
            print('Error: tmp_folder is not a folder or does not exist')
            sys.exit(1);
    
    # check the existance of input and output files
    if not os.path.exists(input_fastq1):
        sys.exit('Error: ' + input_fastq1 + ' does not exist!');
    
    if os.path.isfile(output_bam): 
        if overwrite:
            subprocess.check_call(["rm", output_bam]);
        else:
            sys.exit("error: \'%s\' already exists, remove it first" % output_bam);
    
    # check if can create the output_bam file
    try:
        with open(output_bam, "w") as outfile:
            outfile.write('')
        subprocess.check_call(["rm", output_bam]);
    except IOError:
        print "error: could not create %s, check if the folder exists." % output_bam;
        sys.exit(1)
    
    # check validity of aligner
    aligner = aligner.lower()
    if aligner not in ["bwa", "bowtie", "bowtie2", "minimap2"]:
        sys.exit('Error: only support bwa, bowtie, bowtie2, minimap2');
    
    # default aligner option
    if aligner_options is None:
        if aligner.lower() == "minimap2":
            aligner_options = ["-ax","sr","--secondary=no"]
        elif aligner.lower() == "bowtie":
            aligner_options = ["-S", "-k 1", "-m 1", "--best", "--strata",
                               "--chunkmbs 3072", "-n 1", "-e 100"]
            aligner_options.append("--phred33-quals")
        elif aligner.lower() == "bowtie2": # bowtie2
            aligner_options = []
            aligner_options.append("--phred33-quals")
        elif aligner.lower() == "bwa": # bowtie2
            aligner_options = ["mem"]
    options = aligner_options
    
    # update num_threads if it is given in the aligner_options
    if aligner in ["bowtie", "bowtie2"]:
        if " ".join(options).find(" -p ") == -1:
            options.append("-p "+str(num_threads))
    elif aligner in ["minimap2", "bwa"]:
        if " ".join(options).find(" -t ") == -1:
            options.append("-t "+str(num_threads))
    else:
        sys.exit('Error: only support bwa, bowtie, bowtie2, minimap2');
        
    # if cat_cmd is not given, automatically detect file type and choose cat_cmd
    if read_fastq_command == None:
        if file_type(input_fastq1) == "gz":
            read_fastq_command = "zcat"
        elif file_type(input_fastq1) == "bz2":
            read_fastq_command = "bzcat"
        elif file_type(input_fastq1) == "txt": # .fq or fastq file
            read_fastq_command = "cat"
        else:
            sys.exit('Error: unrecoganized fastq file, supports .fq, .fastq, .gz, .bz2 file');

    # mapping and write the alignments into a temporary file        
    if aligner.lower() == "minimap2":
        args = [path_to_aligner+"minimap2"]
        args.extend(options)
        args.append(input_reference)
        args.append("<(" + read_fastq_command + " " + input_fastq1  + ")")
    elif aligner.lower() == "bowtie":
        args = [path_to_aligner+"bowtie"]
        args.extend(options)
        args.append(input_reference)
        args.append("-1 " + "<(" + read_fastq_command + " " + input_fastq1  + ")")
    elif aligner.lower() == "bowtie2": # bowtie2
        args = [path_to_aligner+"bowtie2"]
        args.extend(options)
        args.append("-x " + input_reference)        
        args.append("-1 " + "<(" + read_fastq_command + " " + input_fastq1  + ")")
    else:
        args = [path_to_aligner+"bwa"]
        args.extend(options)
        args.append(input_reference)
        args.append("<(" + read_fastq_command + " " + input_fastq1  + ")")
    
    ftmp = tempfile.NamedTemporaryFile(delete=False, dir=tmp_folder)
    try:
        subprocess.check_call(" ".join(args), stdout=ftmp, shell=True, executable='/bin/bash');
    except subprocess.CalledProcessError as e:
        sys.exit('Error: failed to run alignment, check if aligner and reference genome is correct!');    
    ftmp.close();    

    if(if_sort):    
        pysam.sort("-n", "-@", str(num_threads), "-o", output_bam, ftmp.name);
    else:
        samfile = pysam.AlignmentFile(ftmp.name, "r")
        fout = pysam.AlignmentFile(output_bam, "wb", template=samfile)
        for read in samfile.fetch():
            fout.write(read)
        fout.close()
        samfile.close()        
    subprocess.check_call(["rm", ftmp.name]);
            
    return 0
