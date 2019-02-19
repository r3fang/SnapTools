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

import sys
import argparse
import snaptools
from snaptools.utilities import str2bool

def parse_args():
    # create the top-level parser
    parser = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description = "Program: snaptools (A module for working with snap files in Python)\n"
         + "Version: " + snaptools.__version__ + "\n"
         + "Contact: Rongxin Fang" + "\n"
         + "E-mail:  r4fang@gmail.com"
    )

    # create the sub-level parser
    subparsers = parser.add_subparsers(
         title="functions",
         dest="command",
         metavar="")

    #add_fastq_dex_subparser(subparsers)
    add_index_ref_subparser(subparsers);
    add_align_pe_subparser(subparsers);
    add_align_se_subparser(subparsers);
    add_snap_pre_subparser(subparsers);
    add_snap_bmat_subparser(subparsers);
    add_snap_pmat_subparser(subparsers);
    add_snap_gmat_subparser(subparsers);
    add_dump_read_subparser(subparsers);
    add_dump_barcode_subparser(subparsers);
    add_call_peak_subparser(subparsers);
    add_louvain_subparser(subparsers);
    
    if len(sys.argv) > 1:
         ## print out version
         if (sys.argv[1] == '--version' or sys.argv[1] == '-v'):
              print(snaptools.__version__)
              exit()
         ## all functions
         args = parser.parse_args()
    else:
         args = parser.parse_args(["-h"])
         exit()

    if args.command == "index-genome":
         from snaptools.alignment import index_ref
         index_ref(input_fasta=args.input_fasta,
                   output_prefix=args.output_prefix,
                   aligner=args.aligner,
                   path_to_aligner=args.path_to_aligner,
                   num_threads=args.num_threads)
                   
    if args.command == "align-paired-end":
         from snaptools.alignment import run_align_pe
         run_align_pe(input_reference=args.input_reference,
                      input_fastq1=args.input_fastq1,
                      input_fastq2=args.input_fastq2,
                      output_bam=args.output_bam,
                      aligner=args.aligner,
                      path_to_aligner=args.path_to_aligner,
                      num_threads=args.num_threads,
                      aligner_options=args.aligner_options,
                      if_sort=args.if_sort,
		              min_cov=args.min_cov,
                      read_fastq_command=args.read_fastq_command,
                      tmp_folder=args.tmp_folder,
                      overwrite=args.overwrite,
                      verbose=args.verbose)
        
    if args.command == "align-single-end":
         from snaptools.alignment import run_align_se
         run_align_se(input_reference=args.input_reference,
                      input_fastq1=args.input_fastq1,
                      output_bam=args.output_bam,
                      aligner=args.aligner,
                      path_to_aligner=args.path_to_aligner,
                      num_threads=args.num_threads,
                      aligner_options=args.aligner_options,
                      if_sort=args.if_sort,
		              min_cov=args.min_cov,
                      read_fastq_command=args.read_fastq_command,
                      tmp_folder=args.tmp_folder,
                      overwrite=args.overwrite)
    
    if args.command == "snap-view-head":
        from snaptools.snap import run_snap_view_head
        run_snap_view_head(input_snap=args.input_snap)
        
    if args.command == "snap-pre":
         from snaptools.snap_pre import snap_pre
         snap_pre(input_file=args.input_file,
                  output_snap=args.output_snap,
                  genome_name=args.genome_name,
                  genome_size=args.genome_size,
                  min_mapq=args.min_mapq,
                  min_flen=args.min_flen,
                  max_flen=args.max_flen,
                  min_cov=args.min_cov,
                  barcode_file=args.barcode_file,
                  keep_chrm=args.keep_chrm,
                  keep_single=args.keep_single,
                  keep_secondary=args.keep_secondary,
                  keep_discordant=args.keep_discordant,
                  tmp_folder=args.tmp_folder,
                  overwrite=args.overwrite,
                  qc_file=args.qc_file,
                  verbose=args.verbose)

    if args.command == "dump-fragment":
         from snaptools.snap import dump_read
         dump_read(snap_file=args.snap_file,
                   output_file=args.output_file,
                   buffer_size=args.buffer_size,
                   barcode_file=args.barcode_file,
                   tmp_folder=args.tmp_folder,
                   overwrite=args.overwrite)

    if args.command == "dump-barcode":
         from snaptools.dump_barcode import dump_barcode
         dump_barcode(snap_file=args.snap_file,
                      output_file=args.output_file,
                      barcode_file=args.barcode_file,
                      tmp_folder=args.tmp_folder,
                      overwrite=args.overwrite)

    if args.command == "snap-add-bmat":
         from snaptools.add_bmat import snap_bmat
         snap_bmat(snap_file=args.snap_file,
                   bin_size_list=args.bin_size_list,
                   tmp_folder=args.tmp_folder,
                   verbose=args.verbose)
    
    if args.command == "snap-add-pmat":
         from snaptools.add_pmat import snap_pmat
         snap_pmat(snap_file=args.snap_file,
                   peak_file=args.peak_file,
                   tmp_folder=args.tmp_folder,
                   verbose=args.verbose)

    if args.command == "snap-add-gmat":
         from snaptools.add_gmat import snap_gmat
         snap_gmat(snap_file=args.snap_file,
                   gene_file=args.gene_file,
                   buffer_size=args.buffer_size,
                   tmp_folder=args.tmp_folder,
                   verbose=args.verbose)

    if args.command == "call-peak":
         from snaptools.call_peak import call_peak
         call_peak(snap_file=args.snap_file,
                   barcode_file=args.barcode_file,
                   gsize=args.gsize,
                   output_prefix=args.output_prefix,
                   path_to_macs=args.path_to_macs,
                   buffer_size=args.buffer_size,
                   macs_options=args.macs_options,
                   tmp_folder=args.tmp_folder)
    
    if args.command == "louvain":
         from snaptools.louvain import louvain
         louvain(edge_file=args.edge_file,
                 output_file=args.output_file,
                 resolution=args.resolution)
                        
def add_fastq_dex_subparser(subparsers):
     # create the parser for the "DMRfind" command
     parser_build = subparsers.add_parser(
          "decomplex-fastq",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Decomplex fastq file.")
     
     # add options
     parser_build_req = parser_build.add_argument_group("required inputs")
     parser_build_req.add_argument("--input-fastq",
                                   type=str,
                                   required=True,
                                   help="fastq file contains the sequencing reads")

     parser_build_req.add_argument("--output-fastq",
                                   type=str,
                                   required=True,
                                   help="output decomplexed fastq file")

     parser_build_req.add_argument("--index1-fastq",
                                   type=str,
                                   required=True,
                                   help="fastq file contains r7,i7 barcodes")

     parser_build_req.add_argument("--index2-fastq",
                                   type=str,
                                   required=True,
                                   help="fastq file contains r5,i5 barcodes")

     parser_build_req.add_argument("--index-list",
                                   type=str,
                                   required=True,
                                   help="pre-defined barcode list contains r7,i7,r5,i5 barcodes")

     
def add_index_ref_subparser(subparsers):
     # create the parser for the "DMRfind" command
     parser_build = subparsers.add_parser(
          "index-genome",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Index reference genome."
          )
     
     # add options
     parser_build_req = parser_build.add_argument_group("required inputs")
     parser_build_req.add_argument("--input-fasta",
                                   type=str,
                                   required=True,
                                   help="genome fasta file to build the index from")

     parser_build_opt = parser_build.add_argument_group("optional inputs")
     parser_build_opt.add_argument("--output-prefix",
                                   type=str,
                                   default=None,
                                   help="prefix of indexed file")

     parser_build_opt.add_argument("--aligner",
                                   type=str,
                                   default="bwa",
                                   help="aligner to use. Currently, snaptools supports bwa, bowtie, bowtie2 and minimap2.")

     parser_build_opt.add_argument("--path-to-aligner",
                                   type=str,
                                   default=None,
                                   help="path to fold that contains bwa")

     parser_build_opt.add_argument("--num-threads",
                                   type=int,
                                   default=3,
                                   help="=number of indexing threads")

def add_align_pe_subparser(subparsers):
     parser_build = subparsers.add_parser(
          "align-paired-end",
          formatter_class=argparse.ArgumentDefaultsHelpFormatter,
          help="Align paired-end reads.")
     
     # add options
     parser_build_req = parser_build.add_argument_group("required inputs")
     parser_build_req.add_argument("--input-reference",
                                   type=str,
                                   required=True,
                                   help="reference genome file contains the reference genome that reads are mapped against, "
                                   + "the genome index must be under the same folder")
    
     parser_build_req.add_argument("--input-fastq1",
                                   type=str,
                                   required=True,
                                   help="fastq file contains R1 reads, currently supports fastq, gz, bz2 file")
     
     parser_build_req.add_argument("--input-fastq2",
                                   type=str,
                                   required=True,
                                   help="fastq file contains R2 reads, currently supports fastq, gz, bz2 file")

     parser_build_req.add_argument("--output-bam",
                                   type=str,
                                   required=True,
                                   help="output bam file contains unfiltered alignments")
          
     parser_build_opt = parser_build.add_argument_group("optional inputs")
     parser_build_opt.add_argument("--aligner",
                                   type=str,
                                   default="bwa",
                                   help="aligner to use. Currently, snaptools supports bwa, bowtie, bowtie2 and minimap2.")
     
     parser_build_opt.add_argument("--path-to-aligner",
                                   type=str,
                                   default=None,
                                   help="path to fold that contains bwa")
     
     parser_build_opt.add_argument("--aligner-options",
                                   type=str,
                                   nargs="+",
                                   help="list of strings indicating options you would like passed to aligner"
                                   +"strongly do not recommand to change unless you know what you are doing. "
                                   +"the default is to align reads without filteration.")
     
     parser_build_opt.add_argument("--read-fastq-command",
                                   type=str,
                                   default=None,
                                   help="command line to execute for each of the input file. This command "
                                   +"should generate FASTQ text and send it to stdout. For example, "
                                   +"--read-fastq-command should be zcat, bzcat and cat for .gz, .bz2 and "
                                   +"plain fastq file respectively")
     
     parser_build_opt.add_argument("--min-cov",
                                   type=int,
                                   default=0,
                                   help="min number of fragments per barcode. barcodes of total fragments less than "
                                   +"--min-cov will be filreted before alingment. Note: though this feature is included, "
                                   +"we found it barely benefit anything, recommand to set it 0.")
     
     parser_build_opt.add_argument("--num-threads",
                                   type=int,
                                   default=1,
                                   help="number of alignment threads, also number of threads for sorting "
                                   +"a bam file.")
                                        
     parser_build_opt.add_argument("--if-sort",
                                   type=str2bool,
                                   default=True,
                                   help="weather to sort the bam file based on the read name")

     parser_build_opt.add_argument("--tmp-folder",
                                   type=str,
                                   default=None,
                                   help="directory to store temporary files. If not given, snaptools will automatically"
                                   + "generate a temporary location to store temporary files")
     
     parser_build_opt.add_argument("--overwrite",
                                   type=str2bool,
                                   default=False,
                                   help="whether to overwrite the output file if it already exists")
    
     parser_build_opt.add_argument("--verbose",
                                   type=str2bool,
                                   default=True,
                                   help="a boolen tag indicates output the progress.")    

def add_align_se_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "align-single-end",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Align single-end reads.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--input-reference",
                                  type=str,
                                  required=True,
                                  help="reference genome file contains the reference genome that reads are mapped against, "
                                  + "the genome index must be under the same folder")
   
    parser_build_req.add_argument("--input-fastq1",
                                  type=str,
                                  required=True,
                                  help="fastq file contains R1 reads, currently supports fastq, gz, bz2 file")
    
    parser_build_req.add_argument("--output-bam",
                                  type=str,
                                  required=True,
                                  help="output bam file contains unfiltered alignments")
         
    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--aligner",
                                  type=str,
                                  default="bwa",
                                  help="aligner to use. Currently, snaptools supports bwa, bowtie, bowtie2 and minimap2.")
    
    parser_build_opt.add_argument("--path-to-aligner",
                                  type=str,
                                  default=None,
                                  help="path to fold that contains bwa")
    
    parser_build_opt.add_argument("--aligner-options",
                                  type=str,
                                  nargs="+",
                                  help="list of strings indicating options you would like passed to aligner"
                                  + "strongly do not recommand to change unless you know what you are doing.")
    
    parser_build_opt.add_argument("--read-fastq-command",
                                  type=str,
                                  default=None,
                                  help="command line to execute for each of the input file. This command"
                                  +"should generate FASTA or FASTQ text and send it to stdout")
    
    parser_build_opt.add_argument("--num-threads",
                                  type=int,
                                  default=1,
                                  help="number of alignment threads")

    parser_build_opt.add_argument("--min-cov",
                                  type=int,
                                  default=0,
                                  help="min number of fragments per barcode. barcodes of total fragments less than "
                                  +"--min-cov will be filreted before alingment. Note: though this feature is included, "
                                  +"we found it barely speed up the process, so recommand to set it to be 0.")

    parser_build_opt.add_argument("--if-sort",
                                  type=str2bool,
                                  default=True,
                                  help="weather to sort the bam file based on the read name")

    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="directory to store temporary files. If not given, snaptools will automatically"
                                  + "generate a temporary location to store temporary files")

    parser_build_opt.add_argument("--overwrite",
                                  type=str2bool,
                                  default=False,
                                  help="whether to overwrite the output file if it already exists")


def add_snap_pre_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "snap-pre",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Create a snap file from bam or bed file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--input-file",
                                  type=str,
                                  required=True,
                                  help="input bam or bed file.")
   
    parser_build_req.add_argument("--output-snap",
                                  type=str,
                                  required=True,
                                  help="output snap file.")
    
    parser_build_req.add_argument("--genome-name",
                                  type=str,
                                  required=True,
                                  help="genome identifier (i.e. hg19, mm10). This tag does not change anything "
                                  +"unless merge or compare multiple snap files.")

    parser_build_req.add_argument("--genome-size",
                                  type=str,
                                  required=True,
                                  help="a txt file contains corresponding genome sizes. It must be in the following "
                                  +"format with the first column the chromsome name and the second column as chromsome "
                                  +"length. This tag does not change anything unless merge or compare multiple snap files.")
    
    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--barcode-file",
                                  type=str,
                                  default=None,
                                  help="a txt file contains pre-selected cell barcodes. "
                                  +"If --barcode-file is given, snaptools will ignore any barcodes not present in the --barcode-file. "
                                  +"If it is None, snaptools will automatically identify barcodes from bam file. "
                                  +"The first column of --barcode-file must be the selected barcodes and the other columns "
                                  +"could be any attributes of the barcode as desired (`ATGCTCTACTAC attr1 att2`). The "
                                  +"attributes, however, will not be kept in the snap file. This tag will be ignored if the output "
                                  +"snap file already exists.")

    parser_build_opt.add_argument("--min-mapq",
                                  type=int,
                                  default=10,
                                  help="min mappability score. Fargments with mappability score less than --min-mapq will be filtered.")
    
    parser_build_opt.add_argument("--min-flen",
                                  type=int,
                                  default=0,
                                  help="min fragment length. Fragments of length shorted than --min-flen will be filtered.")

    parser_build_opt.add_argument("--max-flen",
                                  type=int,
                                  default=1000,
                                  help="max fragment length. Fragments of length longer than --max-flen will be filtered.")
                                      
    parser_build_opt.add_argument("--min-cov",
                                  type=int,
                                  default=0,
                                  help="min number of fragments per barcode. barcodes of total fragments fewer than --min-cov will be considered "
                                  +"when creating the cell x bin count matrix. Note: because the vast majority of barcodes contains very few reads, "
                                  +"we found by setting --min-cov, one can remove barcodes of low coverage without wasting time and storage. "
                                  +"Please note that this is not selection of good barcodes for downstream clustering analysis, it is only filteration"
                                  +"of very low-quality barcodes.")
    
    parser_build_opt.add_argument("--keep-chrm",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates whether to keep fragments mapped to chrM. If set Fasle, fragments "
                                  +"aligned to the mitochondrial sequence will be filtered.")

    parser_build_opt.add_argument("--keep-single",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates whether to keep those reads whose mates are not mapped or missing. "
                                  +"If False, unpaired reads will be filtered. If True, unpaired reads will be simply treated "
                                  +"as a fragment. Note: for single-end such as scTHS-seq, --keep-single must be True.")

    parser_build_opt.add_argument("--keep-secondary",
                                  type=str2bool,
                                  default=False,
                                  help="a boolen tag indicates whether to keep secondary alignments. If False, "
                                  +"secondary alignments will be filtered. If True, a secondary alignments will "
                                  +"be treated as fragments just single-end.")

    parser_build_opt.add_argument("--keep-discordant",
                                  type=str2bool,
                                  default=False,
                                  help="a boolen tag indicates whether to keep discordant read pairs.")

    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")

    parser_build_opt.add_argument("--overwrite",
                                  type=str2bool,
                                  default=False,
                                  help="a boolen tag indicates whether to overwrite the matrix session if it already exists.")

    parser_build_opt.add_argument("--qc-file",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates whether to create a master qc file. This .qc file contains "
                                  +"basic quality control metrics at the bulk level. Quality control is only estimated by "
                                  +"selected barcodes only.")
    
    parser_build_opt.add_argument("--verbose",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates output the progress.")

def add_dump_read_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "dump-fragment",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Dump fragments of selected barcodes from a snap file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="input snap file.")

    parser_build_req.add_argument("--output-file",
                                  type=str,
                                  required=True,
                                  help="output bed file (supports .bed, .gz and .bz2 format).")

    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--barcode-file",
                                  type=str,
                                  required=False,
                                  help="a txt file of selected barcodes. If --barcode-file is given, only fragments of "
                                  + "selected barcodes are writen in the output file. Default None, all fragments "
                                  + "will be dumped into the output bed file.")

    parser_build_opt.add_argument("--buffer-size",
                                  type=int,
                                  default=1000,
                                  required=False,
                                  help="max number of barcodes be stored in the memory.")
    
    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store potential temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    
    parser_build_opt.add_argument("--overwrite",
                                   type=str2bool,
                                   default=False,
                                   help="overwrite the output file if it already exists.")

def add_dump_barcode_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "dump-barcode",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Dump barcodes from a snap file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="input snap file.")

    parser_build_req.add_argument("--output-file",
                                  type=str,
                                  required=True,
                                  help="output txt file.")

    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--barcode-file",
                                  type=str,
                                  required=False,
                                  help="a txt file of selected barcodes. If --barcode-file is given, only metadata of "
                                  + "selected barcodes are writen in the output file. Default None, meta data of all "
                                  + "barcodes will be dumped into the output file.")

    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store potential temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    
    parser_build_opt.add_argument("--overwrite",
                                   type=str2bool,
                                   default=False,
                                   help="overwrite the output file if it already exists.")

def add_snap_bmat_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "snap-add-bmat",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Add cell x bin count matrix to snap file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="snap file.")

    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--bin-size-list",
                                  type=int,
                                  nargs="+",
                                  default=[5000],
                                  help="a list of bin size(s) to create the cell-by-bin count matrix. The genome "
                                  +"will be divided into bins of the equal size of --bin-size-list to create the "
                                  +"cell x bin count matrix. If more than one bin size are given, snaptools will "
                                  +"generate a list of cell x bin matrices of different resolutions and stored in "
                                  +"the same snap file.")
    
    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    
    parser_build_opt.add_argument("--verbose",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates output the progress.")

def add_snap_pmat_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "snap-add-pmat",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Add cell x peak count matrix to snap file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="snap file.")

    parser_build_req.add_argument("--peak-file",
                                  type=str,
                                  required=True,
                                  help="bed file contains peaks.")

    parser_build_opt = parser_build.add_argument_group("optional inputs")    
    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    
    parser_build_opt.add_argument("--verbose",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates output the progress.")

def add_snap_gmat_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "snap-add-gmat",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Add cell x gene count matrix to snap file.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="snap file.")

    parser_build_req.add_argument("--gene-file",
                                  type=str,
                                  required=True,
                                  help="bed file contains genes.")

    parser_build_opt = parser_build.add_argument_group("optional inputs")    
    parser_build_opt.add_argument("--buffer-size",
                                  type=int,
                                  default=1000,
                                  required=False,
                                  help="max number of barcodes be stored in the memory.")

    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    
    parser_build_opt.add_argument("--verbose",
                                  type=str2bool,
                                  default=True,
                                  help="a boolen tag indicates output the progress.")


def add_call_peak_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "call-peak",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Call peak using selected barcodes.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--snap-file",
                                  type=str,
                                  required=True,
                                  help="snap file.")

    parser_build_req.add_argument("--barcode-file",
                                  type=str,
                                  required=True,
                                  help="a file contains selected barcodes.")

    parser_build_req.add_argument("--output-prefix",
                                  type=str,
                                  required=True,
                                  help="prefix of output files.")

    parser_build_req.add_argument("--path-to-macs",
                                  type=str,
                                  required=True,
                                  help="path to fold that contains macs2.")

    parser_build_req.add_argument("--gsize",
                                  type=str,
                                  required=True,
                                  help="effective genome size. 'hs' for human, 'mm' for mouse, 'ce' for C. elegans, 'dm' for fruitfly")

    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--buffer-size",
                                  type=int,
                                  default=1000,
                                  required=False,
                                  help="max number of barcodes be stored in the memory.")
    
    parser_build_opt.add_argument("--macs-options",
                                  type=str,
                                  nargs="+",
                                  help="list of strings indicating options you would like passed to macs2"
                                  +"strongly do not recommand to change unless you know what you are doing. "
                                  +"the default is '--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR -call-summits'.")

    parser_build_opt.add_argument("--tmp-folder",
                                  type=str,
                                  default=None,
                                  help="a directory to store temporary files. If not given, snaptools will automatically "
                                  + "generate a temporary location to store temporary files.")
    

def add_louvain_subparser(subparsers):
    parser_build = subparsers.add_parser(
         "louvain",
         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
         help="Louvain communities finding.")
    
    # add options
    parser_build_req = parser_build.add_argument_group("required inputs")
    parser_build_req.add_argument("--edge-file",
                                  type=str,
                                  required=True,
                                  help="txt file contains edges and weights.")

    parser_build_req.add_argument("--output-file",
                                  type=str,
                                  required=True,
                                  help="output file name.")

    parser_build_opt = parser_build.add_argument_group("optional inputs")
    parser_build_opt.add_argument("--resolution",
                                  type=float,
                                  default=1.0,
                                  required=False,
                                  help="resolution for finding communities.")
    
