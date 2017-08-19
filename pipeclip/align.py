#!/usr/bin/python
# programmer : bbc
# usage: align reads

import sys
import re
import pandas as pd
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from multiprocessing import Pool
import pysam
import pybedtools
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Align reads. For PE, sort by name"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="input fastq file. For pair-end, use comma to separate R1 and R2")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-g","--genome",dest = "genome",type=str, required=True, help="Genome index location")
  argparser.add_argument("--pe",dest="pe", default=False, action="store_true", help = "Set if input are pair-end")
  #argparser.add_argument("--rmdup",dest="rmdup", default=False, action="store_true", help = "Set if want to remove duplicates")
  #argparser.add_argument("--barcode",dest="barcode",type=str,help = "Barcode/randomer file")
  argparser.add_argument("--prog",dest="prog",type=str,default="bowtie",help = "Alignment program to use.['bowtie','tophat2']", choices=['bowtie','tophat2'])
  return(argparser)


def bowtie_align(f,out,genome, pe):
  if pe:
    f1,f2 = f.rstrip().split(",")
    bt_command = "bowtie -p 12 -a -m 1 "+genome+" -1 "+f1+" -2 "+f2+" -S "+out+".bowtie.sam 2>"+out+".bowtieLog"
  else: #single end
    bt_command = "bowtie -p 12 -q -m 1 "+genome+" -q "+f.rstrip()+" -S "+out+".bowtie.sam 2>"+out+".bowtieLog"
  logging.debug(bt_command)
  p1 = subprocess.Popen(bt_command, shell=True)
  p1.communicate()
  if pe:
    st_command = "samtools view -uhS "+out+".bowtie.sam | samtools sort -n - -T "+out+" > "+out+".bowtie.bam" 
  else:
    st_command = "samtools view -uhS "+out+".bowtie.sam | samtools sort - -T "+out+" > "+out+".bowtie.bam" 
  p2 = subprocess.Popen(st_command, shell=True)
  p2.communicate()
  return out+".bowtie.bam"

def bowtie_align_wrapper(args):
  bowtie_align(*args)

def tp2_align(f,out,genome, pe):
  if pe:
    f = f.rstrip().replace(","," ")
    tp2_command = "tophat2 --num-threads 12 --max-multihits 1 --output-dir "+out+"_tp2 "+genome+" "+f
  else: #single end
    tp2_command = "tophat2 --num-threads 12 --max-multihits 1 --output-dir "+out+"_tp2 "+genome+" "+f.rstrip()
  p1 = subprocess.Popen(tp2_command, shell=True)
  p1.communicate()
  if pe:
    st_command = "samtools sort -n "+out+"_tp2/accepted_hits.bam > "+out+".tp2.bam" 
  else:
    st_command = "samtools sort "+out+"_tp2/accepted_hits.bam > "+out+".tp2.bam" 
  p2 = subprocess.Popen(st_command, shell=True)
  p2.communicate()
  return out+".tp2.bam"

def tp2_align_wrapper(args):
  tp2_align(*args)

def run_align(f, out, genome, prog, pe, mode):
  if mode=="single":
    if prog == "bowtie":
      outfile = bowtie_align(f,out,genome,pe)
    elif prog == "tophat2":
      outfile = tp2_align(f,out,genome,pe)
  elif mode == "multiple":#parse design file here
    design = pd.read_table(f.rstrip())
    align_arglist_treat = zip(design['treat_align_input'].tolist(),design['sampleName'].tolist(),[genome]*design.shape[0],[pe]*design.shape[0])
    align_arglist_ctrl = zip(design['ctrl_align_input'].tolist(),design['sampleName'].tolist(),[genome]*design.shape[0],[pe]*design.shape[0])
    align_arglist = align_arglist_treat +  align_arglist_ctrl
    work_pool = Pool(min(12,design.shape[0]))
    if prog == "bowtie":
      resultlist = work_pool.map(bowtie_align_wrapper, align_arglist)
      design['treat_aligned_bam'] = design['sampleName'].str.cat(["bowtie.bam"]*design.shape[0],sep=".")
      design['ctrl_aligned_bam'] = design['sampleName'].str.cat(["bowtie.bam"]*design.shape[0],sep=".")
    else:
      resultlist = work_pool.map(tp2_align_wrapper, align_arglist)
      design['treat_aligned_bam'] = design['sampleName'].str.cat(["tp2.bam"]*design.shape[0],sep=".")
      design['ctrl_aligned_bam'] = design['sampleName'].str.cat(["tp2.bam"]*design.shape[0],sep=".")
    return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    aligned = run_align(args.infile, args.outprefix, args.genome, args.prog, args.pe,"single")
  elif args.dfile:
    newdesign = run_trim(args.dfile, None, args.genome, args.prog, args.pe, "multiple")
    newdesign.to_csv(args.dfile+".alignDesign",sep="\t",index=False)

if __name__=="__main__":
  main()
