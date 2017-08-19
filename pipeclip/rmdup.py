#!/usr/bin/python
# programmer : bbc
# usage: remove PCR duplication using barcode file

import sys
import re
import pandas as pd
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from multiprocessing import Pool
from collections import Counter
import pysam
import pybedtools
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Remove PCR duplicates using barcode file"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="input fastq file. For pair-end, use comma to separate R1 and R2")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("--pe",dest="pe", default=False, action="store_true", help = "Set if input are pair-end")
  argparser.add_argument("--barcode",dest="barcode",type=str,help = "Barcode/randomer file")
  return(argparser)

def run_rmdup_se_wrapper(args):
  rmdup_se(*args)

def rmdup_se(infile,barfile,outfile,metricfile):
  barcode_dict = {}
  removed_count = Counter()
  total_count = Counter()
  results = {}
  mfile = open(metricfile,"w")
  barcode_file = open(barfile,"r")
  for row in barcode_file:
    buf = re.split('\t| ',row.rstrip())
    barcode_dict[buf[0]] = buf[-1]
  #start to parse bam file
  samfile = pysam.Samfile(infile,"rb")
  outbam = pysam.Samfile(outfile,'wb',template=samfile)
  for read in samfile:
    if not read.is_unmapped:    
      try:
        barcode = barcode_dict[read.qname]
      except:
        logging.warning("No barcode for read "+read.qname+", skip.")
      else:
        start = read.pos
        stop = read.positions[-1]
        strand = "-" if read.is_reverse else "+"
        unique_location = (read.rname, start, stop, strand, barcode)
        total_count[barcode] += 1
        if unique_location in results:
          removed_count[barcode] += 1
          continue
        results[unique_location] = read
    
  for ur in results.values():
    outbam.write(ur)
  outbam.close()
  #write metrics file
  #mfile.write("#Concordant pairs: "+str(concordant_count)+"\n")
  mfile.write("#Total count: "+str(sum(total_count.values()))+"\n")
  mfile.write("#Removed count: "+str(sum(removed_count.values()))+"\n")
  for b in total_count.keys():
    mfile.write("\t".join([b,str(total_count[b]),str(removed_count[b])])+"\n")
  mfile.close()


def run_rmdup(f,out,barcode,pe,mode):
    if pe:
      logging.info("Function to be added. Please use single end")
    else:
      if mode == "single":
        rmdup_se(f,barcode,out+"_rmdup.bam",out+"_rmdup.metric")
      elif mode == "multiple":
        design = pd.read_table(f)
        rmdup_arglist_treat = zip(design['treat_aligned_bam'],design['treat_barcode_txt'],design['sampleName'].str.cat(['rmdup_treat.bam']*design.shape[0],sep="_"),design['sampleName'].str.cat(['rmdup_treat.metric']*design.shape[0],sep="_"))
        rmdup_arglist_treat = zip(design['ctrl_aligned_bam'],design['ctrl_barcode_txt'],design['sampleName'].str.cat(['rmdup_ctrl.bam']*design.shape[0],sep="_"),design['sampleName'].str.cat(['rmdup_ctrl.metric']*design.shape[0],sep="_"))
        work_pool = Pool(min(12,design.shape[0]))
        resultlist = work_pool.map(run_rmdup_se_wrapper, rmdup_arglist)
        design['treat_rmdup_bam'] = design['sampleName'].str.cat(['rmdup_treat.bam']*design.shape[0],sep="_")
        design['ctrl_rmdup_bam'] = design['sampleName'].str.cat(['rmdup_ctrl.bam']*design.shape[0],sep="_")
        return design


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_rmdup(args.infile, args.barcode, args.outprefix, args.pe, "single")
  elif args.dfile:
    new_design = run_rmdup(args.dfile, None, None, args.pe, "multiple")
    design.to_csv(args.dfile+".rmdup", sep="\t", index=False)

if __name__=="__main__":
  main()
