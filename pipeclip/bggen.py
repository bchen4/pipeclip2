#!/usr/bin/python
# programmer : bbc
# usage: Generate bedgraph

import sys
import re
import pandas as pd
import subprocess
from multiprocessing import Pool
import pysam
import pybedtools
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Generate bedgraph"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="input fastq file. For pair-end, use comma to separate R1 and R2")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-g","--genome",dest = "genome",type=str, required=True, help="Genome chromsome length file")
  argparser.add_argument("--scale",dest="scale", default=False, action="store_true", help = "Set if need to generate scaled bedgraph")
  argparser.add_argument("--pe",dest="pe", default=False, action="store_true", help = "Set if pair end")
  return(argparser)

def scale_bg(input, outfile, factor=0):
  infile = pd.read_table(input,header=None)
  infile.columns = ['chr','start','stop','score']
  if factor>0:
    infile['scale_score'] = infile['score']*1000000/factor
  else:
    infile['scale_score'] = infile['score']*1000000/sum(infile.score)
  infile = infile.loc[:,['chr','start','stop','scale_score']]
  infile.to_csv(outfile,sep="\t",header=False, index=False)



def bg_gen(f,out,genome,scale, pe):
  if pe:
    logging.info("Not supported right now. Function in development")
    sys.exit()
  else:
    #sort bam file first
    sort_command = "samtools sort -o "+f.replace("bam","sort.bam")+" "+f
    p1 = subprocess.Popen(sort_command, shell=True)
    p1.communicate()
    index_command = "samtools index "+ f.replace("bam","sort.bam")
    p2 = subprocess.Popen(index_command, shell=True)
    p2.communicate()
    sorted_bam  = f.replace("bam","sort")+".bam"
    bg1_command = " ".join(["bedtools genomecov -bg -split -strand + -ibam",sorted_bam,"-g",genome,">",out+"_positive.bg"])
    bg2_command = " ".join(["bedtools genomecov -bg -split -strand - -ibam",sorted_bam,"-g",genome,">",out+"_negative.bg"])
    p3 = subprocess.Popen(bg1_command, shell=True)
    p3.communicate()
    p4 = subprocess.Popen(bg2_command, shell=True)
    p4.communicate()
    if scale:
      #logging.debug("Start to scale bedgraph")
      #calculate scale factor: total mapped reads
      bamfile = pysam.AlignmentFile(sorted_bam,"rb")
      scale_factor = bamfile.mapped
      scale_bg(out+"_positive.bg", out+"_positive.scaled.bg", scale_factor)
      scale_bg(out+"_negative.bg", out+"_negative.scaled.bg", scale_factor)

def bg_gen_wrapper(args):
  bg_gen(*args)

def run_bg_gen(f, out, genome, scale, pe, mode):
  if mode=="single":
    bg_gen(f,out,genome,scale, pe)
  elif mode == "multiple":#parse design file here
    design = pd.read_table(f.rstrip())
    bggen_arglist_treat = zip(design['treat_rmdup_bam'].tolist(),design['sampleName'].tolist(),[genome]*design.shape[0],[scale]*design.shape[0],[pe]*design.shape[0])
    bggen_arglist_ctrl = zip(design['ctrl_rmdup_bam'].tolist(),design['sampleName'].tolist(),[genome]*design.shape[0],[scale]*design.shape[0],[pe]*design.shape[0])
    work_pool = Pool(min(12,design.shape[0]))
    work_pool.map(bg_gen_wrapper, bggen_arglist)
    design['treat_pos_bg'] = design['sampleName'].str.cat(["positive_treat.bg"]*design.shape[0],sep="_")
    design['treat_neg_bg'] = design['sampleName'].str.cat(["negative_treat.bg"]*design.shape[0],sep="_")
    design['ctrl_pos_bg'] = design['sampleName'].str.cat(["positive_ctrl.bg"]*design.shape[0],sep="_")
    design['ctrl_neg_bg'] = design['sampleName'].str.cat(["negative_ctrl.bg"]*design.shape[0],sep="_")
    return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_bg_gen(args.infile,args.outprefix,args.genome, args.scale, args.pe,"single")
  elif args.dfile:
    newdesign = run_bg_gen(args.dfile,None,args.genome, args.scale, args.pe,"multiple")
    newdesign.to_csv(args.dfile+".bgbin",sep="\t",index=False)

if __name__=="__main__":
  main()
