#!/usr/bin/python
# programmer : bbc
# usage: trim adapters as well as randomers

import sys
import re
import pandas as pd
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from multiprocessing import Pool
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Trim adapter and/or randomers for eclip"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="input fastq file. For pair-end, use comma to separate R1 and R2")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("--pe",dest="pe", default=False, action="store_true", help = "Set if input are pair-end")
  argparser.add_argument("--ops",dest="ops",type=str,default="all", choices=["cutadapter","rmbarcode","all"])
  argparser.add_argument("--barcode-len",dest="barcodelen",type=int,default=10,help = "Barcode/randomer length")
  return(argparser)

def run_cutadapter(f, out, pe):
  logging.info("Running cutadapt")
  if pe:
    f1,f2 = f.rstrip().split(",")
    out1 = out+"_R1_cutadapt.fastq"
    out2 = out+"_R2_cutadapt.fastq"
    metric = out+"_cutadapt.metric"
    command = " ".join(["sh cutadapt_skeleton_pe.sh",f1, f2, out1, out2, metric])
  else:
    out1 = out+"_R2_cutadapt.fastq"
    out2 = ""
    metric = out+"_cutadapt.metric"
    command = " ".join(["sh cutadapt_skeleton_se.sh",f.rstrip(), out1, metric])
  logging.info(command)
  p = subprocess.Popen(command, shell=True)
  p.communicate()
  logging.info("Cutadapter finished")
  return ",".join([out1,out2]).rstrip(",")

def run_cutadapter_wrapper(args):
  run_cutadapter(*args)

def run_rmbarcode(f,out,blen,pe):
  if pe:
    f1,f2 = f.rstrip().split(",")
    input_file = f2
  else:
    input_file = f.rstrip()
  fq_out = open(out+"_rmbarcode.fastq","w")
  bc_out = open(out+"_barcode.txt","w")
  for title, seq, qual in FastqGeneralIterator(open(input_file)):
    print >> fq_out, "@"+str(title)
    print >> fq_out, seq[blen:]
    print >> fq_out, "+"
    print >> fq_out, qual[blen:]
    print >> bc_out, "\t".join([str(title),seq[0:blen]])
  fq_out.close()
  bc_out.close()

def run_rmbarcode_wrapper(args):
  run_rmbarcode(*args)

def run_trim(f, op, out, pe, blen, mode):
  infile = f
  if mode=="single":
    if op in ['cutadapter','all']:
      infile = run_cutadapter(infile, out, pe)
    if op in ['rmbarcode','all']:
      run_rmbarcode(infile, out,blen,pe)
  elif mode == "multiple":#parse design file here
    design = pd.read_table(f.rstrip())
    if op in ['cutadapter','all']:
      cutadapter_arglist_treat = zip(design['treat_fastq'].tolist(),design['sampleName'].tolist(),[pe]*design.shape[0])
      cutadapter_arglist_ctrl = zip(design['ctrl_fastq'].tolist(),design['sampleName'].tolist(),[pe]*design.shape[0])
      cutadapter_arglist = cutadapter_arglist_treat + cutadapter_arglist_ctrl
      work_pool = Pool(min(12,design.shape[0]))
      resultlist = work_pool.map(run_cutadapter_wrapper, cutadapter_arglist)
      #modify design file separately since resultlist won't maintain order
      design['treat_align_r1'] = design['sampleName'].str.cat(["R1_cutadapt_treat.fastq"]*design.shape[0],sep="_")
      design['treat_cutadapt_metric'] = design['sampleName'].str.cat(["cutadapt_treat.metric"]*design.shape[0],sep="_")
      design['ctrl_align_r1'] = design['sampleName'].str.cat(["R1_cutadapt_ctrl.fastq"]*design.shape[0],sep="_")
      design['ctrl_cutadapt_metric'] = design['sampleName'].str.cat(["cutadapt_ctrl.metric"]*design.shape[0],sep="_")
      if op=="all":
        design['treat_rmbc_input'] = design['sampleName'].str.cat(["R2_cutadapt_treat.fastq"]*design.shape[0],sep="_")
        design['ctrl_rmbc_input'] = design['sampleName'].str.cat(["R2_cutadapt_ctrl.fastq"]*design.shape[0],sep="_")
      elif op=="cutadapt":#for alignment directly
        design['treat_align_r1'] = design['sampleName'].str.cat(["R2_cutadapt_treat.fastq"]*design.shape[0],sep="_")
        design['ctrl_align_r1'] = design['sampleName'].str.cat(["R2_cutadapt_ctrl.fastq"]*design.shape[0],sep="_")

    if op in ['rmbarcode','all']:
      if op == "all":
        rmbarcode_arglist_treat = zip(design['treat_rmbc_input'].tolist(),design['sampleName'].tolist(), [blen]*design.shape[0], [pe]*design.shape[0])
        rmbarcode_arglist_ctrl = zip(design['ctrl_rmbc_input'].tolist(),design['sampleName'].tolist(), [blen]*design.shape[0], [pe]*design.shape[0])
      else:
        rmbarcode_arglist_treat = zip(design['treat_fastq'].tolist(),design['sampleName'].tolist(), [blen]*design.shape[0], [pe]*design.shape[0])
        rmbarcode_arglist_ctrl = zip(design['ctrl_fastq'].tolist(),design['sampleName'].tolist(), [blen]*design.shape[0], [pe]*design.shape[0])
      rmbarcode_arglist = rmbarcode_arglist_treat + rmbarcode_arglist_ctrl
      work_pool = Pool(min(12,design.shape[0]))
      work_pool.map(run_rmbarcode_wrapper, rmbarcode_arglist)
      design['treat_align_r2'] = design['sampleName'].str.cat(["rmbarcode_treat.fastq"]*design.shape[0],sep="_")
      design['treat_ctrl_r2'] = design['sampleName'].str.cat(["rmbarcode_ctrl.fastq"]*design.shape[0],sep="_")
      design['treat_barcode_txt'] = design['sampleName'].str.cat(["barcode_treat.txt"]*design.shape[0],sep="_")
      design['ctrl_barcode_txt'] = design['sampleName'].str.cat(["barcode_ctrl.txt"]*design.shape[0],sep="_")
      if pe:
        design['treat_align_input'] = design['treat_align_r1'].str.cat(design['treat_align_r2'],sep=",")
        design['ctrl_align_input'] = design['ctrl_align_r1'].str.cat(design['ctrl_align_r2'],sep=",")
      else:
        design['treat_align_input'] = design['treat_align_r2']
        design['ctrl_align_input'] = design['ctrl_align_r2']
      design.drop(['treat_align_r1','treat_align_r2','ctrl_align_r1','ctrl_align_r2'], axis=1, inplace=True)
      return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_trim(args.infile, args.ops, args.outprefix, args.pe, args.barcodelen,"single")
  elif args.dfile:
    newdesign = run_trim(args.dfile, args.ops, None, args.pe, args.barcodelen, "multiple")
    newdesign.to_csv(args.dfile+".trimDesign",sep="\t",index=False)

if __name__=="__main__":
  main()
