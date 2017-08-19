#!/usr/bin/python
# programmer : bbc
# usage: Filter peaks  

import sys
import re
import pandas as pd
import numpy as np
import subprocess
from multiprocessing import Pool
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Filter peaks"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("--rc",dest = "rc",type=int, default=5, help="Reads count in IP cutoff")
  argparser.add_argument("--fc",dest = "fc",type=float, default=-1.0, help="fold change cutoff")
  argparser.add_argument("--fp",dest="fp",type=float, default=-1.0, help = "Fisher p value cutoff")
  argparser.add_argument("--yp",dest="yp",type=float, default=-1.0, help = "Yates p value cutoff")
  #argparser.add_argument("--peak-col",dest="peak_col", default="peakbed", type=str, help = "Column name of the peaks", choices=['peakbed','mergepeak'])
  return(argparser)


def filter(infile, rc,fc, fp, yp,out):
  df = pd.read_table(infile)
  logging.debug(df.columns) 
  df = df[df['treat_count']>= rc]
  if fc > 0:
    df = df[abs(df['norm_log2fc'])>=np.log2(fc)]
  elif fc == 0:
    df = df[abs(df['norm_log2fc'])>0]
  if fp>0:
    df =df[df['fisherP']<=fp]
  if yp>0:
    df = df[df['yatesP']<=yp]
  df = df.iloc[:,range(6)]
  df.to_csv(out+"_normfc_filter.bed",sep="\t",index=False,header=False)
  return out+"_normfc_filter.bed"

def filter_wrapper(args):
  return filter(*args)

def run_filter_rep(design,rc,fc,fp,yp):
  filter_arglist = zip(design['normfc_txt'].tolist(),[rc]*design.shape[0],[fc]*design.shape[0],[fp]*design.shape[0],[yp]*design.shape[0],design["sampleName"])
  work_pool = Pool(min(12,design.shape[0]))
  resultlist = work_pool.map(filter_wrapper, filter_arglist)
  logging.debug(resultlist)
  df = pd.DataFrame({'filter':resultlist},index=design.index)
  design = pd.concat([design,df],axis=1)
  return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.dfile:
    design = pd.read_table(args.dfile)
    newdesign = pd.DataFrame()
    for name, group in design.groupby('group'):
      groupdesign = run_filter_rep(group ,args.rc, args.fc, args.fp, args.yp)
      if newdesign.shape[0]==0:
        newdesign = groupdesign
      else:
        newdesign = newdesign.append(groupdesign)
    newdesign.to_csv(args.dfile+".filter",sep="\t",index=False)

if __name__=="__main__":
  main()
