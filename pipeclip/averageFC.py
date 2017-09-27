#!/usr/bin/python
# programmer : bbc
# usage: Merge replicates and re-count 

import sys
import re
import pandas as pd
import numpy as np
import subprocess
from multiprocessing import Pool
from pybedtools import BedTool
from itertools import combinations
import argparse as ap
import logging
import normfc

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Filter peaks"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
 # argparser.add_argument("-f","--fraction",dest = "fraction",type=float, default=0.001, help="Fraction of the overlap")
  #argparser.add_argument("--fc",dest = "fc",type=float, default=-1.0, help="fold change cutoff")
  #argparser.add_argument("--fp",dest="fp",type=float, default=-1.0, help = "Fisher p value cutoff")
  #argparser.add_argument("--yp",dest="yp",type=float, default=-1.0, help = "Yates p value cutoff")
  #argparser.add_argument("--overlap-level",dest="level", default="two", type=str, help = "Peak overlap level", choices=['two','all'])
  return(argparser)

def add_bed_name(fn, name):
  infile = open(fn,"r")
  count = 1
  s = ""
  for row in fn:
    buf = row.split("\t")
    buf.insert(3, name+"_"+str(count)) 
    buf.insert(4, "0")
    count += 1
    s += "\t".join(buf)
  return BedTool(s, from_string=True)


def average_fc(df):
  peaks_list = []
  for name,group in df.groupby('group'):
    logging.debug(name)
    peak = []
    treat_cols = []
    ctrl_cols = []
    newdf = pd.DataFrame()
    for f, rep in zip(group['mergepeak_normfc_txt'],group['replicates']):
      fcfile = pd.read_table(f)
      fcfile['count_treat_rep'+str(rep)] = fcfile['treat_count']
      fcfile['count_ctrl_rep'+str(rep)] = fcfile['ctrl_count']
      fcfile['total_treat_rep'+str(rep)] = fcfile['treat_total']
      fcfile['total_ctrl_rep'+str(rep)] = fcfile['ctrl_total']
      fcfile['norm_count_treat_rep'+str(rep)] = (fcfile['treat_count']+1)*1000000/fcfile['treat_total']
      fcfile['norm_count_ctrl_rep'+str(rep)] = (fcfile['ctrl_count']+1)*1000000/fcfile['ctrl_total']
      treat_cols.append('norm_count_treat_rep'+str(rep))
      ctrl_cols.append('norm_count_ctrl_rep'+str(rep))
      if newdf.shape[0]==0:
        newdf = fcfile.loc[:,['chr','start','stop','id','summit_cov','strand','count_treat_rep'+str(rep),'count_ctrl_rep'+str(rep),'total_treat_rep'+str(rep),'total_ctrl_rep'+str(rep),'norm_count_treat_rep'+str(rep),'norm_count_ctrl_rep'+str(rep)]]
      else:
        newdf = newdf.merge(fcfile.loc[:,['id','count_treat_rep'+str(rep),'count_ctrl_rep'+str(rep),'total_treat_rep'+str(rep),'total_ctrl_rep'+str(rep),'norm_count_treat_rep'+str(rep),'norm_count_ctrl_rep'+str(rep)]])
    #calculate average
    newdf['ave_log2fc'] = newdf.apply(lambda x:np.log2(x[treat_cols].sum()/x[ctrl_cols].sum()),axis=1)
    newdf = newdf.drop_duplicates()
    newdf.to_csv(name+"_finalpeak.xls",sep="\t",index=False)
    newdf.iloc[:,range(6)].to_csv(name+"_finalpeak.bed",sep="\t",index=False,header=False)
    peaks_list.append(name+"_finalpeak")

  return peaks_list


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  if args.dfile:
    recount_design = pd.read_table(args.dfile)
    final_peaks = pd.DataFrame({'final_peak':average_fc(recount_design)})
    final_peaks.to_csv(args.outprefix,index=False)
    

if __name__=="__main__":
  main()
