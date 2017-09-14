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
  argparser.add_argument("-f","--fraction",dest = "fraction",type=float, default=0.001, help="Fraction of the overlap")
  #argparser.add_argument("--fc",dest = "fc",type=float, default=-1.0, help="fold change cutoff")
  #argparser.add_argument("--fp",dest="fp",type=float, default=-1.0, help = "Fisher p value cutoff")
  #argparser.add_argument("--yp",dest="yp",type=float, default=-1.0, help = "Yates p value cutoff")
  argparser.add_argument("--overlap-level",dest="level", default="two", type=str, help = "Peak overlap level", choices=['two','all'])
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


def merge_intersect(feature):
  feature.start = min(int(feature[1]),int(feature[7]))#smaller start
  feature.stop = max(int(feature[2]),int(feature[8])) #smaller stop
  return feature[0:6]

def merge_rep(infile,level,out):
  files = infile['filter'] 
  bedfiles = []
  for f in files:
    bedfiles.append(BedTool(f))
  if level=="all":
    intersect_bed = bedfiles[0]
    for b in bedfiles[1:]:
      intersect_bed = intersect_bed.intersect(b,s=True,wo=True,f=args.fraction, F=args.fraction,e=True).each(merge_intersect)
    
    intersect_bed.saveas(out+"_mergerep.bed") 
  elif level == "two":
    intersect_list = []
    count = 1
    for a,b in combinations(bedfiles,2):
      a.intersect(b,s=True,wo=True,f=args.fraction, F=args.fraction,e=True).each(merge_intersect).saveas(out+".merge.tmp"+str(count))
      intersect_list.append(out+".merge.tmp"+str(count))
      count += 1
    intersect_bed = BedTool(intersect_list[0])  
    for f in intersect_list[1:]:
      intersect_bed = intersect_bed.cat(f,s=True)#.saveas(out+".merge.tmp")
    intersect_bed.saveas(out+"_mergerep.bed")
  #mergerep_bed = add_bed_name(out+"_mergerep.tmp", out)
  #mergerep_bed.saveas(out+"_mergerep.bed")
  return out+"_mergerep.bed"


def average_fc(df):
  peaks_list = []
  for name,group in df.groupby('group'):
    logging.debug(name)
    peak = []
    treat_cols = []
    ctrl_cols = []
    newdf = pd.DataFrame()
    for f, rep in zip(group['normfc_txt'],group['replicates']):
      logging.debug(f)
      logging.debug(rep)
      fcfile = pd.read_table(f)
      fcfile['norm_count_treat_rep'+str(rep)] = (fcfile['treat_count']+1)*1000000/fcfile['ctrl_total']
      fcfile['norm_count_ctrl_rep'+str(rep)] = (fcfile['ctrl_count']+1)*1000000/fcfile['ctrl_total']
      treat_cols.append('norm_count_treat_rep'+str(rep))
      ctrl_cols.append('norm_count_ctrl_rep'+str(rep))
      if newdf.shape[0]==0:
        newdf = fcfile.loc[:,['chr','start','stop','id','summit_cov','strand','norm_count_treat_rep'+str(rep),'norm_count_ctrl_rep'+str(rep)]]
      else:
        newdf = newdf.merge(fcfile.loc[:,['id','norm_count_treat_rep'+str(rep),'norm_count_ctrl_rep'+str(rep)]])
    #calculate average
    logging.debug(treat_cols)
    logging.debug(newdf.columns)
    newdf['ave_log2fc'] = newdf.apply(lambda x:np.log2(x[treat_cols].sum()/x[ctrl_cols].sum()),axis=1)
    newdf.to_csv(name+"_finalpeak.xls",sep="\t",index=False)
    newdf.iloc[:,range(6)].to_csv(name+"_finalpeak.bed",sep="\t",index=False)
    peaks_list.append(name+"_finalpeak")

  return peaks_list


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.dfile:
    design = pd.read_table(args.dfile)
    newdesign = pd.DataFrame()
    for name, group in design.groupby('group'):
      merged_fn = merge_rep(group ,args.level, name)
      group['mergepeak'] = [merged_fn]*group.shape[0]
      newdesign = newdesign.append(group)
    #remove previour treat count, total, control count, total columns
    try:
      newdesign.drop(['treat_count','treat_total','ctrl_count','ctrl_total','normfc_txt'],axis=1, inplace=True)
    except:
      logging.warning("Failed to drop the columns")
    newdesign.to_csv(args.dfile+".mergerep",sep="\t",index=False)
  #re-count and calculate average fold change
    logging.debug("Start to recount")
    recount_design = normfc.run_normfc(args.dfile+".mergerep",None,None,None,"mergepeak","multiple")
  # Get average fc
    recount_design.to_csv("mergerep_recount.design",sep="\t",index=False)
    #recount_design = pd.read_table("mergerep_recount.design")
    final_peaks = pd.DataFrame({'final_peak':average_fc(recount_design)})
    final_peaks.to_csv(args.outprefix,index=False)
    

if __name__=="__main__":
  main()
