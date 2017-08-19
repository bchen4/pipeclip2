#!/usr/bin/python
# programmer : bbc
# usage: Merge replicates and re-count 

import sys
import re
import pandas as pd
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
  #argparser.add_argument("--rc",dest = "rc",type=int, default=5, help="Reads count in IP cutoff")
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
  fields = feature.fields
  record[1] = min(int(fields[1]),int(fields[6]))#smaller start
  record[2] = max(int(fields[2]),int(fields[7])) #smaller stop
  return record

def merge_rep(infile,level,out):
  files = infile['filter'] 
  bedfiles = []
  for f in files:
    bedfiles.append(BedTool(f))
  logging.debug(bedfiles)
  if level=="all":
    intersect_bed = bedfiles[0]
    for b in bedfiles[1:]:
      intersect_bed = intersect_bed.intersect(b,s=True,wo=True).each(merge_intersect)
    
    intersect_bed.saveas(out+"_mergerep.tmp") 
  elif level == "two":
    intersect_list = []
    count = 1
    for a,b in combinations(bedfiles,2):
      a.intersect(b,s=True,wo=True).each(merge_intersect).saveas(out+".merge.tmp"+str(count))
      intersect_list.append(out+".merge.tmp"+str(count))
      count += 1
    intersect_bed = BedTool(intersect_list[0])  
    for f in intersect_list[1:]:
      intersect_bed = intersect_bed.cat(f,s=True)#.saveas(out+".merge.tmp")
    intersect_bed.saveas(out+"_mergerep.tmp")
  #mergerep_bed = add_bed_name(out+"_mergerep.tmp", out)
  #mergerep_bed.saveas(out+"_mergerep.bed")
  return out+"_mergerep.bed"


def average_fc(df):
  peaks_list = []
  for name,group in df.groupby('group'):
    peak = []
    treat_cols = []
    ctrl_cols = []
    newdf = pd.DataFrame()
    for f, rep in zip(group['normfc_txt'],group['replicates']):
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
    newdf['ave_log2fc'] = newdf.apply(lambda x: sum(x.loc[:,treat_cols])/sum(x.loc[:,ctrl_cols]),axis=1)
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
    newdesign.to_csv(args.dfile+".mergerep",sep="\t",index=False)
  #re-count and calculate average fold change
    #BC#logging.debug("Start to recount")
    #BC#recount_design = normfc.run_normfc(args.dfile+".mergerep",None,None,None,"mergepeak","multiple")
  # Get average fc
    #BC#final_peaks = pd.DataFrame({'final_peak':average_fc(recount_design)})
    #BC#final_peaks.to_csv(args.out,index=False)
    

if __name__=="__main__":
  main()
