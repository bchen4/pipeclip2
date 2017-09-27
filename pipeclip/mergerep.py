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

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Filter peaks"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-f","--fraction",dest = "fraction",type=float, default=0.001, help="Fraction of the overlap")
  argparser.add_argument("--overlap-level",dest="level", default="two", type=str, help = "Peak overlap level", choices=['two','all'])
  return(argparser)

def add_bed_name(fn, name):
  logging.debug(fn+" "+name)
  infile = open(fn,"r")
  count = 1
  s = ""
  for row in infile:
    buf = row.split("\t")
    buf.insert(3, name+"_"+str(count)) 
    buf.insert(4, "0")
    count += 1
    s += "\t".join(buf)
  return BedTool(s,from_string=True)


def merge_intersect(feature):
  feature.start = min(int(feature[1]),int(feature[7]))#smaller start
  feature.stop = max(int(feature[2]),int(feature[8])) #smaller stop
  return feature[0:6]

def merge_rep(infile,level,out,fraction):
  files = infile['filter'] 
  logging.debug(files)
  bedfiles = []
  for f in files:
    bedfiles.append(BedTool(f))
  if level=="all":
    intersect_bed = bedfiles[0]
    for b in bedfiles[1:]:
      intersect_bed = intersect_bed.intersect(b,s=True,wo=True,f=fraction, F=fraction,e=True).each(merge_intersect)
    
    intersect_bed.saveas(out+"_mergerep.bed") 
  elif level == "two":
    intersect_list = []
    count = 1
    for a,b in combinations(bedfiles,2):
      a.intersect(b,s=True,wo=True,f=fraction, F=fraction,e=True).each(merge_intersect).saveas(out+".merge.tmp"+str(count))
      intersect_list.append(out+".merge.tmp"+str(count))
      count += 1
    intersect_bed = BedTool(intersect_list[0])  
    for f in intersect_list[0:]:#make sure it will merge
      intersect_bed = intersect_bed.cat(f,s=True).saveas(out+".merge.tmp")
    mergerep_bed =  add_bed_name(out+".merge.tmp", out)
    mergerep_bed.saveas(out+"_mergerep.bed")
  return out+"_mergerep.bed"

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.dfile:
    design = pd.read_table(args.dfile)
    newdesign = pd.DataFrame()
    for name, group in design.groupby('group'):
      merged_fn = merge_rep(group ,args.level, name, args.fraction)
      group['mergepeak'] = [merged_fn]*group.shape[0]
      newdesign = newdesign.append(group)
    #remove previour treat count, total, control count, total columns
    

if __name__=="__main__":
  main()
