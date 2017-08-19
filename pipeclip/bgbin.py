#!/usr/bin/python
# programmer : bbc
# usage: Generate bedgraph bin

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
  description = "Generate bedgraph bin, output bed format"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="input bedgraph file")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-s","--strand",dest = "strand",type=str, default=".", help="Strand for the bedgraph")
  argparser.add_argument("--peak-name",dest="pname", default="peak",  help = "Used for bed file id generation")
  #argparser.add_argument("--pe",dest="pe", default=False, action="store_true", help = "Set if pair end")
  argparser.add_argument("--bin-min",dest="binmin", type=int, default=3, help = "Minimal coverage to trigger bin start/stop")
  argparser.add_argument("--bin-peak",dest="binpeak",type=int, default=10, help = "Minimal peak coverage")
  return(argparser)



def bg_bin(f,out,binmin, binpeak, pname, strand):
  #initial = infile.readline().rstrip().split("\t")
  try:
    infile = open(f,"r")
    outfile = open(out,"a") #append mode for multiple mode
  except IOError:
    logging.error("Cannot open bedgraph file, exit.")
    sys.exit(1)
  else:
    #current_chr = ""#initial[0]
    #current_start = ""#initial[1]
    #current_stop = ""#initial[2]
    #peak = 0#int(initial[3])
    nameprefix = pname
    count = 0
    ini = infile.readline().rstrip().split("\t")
    #logging.debug(ini)
    current_chr = ini[0]
    current_start = ini[1]
    current_stop = ini[2]
    previous_stop = ini[1]
    flag = True
    peak = int(ini[3])
    for row in infile:
      #logging.debug(row)
      buf = row.rstrip().split("\t")
      if not flag:#there is no current bin
        #logging.debug("create a bin")
        if int(buf[3])>=binmin:#start a bin
          #logging.debug("start a bin")
          current_chr=buf[0]
          current_start = buf[1]
          previous_stop = buf[2]#iniciatie: need to compare to itself
          current_stop = buf[2]
          flag = True
          peak = int(buf[3])
      else:#For existing bin
        #check if the chr is still the same
        #logging.debug("compare to exsiting bin")
        if buf[0]==current_chr and buf[1]==previous_stop:
         # logging.debug("Same chr")
          if int(buf[3])>peak:
            peak = int(buf[3])#update peak value
          previous_stop = buf[2]
          if int(buf[3])<binmin:#stop of a peak, elif to if 8/11/2017
            if peak>=binpeak:#if peak value passed the threshold
            #print bin
              count += 1
              b_name = pname+"_"+str(count)
              #logging.debug("output peak")
              print >>outfile,"%s\t%s\t%s\t%s\t%d\t%s" % (current_chr,current_start,current_stop,b_name,peak,strand)
            #initialize current_chr and flag
            current_chr = ""
            flag = False
          else:#value is between 3 and peak, update bed stop
            current_stop = buf[2]
        else:#chr is not the same, check and print previous bin and start a new bin
          #logging.debug("not connected with previous record")
          if peak >=binpeak:
            count += 1
            b_name = nameprefix+"_"+str(count)
            print >>outfile,"%s\t%s\t%s\t%s\t%d\t%s" % (current_chr, current_start,current_stop, b_name,peak,strand)
          if int(buf[3])>=binmin:
            #start a new bin
            current_chr = buf[0]
            current_start = buf[1]
            current_stop = buf[2]
            previous_stop = buf[2]
            peak = int(buf[3])
            flag = True
          else:#set everything to empty
            current_chr = ""
            current_start = ""
            current_stop = ""
            peak = 0
            flag = False
    #check if there is still a bin
    if current_chr != "":
      if peak >= binpeak:
        count += 1
        b_name = pname+"_"+str(count)
        print >>outfile,"%s\t%s\t%s\t%s\t%d\t%s" % (current_chr,current_start,current_stop,b_name,peak,strand)
    outfile.close()

def bg_bin_wrapper(args):
  bg_bin(*args)

def run_bg_bin(f, out, bmin, bpeak, pname, strand, mode):
  if mode=="single":
    bg_bin(f,out,bmin,bpeak, pname,strand)
  elif mode == "multiple":#parse design file here
    logging.debug("Run bgBin multiple")
    design = pd.read_table(f.rstrip())
    bgbin_arglist_pos = zip(design['treat_pos_bg'].tolist(),design['sampleName'].str.cat(['_bgbin.bed']*design.shape[0],sep=""),[bmin]*design.shape[0],[bpeak]*design.shape[0],design['sampleName'].str.cat(["_pos"]*design.shape[0],sep=""),['+']*design.shape[0])
    bgbin_arglist_neg = zip(design['treat_neg_bg'].tolist(),design['sampleName'].str.cat(['_bgbin.bed']*design.shape[0],sep=""),[bmin]*design.shape[0],[bpeak]*design.shape[0],design['sampleName'].str.cat(["_neg"]*design.shape[0],sep=""),['-']*design.shape[0])
    #bgbin_arglist = bgbin_arglist_pos + bgbin_arglist_neg
    work_pool = Pool(min(12,design.shape[0]))
    resultlist = work_pool.map(bg_bin_wrapper, bgbin_arglist_pos)
    work_pool = Pool(min(12,design.shape[0]))
    resultlist = work_pool.map(bg_bin_wrapper, bgbin_arglist_neg)

    design['peakbed'] = design['sampleName'].str.cat(["_bgbin.bed"]*design.shape[0],sep="")
    return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_bg_bin(args.infile,args.outprefix,args.binmin,args.binpeak, args.pname, args.strand,"single")
  elif args.dfile:
    newdesign = run_bg_bin(args.dfile,None,args.binmin,args.binpeak, None, None, "multiple")
    newdesign.to_csv(args.dfile+".bgbin",sep="\t",index=False)

if __name__=="__main__":
  main()
