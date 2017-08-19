#!/usr/bin/python
# programmer : bbc
# usage: Calculate normalize FC(treat/ctrl) on given peak bed

import sys
import re
import pandas as pd
import subprocess
from multiprocessing import Pool
import pysam
import scipy.stats
import pybedtools
import numpy as np
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Generate bedgraph"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group(required=True)
  group.add_argument("-i","--input",dest = "infile",type=str, help="peak bam file")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-t","--treat",dest = "treatbam",type=str, help="treatment bam file")
  argparser.add_argument("-c","--ctrl",dest="ctrlbam",type=str, help = "Set if need to generate scaled bedgraph")
  argparser.add_argument("--peak-col",dest="peak_col", default="peakbed", type=str, help = "Column name of the peaks", choices=['peakbed','mergepeak'])
  return(argparser)

def normalized_log2_fc(tc,tn,cc,cn):
  return np.log2(((tc+1)/tn)/((cc+1)/cn))

def yatesP(a,b,c,d):
  test = [[a,b-a],[c+1,d-c-1]]
  chi2, p, ddof, expected = scipy.stats.chi2_contingency( test )
  return

def fisherP(a,b,c,d):
  test = [[a,b-a],[c+1,d-c-1]]
  odds, p = scipy.stats.fisher_exact( test )
  return p


def readcount(peak, bam,out):
  #Use bedtools to get readscount
  #logging.debug(bam)
  peakfile = pybedtools.BedTool(peak)
  bamfile = pybedtools.BedTool(bam)
  intersect = peakfile.intersect(bamfile, wa=True, s=True, F=0.51, c=True).remove_invalid().saveas(out)
  #logging.debug((out,bamfile.count()))
  return (out,bamfile.count())

def readcount_wrapper(args):
  return readcount(*args)


def normfc(tfile, cfile, t_total, c_total, out):
  treatcount = pd.read_table(tfile, header=None)
  ctrlcount = pd.read_table(cfile, header=None)
  treatcount.columns = ['chr','start','stop','id','summit_cov','strand','treat_count']
  ctrlcount.columns = ['chr','start','stop','id','summit_cov','strand','ctrl_count']
  ctrlcount = ctrlcount.loc[:,['id','ctrl_count']]
  df = treatcount.merge(ctrlcount, how="outer", on="id")
  df['treat_total'] = [t_total*1.0]*df.shape[0]
  df['ctrl_total'] = [c_total*1.0]*df.shape[0]
  df['norm_log2fc'] = df.apply(lambda x: normalized_log2_fc(x['treat_count'],x['treat_total'],x['ctrl_count'],x['ctrl_total']),axis=1)
  df['yatesP'] = df.apply(lambda x : yatesP(x['treat_count'],x['treat_total'],x['ctrl_count'],x['ctrl_total']),axis=1)
  df['fisherP'] = df.apply(lambda x : fisherP(x['treat_count'],x['treat_total'],x['ctrl_count'],x['ctrl_total']),axis=1)
  df.to_csv(out+"_normfc.txt",sep="\t",index=False)
  return out+"_normfc.txt"

def normfc_wrapper(args):
  return normfc(*args)

def run_normfc(f, out, treat, ctrl, peakcol, mode):
  if mode=="single":
    treat_count, treat_total = readcount(f,treat,out+"_treat.readscount")
    ctrl_count, ctrl_total = readcount(f,ctrl,out+"_ctrl.readscount")
    normfc(out+"_treat.readscount", out+"_ctrl.readscount", treat_total, ctrl_total,out)
  elif mode == "multiple":#parse design file here
    design = pd.read_table(f.rstrip())
    count_arglist_treat = zip(design[peakcol].tolist(),design['treat_rmdup_bam'],design['sampleName'].str.cat(["_treat.readscount"]*design.shape[0],sep=""))
    count_arglist_ctrl = zip(design[peakcol].tolist(),design['ctrl_rmdup_bam'],design['sampleName'].str.cat(["_ctrl.readscount"]*design.shape[0],sep=""))
    #logging.debug(count_arglist_treat)
    #logging.debug(count_arglist_ctrl)  
    work_pool_treat = Pool(min(12,design.shape[0]))
    resultlist_treat = work_pool_treat.map(readcount_wrapper, count_arglist_treat)
    #logging.debug(resultlist_treat)
    df1 = pd.DataFrame(resultlist_treat,index=design.index,columns=["treat_count","treat_total"])
    #logging.debug("treatment reads count finished")  
    work_pool_ctrl = Pool(min(12,design.shape[0]))
    resultlist_ctrl = work_pool_ctrl.map(readcount_wrapper, count_arglist_ctrl)
    #logging.debug(resultlist_ctrl)
    df2 = pd.DataFrame(resultlist_ctrl,index=design.index,columns=["ctrl_count","ctrl_total"])
    #logging.debug("control reads count finished")
    design = pd.concat([design, df1, df2],axis=1)
    ##logging.debug(design)
    #start to calculate fold change
    normfc_arglist = zip(design['treat_readscount'], design['ctrl_readscount'],design['treat_total'],design['ctrl_total'],design['sampleName'])
    #logging.debug(normfc_arglist)
    work_pool = Pool(min(12,design.shape[0]))
    outfile_list = work_pool.map(normfc_wrapper, normfc_arglist)
    design = pd.concat([design, pd.DataFrame({"normfc_txt":outfile_list}, index=design.index)],axis=1)
    return design

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_normfc(args.infile,args.outprefix,args.treatbam, args.ctrlbam,args.peak_col,"single")
  elif args.dfile:
    newdesign = run_normfc(args.dfile,None,None,None,args.peak_col,"multiple")
    newdesign.to_csv(args.dfile+".normfc",sep="\t",index=False)

if __name__=="__main__":
  main()
