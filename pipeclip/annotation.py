#!/usr/bin/python
# programmer : bbc
# usage: Annotate peaks

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
  group.add_argument("-i","--input",dest = "infile",type=str, help="peak file prefix. there should be two files for each prefix, .xls and .bed")
  group.add_argument("-d","--design",dest = "dfile",type=str, help="design file")
  argparser.add_argument("-o","--output",dest = "outprefix",type=str, help="output prefix for single input. For design file, use sample name to generate output file names")
  argparser.add_argument("-g","--genome",dest = "genome",type=str, required=True, help="Detailed bed annotation")
  #argparser.add_argument("-c","--ctrl",dest="ctrlbam",type=str, required=True, help = "Set if need to generate scaled bedgraph")
  #argparser.add_argument("--peak-col",dest="peak_col", default="peakbed", type=str, help = "Column name of the peaks", choices=['peakbed','mergepeak'])
  return(argparser)





def prioritizeAnnotation(gdf):
  sorter = ['Coding_cdsexon','UTR3','UTR5','Coding_intron','Noncoding_exon','Noncoding_intron','Intergenic']
  gdf['annotation_type'] = gdf['annotation_type'].astype('category')
  gdf['annotation_type'].cat.set_categories(sorter, inplace=True)
  gdf = gdf.sort_values('annotation_type')
  return gdf.iloc[0,:]

def anno_filter(feature):
  #logging.debug(feature)
  f_len = len(feature)
  new_feature = feature[0:6]+feature[f_len-3:f_len]
  #logging.debug(new_feature)
  return new_feature

def annotate(bfile,xfile, out, genome):
  xls = pd.read_table(xfile)
  bed = pybedtools.BedTool(bfile)
  anno = pybedtools.BedTool(genome)
  anno_bed = bed.intersect(anno, s=True, f=0.55, wao=True).saveas(out+"_anno_intersect")
  logging.debug(anno_bed.count())
  annobed_df = pd.read_table(out+"_anno_intersect",header=None)
  logging.debug("read data frame finished")
  annobed_df = annobed_df.iloc[:,[0,1,2,3,4,5,12,13,14]]
  header = ['chr','start','stop','id','summit_cov','strand','gene_id','annotation_type','overlap_len']
  annobed_df.columns = header
  annobed_df = annobed_df.replace(".","Intergenic")
  df = pd.DataFrame(columns=header)
  for name,group in annobed_df.groupby('id'):
    if group.shape[0]>1:
      newdf = prioritizeAnnotation(group)
      df = df.append(newdf)
    else:
      df = df.append(group)
  result_df = xls.merge(df.loc[:,["id","gene_id","annotation_type","overlap_len"]],on="id",how="left").fillna("None")
  result_df.to_csv(out+".annotation.xls",sep="\t",index=False)
  typecount = result_df['annotation_type'].value_counts()
  countdf = pd.DataFrame({'annotation':typecount.index,'counts':typecount.values})
  countdf.to_csv(out+".annoStat.xls",sep="\t",index=False)


def annotate_wrapper(args):
  annotate(*args)


def  run_annotate(f, out, genome, mode):
  if mode=="single":
   annotate(f+".bed",f+".xls",out,genome)
  elif mode == "multiple":#parse design file here
    design = pd.read_table(f.rstrip())
    logging.debug(design.shape)
    annotation_arglist = zip(design['final_peak'].str.cat([".bed"]*design.shape[0],sep=""),design['final_peak'].str.cat([".xls"]*design.shape[0],sep=""),design['final_peak'].tolist(),[genome]*design.shape[0])
    logging.debug(annotation_arglist) 
    work_pool = Pool(min(12,design.shape[0]))
    work_pool.map(annotate_wrapper, annotation_arglist)
    

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  if args.infile:
    run_annotate(args.infile,args.outprefix,args.genome, "single")
  elif args.dfile:
    run_annotate(args.dfile,None,args.genome,"multiple")

if __name__=="__main__":
  main()
