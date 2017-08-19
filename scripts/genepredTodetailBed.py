#!/usr/bin/python
# programmer : bbc
# usage:Convert transcript pred to bed file, including UTR, exon, intron. Only for protein coding genes

import sys
import re
import random
import string

class transcriptpred:
  def __init__(self,info):
    self.transcript = info[0]
    self.chrom = info[1]
    self.strand = info[2]
    self.tss = int(info[3])
    self.tts = int(info[4])
    self.cds_start = int(info[5])
    self.cds_stop = int(info[6])
    self.exon_count = int(info[7])
    self.exon_starts = info[8].rstrip(",").split(",")
    self.exon_stops = info[9].rstrip(",").split(",")
    self.intron_starts = []
    self.intron_stops = []
    self.cds_exon_starts = []
    self.cds_exon_stops = []
    self.utr3_starts = []
    self.utr3_stops = []
    self.utr5_starts = []
    self.utr5_stops = []
    self.cds_exon_count = 0
    self.transcriptType = ""
    for i in range(len(self.exon_starts)):
      self.exon_starts[i] = int(self.exon_starts[i])
      self.exon_stops[i] = int(self.exon_stops[i])
    self.gene = info[10]
    if self.cds_start == self.cds_stop:
      self.transcriptType = "Noncoding"
    else:
      self.transcriptType = "Coding"
  
  def getIntron(self):
    for (i,j) in zip(self.exon_stops[0:self.exon_count-1],self.exon_starts[1:]):
      self.intron_starts.append(i)
      self.intron_stops.append(j)


  def getCDSExonInfo(self):
    #Find TSS, TTS position on exons
    #start_insert_index = self.exon_count - 1
    for i in reversed(range(self.exon_count)):
      if self.cds_start >=  self.exon_starts[i]:
        start_insert_index = i # start_insert_index should >=1
        break

    for j in range(self.exon_count):
      if self.cds_stop <= self.exon_stops[j]:
        stop_insert_index = j
        break
    #compare tss to exon start in front of it
    #print start_insert_index,stop_insert_index
    #print self.cds_start, self.exon_starts[0]
    for scan_index in range(start_insert_index, stop_insert_index+1):
      #append cds starts
      if self.cds_start <= self.exon_starts[scan_index]:
        self.cds_exon_starts.append(self.exon_starts[scan_index])
      else:
        self.cds_exon_starts.append(self.cds_start)
      #append cds stops
      if self.cds_stop >= self.exon_stops[scan_index]:
        self.cds_exon_stops.append(self.exon_stops[scan_index])
      else:
        self.cds_exon_stops.append(self.cds_stop)
      #check if lengths are the same
      if len(self.cds_exon_starts) == len(self.cds_exon_stops):
        self.cds_exon_count = len(self.cds_exon_starts)
      else:
        print >> sys.stderr,"CDS starts and stops don't match"
        self.cds_exon_count  = -1


  def getUTR5Info(self):
    #Find TSS, TTS position on exons
    #start_insert_index = self.exon_count - 1
    for j in range(self.exon_count):
      if self.cds_start > self.exon_starts[j]:
        self.utr5_starts.append(self.exon_starts[j])
        if self.cds_start <= self.exon_stops[j]:
          self.utr5_stops.append(self.cds_start)
          break
        else:
          self.utr5_stops.append(self.exon_stops[j])

  def getUTR3Info(self):
    #Find TSS, TTS position on exons
    #start_insert_index = self.exon_count - 1
    #print self.cds_stop
    for j in reversed(range(self.exon_count)):
      #print self.exon_stops[j]
      if self.cds_stop < self.exon_stops[j]:
        self.utr3_stops.insert(0,self.exon_stops[j])
        if self.cds_stop >= self.exon_starts[j]:
          self.utr3_starts.insert(0,self.cds_stop)
          break
        else:
          self.utr3_starts.insert(0,self.exon_starts[j])
    #print self.utr3_starts
    #print self.utr3_stops

  def printCDSpred(self):
    st = "\t".join([self.transcript,self.chrom,self.strand,str(self.tss),str(self.tts),str(self.cds_start),str(self.cds_stop), str(self.cds_exon_count)] )
    starts_str = []
    stops_str = []
    for i in range(self.cds_exon_count):
      starts_str.append(str(self.cds_exon_starts[i]))
      stops_str.append(str(self.cds_exon_stops[i]))
    st += "\t" + ",".join(starts_str)+","
    st += "\t" + ",".join(stops_str)+","
    st += "\t"+self.gene
    return st

  def printCDSExonBed(self):
    for i in range(self.cds_exon_count):
      if self.strand == "+":
        yield "\t".join([self.chrom,str(self.cds_exon_starts[i]),str(self.cds_exon_stops[i]),self.transcript,str(i+1),self.strand,self.gene])
      elif self.strand == "-":
        yield "\t".join([self.chrom,str(self.cds_exon_starts[i]),str(self.cds_exon_stops[i]),self.transcript,str(self.cds_exon_count - i),self.strand,self.gene])


  def printUTRBed(self):
    if self.strand=="+":
      UTR3_starts = self.utr3_starts
      UTR3_stops = self.utr3_stops
      UTR5_starts = self.utr5_starts
      UTR5_stops = self.utr5_stops
    else:
      UTR3_starts = self.utr5_starts
      UTR3_stops = self.utr5_stops
      UTR5_starts = self.utr3_starts
      UTR5_stops = self.utr3_stops

    for i in range(len(UTR3_starts)):
      if (UTR3_stops[i]-UTR3_starts[i])>0:
        yield "\t".join([self.chrom,str(UTR3_starts[i]),str(UTR3_stops[i]),self.transcript,str(i+1),self.strand,self.gene,"UTR3"])
 
    for i in range(len(UTR5_starts)):
      if (UTR5_stops[i]-UTR5_starts[i])>0:
        yield "\t".join([self.chrom,str(UTR5_starts[i]),str(UTR5_stops[i]),self.transcript,str(i+1),self.strand,self.gene,"UTR5"])

  def printBed12(self, tag):
    if tag == 'UTR3':
      if self.strand=="+":
        UTR3_starts = self.utr3_starts
        UTR3_stops = self.utr3_stops
      else:
        UTR3_starts = self.utr5_starts
        UTR3_stops = self.utr5_stops
      starts = UTR3_starts
      stops = UTR3_stop
    elif tag == "UTR5":
      if self.strand=="+":
        UTR5_starts = self.utr5_starts
        UTR5_stops = self.utr5_stops
      else:
        UTR5_starts = self.utr3_starts
        UTR5_stops = self.utr3_stops
      starts = UTR5_starts
      stops = UTR5_stops
    elif tag =="CDS":
      starts = self.cds_starts
      stops = self.cds_stops
    elif tag == "intron":
      starts = self.intron_starts
      stops = self.intron_stops
    elif tag == "exon":
      starts = self.exon_stars
      stops = self.exon_stops
    else:
      logging.error("Tag must be: UTR3/UTR5/CDS/exon/intron")
      sys.exit(0)
    st = "\t".join([self.chrome,str(starts[0]),str(stops[-1]),
      self.transcript],"0",self.strand,str(starts[0]),str(stops[-1]),
      "255,0,0") 
    sizes = []
    length = 0
    offsets = []
    for i in range(len(starts)):
      if (UTR3_stops[i]-UTR3_starts[i])>0:
        s = UTR3_stops[i]-UTR3_starts[i]
        length += s
        sizes.append(str(s))
        offsets.append(str(starts[i]-starts[0]))
    st += "\t"+str(len(sizes))
    st += "\t"+",".join(sizes)
    st += "\t"+",".join(offsets)
    return st
  

  def printListBed(self,tag):
    if tag == "cdsexon":
      starts = self.cds_exon_starts
      stops = self.cds_exon_stops
    elif tag =="intron":
      starts = self.intron_starts
      stops = self.intron_stops
    elif tag == "exon":
      starts = self.exon_starts
      stops = self.exon_stops
    else:
      #print >> sys.stderr("Tag error, please use \'exon\' or \'intron\'")
      sys.exit(0)
    for i in range(len(starts)):
      if self.strand == "+":
        yield "\t".join([self.chrom,str(starts[i]),str(stops[i]),self.transcript,str(i+1),self.strand,self.gene, self.transcriptType+"_"+tag])
      elif self.strand == "-":
        yield "\t".join([self.chrom,str(starts[i]),str(stops[i]),self.transcript,str(self.cds_exon_count - i),self.strand,self.gene, self.transcriptType+"_"+tag])
  
  def printGenepred(self, tag):
    if tag == 'UTR3':
      if self.strand=="+":
        UTR3_starts = self.utr3_starts
        UTR3_stops = self.utr3_stops
        try:
          UTR3_starts[0]+=3 #remove the stop codon
        except:#there is no UTR3
          return ""
      else:
        UTR3_starts = self.utr5_starts
        UTR3_stops = self.utr5_stops
        try:
          UTR3_stops[-1]-=3 #remove stop codon
        except:
          return ""
      starts = UTR3_starts
      stops = UTR3_stops
    elif tag == "UTR5":
      if self.strand=="+":
        UTR5_starts = self.utr5_starts
        UTR5_stops = self.utr5_stops
      else:
        UTR5_starts = self.utr3_starts
        UTR5_stops = self.utr3_stops
      starts = UTR5_starts
      stops = UTR5_stops
    elif tag =="CDS":
      starts = self.cds_starts
      stops = self.cds_stops
    elif tag == "intron":
      starts = self.intron_starts
      stops = self.intron_stops
    elif tag == "exon":
      starts = self.exon_stars
      stops = self.exon_stops
    else:
      logging.error("Tag must be: UTR3/UTR5/CDS/exon/intron")
      sys.exit(0)
    st = "\t".join([self.transcript, self.chrom,self.strand, str(starts[0]),str(stops[-1]),
      str(starts[0]),str(starts[0])])
    starts_str = []
    stops_str = []
    for i in range(len(starts)):
      if stops[i]>starts[i]:
        starts_str.append(str(starts[i]))
        stops_str.append(str(stops[i]))
    if len(starts_str)==0:
      return ""
    else:
      st += "\t"+str(len(starts_str))
      st += "\t"+",".join(starts_str)+","
      st += "\t"+",".join(stops_str)+","
      st += "\t"+self.gene
      return st
  
  def sumLen(self,starts,stops):
    total = 0
    for i in range(len(starts)):
      total += (stops[i]-starts[i])
    return total

  def getPartsLen(self):
    '''Calculate total lenght of CDS, UTR3, UTR5 and intron'''
    transcript_span = self.tts - self.tss
    cds_span = self.cds_stop - self.cds_start
    total_utr = transcript_span - cds_span
    total_cds = self.sumLen(self.cds_exon_starts, self.cds_exon_stops)
    total_utr3 = self.sumLen(self.utr3_starts, self.utr3_stops)
    total_utr5 = self.sumLen(self.utr5_starts, self.utr5_stops)
    total_intron = self.sumLen(self.intron_starts, self.intron_stops)
    calculate_intron = cds_span - total_cds
    return [str(transcript_span),str(cds_span),str(total_utr),
        str(total_cds),str(total_utr3), str(total_utr5), str(total_intron),
        str(calculate_intron)]
  

def main():
  try:
    infile = open(sys.argv[1],"r+")
  except IOError,message:
    print >> sys.stderr, "cannot open file",message
    sys.exit(1)
  
  #print "\t".join(["transcript_span","cds_span","total_utr","total_cds",
   #       "total_utr3", "total_utr5", "total_intron","calculate_intron"])
  for item in infile:
    buf = item.rstrip().split("\t")
    new_isoform  =  transcriptpred(buf)
    if new_isoform.transcriptType == "Coding":
      new_isoform.getCDSExonInfo()
      new_isoform.getIntron()
      new_isoform.getUTR3Info()
      new_isoform.getUTR5Info()
      #utr3_genepred_str = new_isoform.printGenepred("UTR3")
      #if len(utr3_genepred_str)>0:
      #  print utr3_genepred_str
      #    "total_utr3", "total_utr5", "total_intron","calculate_intron"])
    #  print "\t".join([new_isoform.transcript]+new_isoform.getPartsLen())
      #print new_isoform.printCDSpred()
      for exon in new_isoform.printListBed("cdsexon"):
        print exon
      for intron in new_isoform.printListBed("intron"):
        print intron
      for utr in new_isoform.printUTRBed():
        print utr
    else: #noncoding or does not have a type somehow
      new_isoform.getIntron()
      for exon in new_isoform.printListBed("exon"):
        print exon
      for intron in new_isoform.printListBed("intron"):
        print intron

if __name__=="__main__":
  main()
