#!/usr/bin/python
# programmer : bbc
# usage:

import sys
from scipy.stats import poisson
import argparse as ap
import logging

logging.basicConfig(level=10)


def prepare_argparser():
	description = "Calculate Poisson p value"
	epilog = "For command line options of each command, type %(prog)% COMMAND -h"
	argparser = ap.ArgumentParser(description=description, epilog = epilog)
	argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help="files with total coverage")
	argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
	argparser.add_argument("-t","--total",dest="total",type=float,required=True, help = "Total effective length")
	argparser.add_argument("-p","--pval",dest="pval",type=float,default=0.01)
	#argparser.add_argument("-n","--name",dest="trackName",type=str,default="UserTrack",help = "track name for bedgraph header")
	return(argparser)



def main():
	argparser = prepare_argparser()
	args = argparser.parse_args()

	try:
		infile = open(args.infile,"r")
		outfile = open(args.outfile,"w")
	except IOError,message:
		print >> sys.stderr, "cannot open file",message
		sys.exit(1)
	reads = [1,2,3,4,5,6,7,8,9,10]
	header = "file\ttotalcount\tgenomeLen\tbg0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\tmin_read_fit_p_"+str(args.pval)
	print >> outfile, header
	for row in infile:
		buf = row.rstrip().split("\t")
		bg0 = float(buf[1])/args.total
		prob = []
		min_read = 100
		for n in reads:
			p = 1-poisson.cdf(n,bg0)
			prob.append(str(p))
			if p<=args.pval and n<min_read:
				min_read = n
		print >>outfile,"\t".join(buf+[str(args.total),str(bg0)]+prob+[str(min_read)])

if __name__=="__main__":
	main()
