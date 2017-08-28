# PIPE-CLIP2

This is the updated version for PIPE-CLIP, designed for CLIP-seq experiment data analysis.

# Major changes

1. Included eclip support;
2. Included control file support;
3. Included multi-processing
4. Included mapping

# Version Limits

## Alpha

- Can only deal with single-end sequencing
- Best for eclip data
- Assume the reads are the same strand as genes
- Have to manually use each scripts and provide a design file each time

# To-dos

- Make the connection (push-button) scripts to make it easier to use and check points for each step
- Change the design file foramtting so the output of previous script can be used directly for the next script
- Make setup scripts
- Add pair-end support
- Add old PIPE-CLIP functionality

# Report bugs

Please DO NOT email me about bugs. Please either:

1. Submit an issue on github (so it will be tracked and documented)
2. Talk to me directly

# Workflow explain

The whole workflow includes following modules:

## Fastq processing

![](https://static.notion-static.com/0bfd2227371c4c648ffc62831a0dcc2e/Screen_Shot_2017-08-25_at_2.46.31_PM.png)

- Remove sequencing adapter for the original fastq file. Adapter sequences can be changed in the cutadapt_skelenton.sh scripts;
- For eCLIP and iCLIP data, after remove adapters, barcode of the sequences need to be removed and recorded. Make sure you know the length of the barcode;
- Align the reads to the genome using bowtie or tophat. Other aligners could be used.
- For HITS-CLIP and PAR-CLIP data, the duplicates could be removed using
  - The start and the strand of the reads;
  - The start, stop, strand and the sequences of the reads;
  - The will be deprecated and replaced by picard;
- For eCLIP and iCLIP data, the barcode file and the aligned BAM file used to get rid of PCR duplicates

## Peak calling

![](https://static.notion-static.com/bd4594b11eaf45d68da14645cc2282b2/Screen_Shot_2017-08-25_at_4.36.32_PM.png)

- IP BAM is used for peak calling in two different ways
  - Bedgraph moving window (recommended): a trigger coverage and summit coverage will be used to call peaks. The peak will start if the coverage reaches the trigger coverage and stop when the coverage falls below the trigger. Then the program will check if the summit of the region reaches the given summit cutoff. If reaches, the peak will be called
  - Zero-Truncated-Negative-Binomial: program will merge all the reads into bins and find out the significant bins fitting a ZTNB distribution and take the bin length into account at the same time
- IP BAM is used for mutation calling (optional):
  - PAR-CLIP: 4SU and 6xx will induce T to C or X to X mutations in the sequences and this information could be used to varify the peaks
  - The program will find all mutations in the genome and use them to calculate the background mutation rate. Then for each matching mutation, the coverage on that location(k) and number of reads has that mutation (m) will be counted as used to calculate the significance against background using negative binomial distribution
- Peak determination
  - For each peak candidates, the reads count of IP and control samples will be obtained and normalized (to total usable reads) fold change, Fisher's P value and Yate's P value will be calculated for further filtering.
    - From previous experience, the fold change cut off makes the cleanest peaks
  - If there is no control sample, either the peaks called or the peaks containing reliable mutations will be reported

## Merge replicates

![](https://static.notion-static.com/bd1d9a86d0624b1ca911424dac8ec8f7/Screen_Shot_2017-08-25_at_4.36.40_PM.png)

- When there is multiple replicates for one condition, the peaks will be merged if
  - They appear in at least two replicates
  - They appear in all replicates (default)
- The count of each replicates' IP and control will be re-counted since the length of the peaks may change
- Average normalized fold change will be calculatedt

# PIPE-CLIP2 workflow (eCLIP as an example)

1.  **Trim the adapter and/or barcode.** 
```python
    usage: trim.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] [--pe]
     [--ops {cutadapter,rmbarcode,all}] [--barcode-len BARCODELEN]
    
    Trim adapter and/or randomers for eclip
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     input fastq file. For pair-end, use comma to separate
     R1 and R2
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use
     sample name to generate output file names
     --pe Set if input are pair-end
     --ops {cutadapter,rmbarcode,all}
     --barcode-len BARCODELEN
     Barcode/randomer length
```

Use -i to provide fastq file if it is single sample mode. Use -d to provide design file when there are multiple files.

Example:

    # Trim adapter and barcode
    python [trim.py](http://trim.py) -d design.txt --ops all --barcode-len 10
    
    # Trim adapter only
    python [trim.py](http://trim.py) -d design.txt --ops cutadapt
    
    # Trim barcode only
    python [trim.py](http://trim.py) -d design.txt --ops barcode --barcode-len 10

Design file for trimming

- TXT format (tab delimited)
- Required columns
  - fastqs: the path to your fastq files
  - sampleName: the names of your sample. Output file names will use this

 **2.Align the reads** 

    usage: align.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] -g GENOME [--pe]
     [--prog {bowtie,tophat2}]
    
    Align reads. For PE, sort by name
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     input fastq file. For pair-end, use comma to separate
     R1 and R2
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use
     sample name to generate output file names
     -g GENOME, --genome GENOME
     Genome index location
     --pe Set if input are pair-end
     --prog {bowtie,tophat2}
     Alignment program to use.['bowtie','tophat2']

Use -i to provide fastq file if it is single sample mode. Use -d to provide design file when there are multiple files.

-g is pre-made genome index. Please make sure you use the one fits the aligner you chose

Example:

    python [align.py](http://align.py) -d design.txt -g /path-to-genome-index/mm10 --prog bowtie

Design file for alignment

- TXT format (tab delimited)
- Required columns
  - fastqs: the path to your fastq files. If the input is pair-end, use comma to separte R1 and R2
  - sampleName: the names of your sample. Output file names will use this

 **3.Remove PCR duplicates** 

    usage: rmdup.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] [--pe]
     [--barcode BARCODE]
    
    Remove PCR duplicates using barcode file
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     input fastq file. For pair-end, use comma to separate
     R1 and R2
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use sample name to generate output file names
     --pe Set if input are pair-end
     --barcode BARCODE Barcode/randomer file

Design file for PCR duplicates removal

- TXT format (tab delimited)
- Required columns
  - aligned_bam: the path to aligned BAM files
  - sampleName: the names of your sample. Output file names will use this
  - barcode_txt: the path to the barcode txt file generated by trim.py

 **4.Make bedgraph file** 

    usage: bggen.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] -g GENOME [--scale]
     [--pe]
    
    Generate bedgraph
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     input fastq file. For pair-end, use comma to separate
     R1 and R2
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use
     sample name to generate output file names
     -g GENOME, --genome GENOME
     Genome chromsome length file
     --scale Set if need to generate scaled bedgraph
     --pe Set if pair end

This script will separate the positive and negative strand and generates two bedgraph files respectively.

Design file for bedgraph file generation

- TXT format (tab delimited)
- Required columns
  - rmdup_bam: the path to your usable BAM file
  - sampleName: the names of your sample. Output file names will use this

 **5.Generate bins** 

    usage: bgbin.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] [-s STRAND]
     [--peak-name PNAME] [--bin-min BINMIN] [--bin-peak BINPEAK]
    
    Generate bedgraph bin, output bed format
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     input bedgraph file
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use sample name to generate output file names
     -s STRAND, --strand STRAND
     Strand for the bedgraph
     --peak-name PNAME Used for bed file id generation
     --bin-min BINMIN Minimal coverage to trigger bin start/stop
     --bin-peak BINPEAK Minimal peak coverage

strand: since the bedgraph file does not contain strand information and CLIP-seq's strand is important, please indicate the bedgraph of the strand so the bed file could include the information

binmin could be calculated using poissonCutoff.py in the scripts folder. The minimal is set to be 3. When the coverage of the experiments is poor, 1 read coverage will be significant. The program will use a minimum of 3 reads as trigger.

Default of binpeak is 10.

Design file for bedgraph bin generation

- TXT format (tab delimited)
- Required columns
  - sampleName: the names of your sample. Output file names will use this
  - pos_bg: the path to your bedgraph for positive strand
  - neg_bg: the path to your bedgraph for negative strand

 **6.Count reads and calculate normalized fold change** 

    usage: normfc.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] [-t TREATBAM]
     [-c CTRLBAM] [--peak-col {peakbed,mergepeak}]
    
    Generate bedgraph
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     peak bam file
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use sample name to generate output file names
     -t TREATBAM, --treat TREATBAM
     treatment bam file
     -c CTRLBAM, --ctrl CTRLBAM
     Set if need to generate scaled bedgraph
     --peak-col {peakbed,mergepeak}
     Column name of the peaks

This scripts will be used for initial normalization and for the merged peaks for conditions. If use design file, need to tell the program which column to look at.

Design file for bedgraph file generation (format will change dramatically, future updates will use this format)

- TXT format (tab delimited)
- Required columns for initial count
  - sampleName: the names of your sample. Output file names will use this
  - treat_rmdup_bam: the path to your IP usable BAM file
  - ctrl_rmdup_bam: the path to your control usable BAM file
  - peakbed: peak file for this sample

- Required columns for merge peak count
  - sampleName: the names of your sample. Output file names will use this
  - treat_rmdup_bam: the path to your IP usable BAM file
  - ctrl_rmdup_bam: the path to your control usable BAM file
  - mergepeak: merged peak file for this condition
  - group: group information
  - replicates: replicates number

 **7.Filter peaks** 

    usage: filterfc.py [-h] [-d DFILE] [-o OUTPREFIX] [--rc RC] [--fc FC]
     [--fp FP] [--yp YP]
    
    Filter peaks
    
    optional arguments:
     -h, --help show this help message and exit
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use
     sample name to generate output file names
     --rc RC Reads count in IP cutoff
     --fc FC fold change cutoff
     --fp FP Fisher p value cutoff
     --yp YP Yates p value cutoff

Read cutoff has a default minimum of 5.

For other cutoffs, it not set, it won't be applied.

Example:

    python [filterfc.py](http://filterfc.py) -d design.txt --rc 3 --fc 2

Design file for bedgraph bin generation

- TXT format (tab delimited)
- Required columns
  - sampleName: the names of your sample. Output file names will use this
  - normfc_txt: the path to your file output containing normalized fold change

 **8.Merge replicates** 

    usage: mergerep.py [-h] [-d DFILE] [-o OUTPREFIX] [--overlap-level {two,all}]
    
    Filter peaks
    
    optional arguments:
     -h, --help show this help message and exit
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use
     sample name to generate output file names
     --overlap-level {two,all}
     Peak overlap level

This scripts only accepts design file. Each line of the design file represents a IP-ctrl pair.

Overlap level:

- They appear in at least two replicates
- They appear in all replicates (default)

Final output will be merged peaks with averge fold change.

Design file for bedgraph bin generation

- TXT format (tab delimited)
- Required columns
  - sampleName: the names of your sample. Output file names will use this
  - filter: the path to your file after filtered by fold change
  - group: group information
  - replicates: replicates number

 **10.Annotation** 

    usage: annotation.py [-h] (-i INFILE | -d DFILE) [-o OUTPREFIX] -g GENOME
    
    Annotate peak file
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     peak file prefix. there should be two files for each
     prefix, .xls and .bed
     -d DFILE, --design DFILE
     design file
     -o OUTPREFIX, --output OUTPREFIX
     output prefix for single input. For design file, use sample name to generate output file names
     -g GENOME, --genome GENOME
     Detailed bed annotation

All peaks will be annotated and if a peak is annotated to different genomic regions. Half of the peak have to be in that genomic region to be annotated. If a peak is annotated to different regions due to alternative splicing and overlapping genes, it will be sorted in the following order:

Coding_cdsexon > UTR3 > UTR5 > Coding_intron > Noncoding_exon > Noncoding_intron> Intergenic

Genome annotation is a extended BED format, which is the basic 6-col BED plus gene name and genomic feature name (coding_exon, UTR, intron, etc.) This could be generated using the script genepredTodetailBed.py

# Additional scripts

poissonCutoff.py

    usage: poissonCutoff.py [-h] -i INFILE -o OUTFILE -t TOTAL [-p PVAL]
    
    Calculate Poisson p value
    
    optional arguments:
     -h, --help show this help message and exit
     -i INFILE, --input INFILE
     files with total coverage
     -o OUTFILE, --output OUTFILE
     output
     -t TOTAL, --total TOTAL
     Total effective length
     -p PVAL, --pval PVAL

Input file should have 2 columns, first one is the sample name and the second one is the total usable reads count of that sample.

Total effective length is the genome length.

This program will calculate Poisson p value for coverage 1 to 10, and the smallest reads to reach the p value set by user.

genepredTodetailBed.py

This program takes [Genepred Table Format](https://genome.ucsc.edu/FAQ/FAQformat.html#format9) and make it into detailed BED format required by the annotation script.

Usage:

    python [genepredTodetailBed.py](http://genepredtodetailbed.py) genepred_input > detailbed_outPK 
     ��K�X 
�=  �=  .                 PIPE-CLIP2-b4ab87b14384438eab4977dc003b1c30.mdPK      \   E>   