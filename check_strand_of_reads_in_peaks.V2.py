#!/bin/python


'''
Author: Kikuye Koyano
Email: kkoyano@ucla.edu
Date: 11/01/2020

Description: Run for each sample. For each macs2 peak called, calculate the percent of reads sense or antisense. Concerned of DNA contamination for these samples. If most peaks are 50% plus and negative strand, then DNA contamination may be likely. If no contamination, expect a significant strand bias >98%. Note: Assumes BAM file is uniquely mapped.

Run chi squared test to test for significant strand bias.

input - 1) bed file contaning peak coordinates
        2) bam file of alinment
        3) output filename
        4) --sense or --antisense flag to indicate strandedness of library prep

output - 1) bed file of peak coordinates with several other stats such as: total_read_count,  alignment_gene_sense, alignment_gene_antisense, chisq


Usage: python check_strand_of_reads_in_peaks.py -i /u/home/k/kkoyano/gxxiao3/projects/long_csf_2020_03/analysis/peak_calling/macs2/data/broad_nomodel/fdr_5perc/fold_change_2/DKO1_Exo_High_longRNA_3.bed -b /u/home/k/kkoyano/gxxiao3/projects/long_csf_2020_03/data/hisat2/score-min_L,0,-0.2/DKO1_Exo_High_longRNA_3.uniq.sorted.bam -o /u/home/k/kkoyano/gxxiao3/projects/long_csf_2020_03/analysis/peak_calling/macs2/data/broad_nomodel/fdr_5perc/fold_change_2/DKO1_Exo_High_longRNA_3.chosen_anno.check_gene_strand.bed --antisense

--sense: the experiment is sense, so we expect that read1 will be antisense (opposite) strand to the gene
--antisense: the experiment is antisense, so we expect that read1 will be sense to the gene
input:
(1) bedtools intersect file of overlapped peaks and reads ex:
/u/home/k/kkoyano/gxxiao3/projects/long_csf_2020_03/analysis/peak_calling/macs2/data/broad_nomodel/fdr_5perc/fold_change_2/ DKO1_Exo_Low_longRNA_2.chosen_anno.bed
    19  56037808    56038099    DKO1_Exo_Low_longRNA_2.uniq_peak_13699  .   .   .   .
    19  56030657    56031294    DKO1_Exo_Low_longRNA_2.uniq_peak_13698  .   .   .   .
    5   134559913   134560113   DKO1_Exo_Low_longRNA_2.uniq_peak_24208  .   +   protein_coding  C5orf66
    5   134670114   134670938   DKO1_Exo_Low_longRNA_2.uniq_peak_24209  .   -   protein_coding  H2AFY
    5   133960883   133961480   DKO1_Exo_Low_longRNA_2.uniq_peak_24200  .   -   protein_coding  SAR1B
    5   134086353   134087015   DKO1_Exo_Low_longRNA_2.uniq_peak_24201  .   +   protein_coding  CAMLG
    5   134147336   134147554   DKO1_Exo_Low_longRNA_2.uniq_peak_24202  .   +   protein_coding  DDX46

(2) bam file
/u/home/k/kkoyano/gxxiao3/projects/long_csf_2020_03/data/hisat2/score-min_L,0,-0.2/DKO1_Exo_Low_longRNA_2.uniq.sorted.bam
samtools view DKO1_Exo_Low_longRNA_2.uniq.sorted.bam | head
    SRR8507532.1279190  99  1   10005   60  100M    =   10010   105 CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA    BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFFFFFFFFFFFFFFFF<FFFFFFFFFFFBFFFFBB<FFFF/<FFF/F/FFF/B    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100    YS:i:0  YT:Z:CP NH:i:1
    SRR8507532.27411156 99  1   10005   60  99M =   10046   138 CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT BBBBBFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFF<F<FFFF/FFFFFBFFFFF/FFFFF/FFFFFBBFF/FFFF<<FB<FFF/BF<BB/FBFFB AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:99 YS:i:-5 YT:Z:CP NH:i:1
    SRR8507532.30464730 99  1   10005   60  100M    =   10022   107 CCTAACCCTGACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA    BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFFFFFFFFFF/F<FFFFFBBB<B<BFFBFF</FFFF/F    AS:i:-5 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:9A90   YS:i:0  YT:Z:CP NH:i:1

output:
for each peak count how many reads overlap the peak and what is the proportion of reads mapping on the positive strand vs. the negative strand.
header:
    chrom start end peak_name . strand gene_chrom gene_start gene_end gene_name . gene_strand total_count alignment_gene_sense alignment_gene_antisense chisq



Version 2:

1. add the flag --sense in the input. previously we only had --antisense, and assumed if no flag that it was sense. just specified.

2. handle the exception if the sequencing library has no reads for a specified gene, then handle zero division (lines) line 192


'''

import argparse
import os
import sys
import itertools
from subprocess import Popen, PIPE
import subprocess
import time
from datetime import datetime
import gzip
import pysam
from scipy.stats import chisquare

def __unicode__(self):
   return unicode(self.some_field) or u''

def get_arguments():

    parser = argparse.ArgumentParser(description=' \
        Usage: python check_strand_of_reads_in_peaks.py -i <sample_name>.bed -b <sample_name>.bam -o <sample_name>.check_reads_gene_sense.bed --antisense \n\
         script performs the following,\
    1) for each peak, determine how many of the reads are on the same strand as the annotated gene for the peak \
    2) output bed file that tells the proportion of reads mapping sense to the gene. \
    experiment is assumed to be sense (read2 is sense to the gene), use --antisense flag if experiment is antisense. also assumes paired end read. only looks at read1. \
    ')

    parser.add_argument('-i ', dest='input_file', type=str,
                    help='name of the <sample_name>.chosen_anno.bed file of the annotations and the peaks. from annotate_intersect_bed.plus_ncgenes.V2_rank.macs2_peaks.py')
    parser.add_argument('-b ', dest='bam_file', type=str,
                    help='name of sample bam file that is sorted')
    parser.add_argument('-o ', dest='output_file', type=str,
                    help='path of directory that contains input file'),

    parser.add_argument('--antisense ', dest='antisense', action='store_true',
                    help='if flag is up then the experiment is antisense, so read1 is sense to the gene, otherwise will assume that the experiment is sense and read2 is sense to the gene'),
    parser.add_argument('--sense ', dest='sense', action='store_true',
                    help='if flag is up then the experiment is sense, so read1 is antisense to the gene.')


    args = parser.parse_args()


    if len(sys.argv) != 3:
        parser.print_help()
        exit()

    fin = args.input_file
    bam_filename = args.bam_file
    output_file = args.output_file
    antisense_flag = args.antisense
    sense_flag = args.sense



    antisense = False
    if antisense_flag == True:
        antisense = True
        print "The experiment is antisense, read1 strand will be on the same strand as the gene"
    elif (sense_flag == True):
        antisense = False
        print "The experiment is sense, read2 strand will be on the same strand as the gene"
    else:
        print "No strandedness specified, assuming sense. "
    return fin, output_file, bam_filename, antisense

def assign_read_strand(read_flag):
    read_strand = ''
    if read_flag == 99:
        read_strand = '+' # read1 is on the + strand
    elif read_flag == 147: # read2 is on the - strand
        read_strand = '-'
    elif read_flag == 83: # read1 is on the - strand
        read_strand = '-'
    elif read_flag == 163: # read2 is on the + strand
        read_strand = '+'
    else: # read is
        print 'unknown flag {}'.format(read_flag)

    return read_strand

def check_peak_reads_gene_sense(input_file, bam_filename, antisense):
    # print antisense
    input_obj = ''
    if input_file.endswith(".gz"):
        input_obj = gzip.open(input_file, mode='rb')
    else:
        input_obj = open(input_file, mode='r')

    samfile = pysam.AlignmentFile(bam_filename, "rb")

    output_peak_file_lines = []
    # go through each peak and find the reads that overlap with those regions
    for peak in input_obj:
        peak = peak.strip().split("\t")
        chrom, p_start, p_end, p_name, g_strand, g_biotype, g_name, g_feature = peak[0], int(peak[1]), int(peak[2]), peak[3], peak[5], peak[6], peak[7], peak[8]
        if g_strand == '.': # occurs when the peak is not annotated, assign the strand to the positive strand.
            g_strand = '+'
            g_biotype = 'unannotated'
            g_name = p_name
            g_feature = '.'
        # print chrom, p_start, p_end
        reads_within_peak = samfile.fetch(str(chrom), p_start, p_end)

        # initializing counts
        alignment_gene_sense, alignment_gene_antisense = 0, 0
        counted_readIDs = ()
        # count = 0
        for read in reads_within_peak:
            read_strand = ''
            readID, read_flag = read.query_name, read.flag


            # READ 1
            if read_flag == 99 or read_flag == 83:
                read1_strand = assign_read_strand(read_flag)
                # print read_flag, read_strand
                if antisense: # experiment is stranded and experiment is antisense
                    if read1_strand == g_strand:
                        alignment_gene_sense += 0.5
                    else: alignment_gene_antisense +=0.5
                else: # experiment is stranded and experiment is sense, read1 will have the opposite strand of the gene.
                    if read1_strand != g_strand:
                        alignment_gene_sense += 0.5
                    else: alignment_gene_antisense +=0.5
            # READ 2
            elif read_flag == 147 or read_flag == 163: # flags for read2
                read2_strand = assign_read_strand(read_flag)

                if antisense: # experiment is antisense, read2 should be opposite to gene
                    if read2_strand != g_strand:
                        alignment_gene_sense += 0.5
                    else: alignment_gene_antisense +=0.5
                else: # experiment is sense, read2 will have the same strand of the gene.
                    if read2_strand == g_strand:
                        alignment_gene_sense += 0.5
                    else: alignment_gene_antisense +=0.5
            else:
                continue
            # print readID, chrom, p_start, p_end, read_flag, read_strand


        total_count = alignment_gene_sense + alignment_gene_antisense
        chi_pvalue = "NA"
        if total_count == 0:
            alignment_gene_sense_perc = "NA"
        else:
            alignment_gene_sense_perc = alignment_gene_sense/float(total_count)
            obs = chisquare([alignment_gene_sense, alignment_gene_antisense], f_exp=[total_count/2, total_count/2])
            chi_pvalue = obs.pvalue


        peak.append(str(total_count))
        peak.append(str(alignment_gene_sense))
        peak.append(str(alignment_gene_sense_perc))
        peak.append(str(alignment_gene_antisense))


        peak.append(str(chi_pvalue))

        output_line = "\t".join(peak)


        output_peak_file_lines.append(output_line)
        # print peak
    output = "\n".join(output_peak_file_lines)

    return output

def main():
    start_time = datetime.now()
    input_file, output_file, bam_filename, antisense = get_arguments()
    output = check_peak_reads_gene_sense(input_file, bam_filename, antisense)
    end_time=datetime.now()

    header = ["chrom", "start", "end", "name", ".", "strand", "biotype", "gene_name", "gene_feature", "total_count", "reads_sense_to_gene","reads_sense_to_gene_perc", "reads_antisense_to_gene", "chisquare\n"]
    header_string = "\t".join(header)
    with open(output_file, 'w') as outf:
        outf.write(header_string)
        outf.write(output)

    print 'Duration: {}'.format(end_time-start_time)
    print 'Job Completed!'
    return

def test():
    input_file = ''
    output_file = 'test_output.txt'

    return

if __name__ == "__main__":
    # test()
    main()

