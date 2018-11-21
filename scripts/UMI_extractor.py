"""
Author: Will Chen
Date: 10/03/2017
Input: workdir fastq.gz 
Output: fastq file + stats
few things to input: 
1. UMIs: one UMI or two? and UMI patterns
2. Merge reads or not?
3. fwd read and rev read
4. 
Useage: python /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/scripts/UMI_extractor.py --workdir ** -S * -F * -R * -O * --bc_fwd * --bc_rev
        python /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/scripts/pipeline/UMI_extractor.py -S OT3-WT-E9 -F /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/data/OT/2017_09_18_OT3-WT_transient_transfect/OT3-WT-E9_1_S3_R2_001.fastq.gz -R /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/data/OT/2017_09_18_OT3-WT_transient_transfect/OT3-WT-E9_1_S3_R1_001.fastq.gz -O /net/shendure/vol10/projects/CRISPR.OT.Homing/nobackup/data/OT/2017_09_18_OT3-WT_transient_transfect/pipeline_output/ --bc_fwd 'NNNNNNNNNNNNNNNAGGA' --bc_rev 'NNNNNNNNNNNNNNNACGG' --umi_correction bc_fwd --threshold 3
"""
#System tools 
import gzip
import subprocess
import os,sys,csv,re
import collections
import itertools
from collections import Counter

# Usefule modules
from Seq_analysis import *
import numpy as np
import edlib

# Plot tools 
import pandas as pd

from optparse import OptionParser,OptionGroup
parser = OptionParser("%prog [options] ")
#parser.add_option("--workdir", dest="workdir")
parser.add_option("-S", "--sample",dest="sample")
parser.add_option("-F", "--fwdread", dest="forward")
parser.add_option("-R", "--revread", dest="reverse")
parser.add_option("-O", "--Outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("--bc_fwd", dest="bc_fwd", default=None)
parser.add_option("--bc_rev", dest="bc_rev", default=None)
parser.add_option("--umi_correction", dest="umi_correction", default="False")
parser.add_option("--threshold", dest="threshold", default=3)
(options, args) = parser.parse_args()

#options.workdir += '/' if options.workdir[-1] != '/' else ''
options.outdir  += '/' if options.outdir[-1]  != '/' else ''
f1 = options.forward #fwd read nextera
f2 = options.reverse #rev read truseq
outdir = options.outdir
threshold = int(options.threshold) if options.threshold else 3 # N read to count as a cell

if outdir != None and not os.path.isdir(outdir):
    os.makedirs(outdir)

# Generate regex function for barcode 1 and barcode 2.
bc_fwd = options.bc_fwd    
bc_rev = options.bc_rev
if bc_fwd != "None":
    p, umi1_length = parse_patterns(bc_fwd)
    regex = r"^" + p[0] + "([ATCG]{" + str(umi1_length) + "})(" + p[1] + ".*)"
    pattern1 = re.compile(regex)
else:
    pattern1 = re.compile(r".*")
    print ("This sample only has reverse barcode")
    
if bc_rev != "None":
    p, umi2_length = parse_patterns(bc_rev)
    regex = r"^" + p[0] + "([ATCG]{" + str(umi2_length) + "})(" + p[1] + ".*)"
    pattern2 = re.compile(regex)
else:
    pattern2 = re.compile(r".*")
    print ("This sample only has forward barcode")

#Read fwd read and rev read
print ("reading fwd_reads")
table = {}
count_line=0
cmd = "zcat "+f1
with os.popen(cmd) as handle: #always deal with the forward read first\
    handle = iter(handle)
    for line in handle:
        if line[0] == '@':
            count_line += 1
            read = next(handle).rstrip('\n')
            s = pattern1.match(read)
            if s:
                Rname = re.split('[@ ]+',line)[1]
                table[Rname] = [None]*15
                table[Rname][2] = Rname
                table[Rname][0] = s.groups()[0] if s.groups() else None
                table[Rname][3] = s.groups()[1] if s.groups() else s.group(0)
handle.close()
print ("reading rev_reads")
cmd = "zcat "+f2
with os.popen(cmd) as handle: #always deal with the forward read first\
    handle = iter(handle)
    for line in handle:
        if line[0] == '@':
            read = next(handle).rstrip('\n')
            s = pattern2.match(read)
            if s:
                Rname = re.split('[@ ]+',line)[1]
                try:
                    table[Rname][1] = s.groups()[0] if s.groups() else None
                    table[Rname][4] = s.groups()[1] if s.groups() else s.group(0)
                except KeyError:
                    pass
handle.close()

#UMI correction
UMI_count ={}
l_temp = 0 
for read in table.values():
    if read[4]!=None:
        l_temp +=1
        try:
            UMI_count[read[0]] += 1
        except KeyError:
            UMI_count[read[0]] = 1
print ("assigning cells")

cells = {}
l = len(table)
for s in table.values():
    if s[4]!=None:
        s[6] = s[0]
        if s[6] in cells:
            array =[]
            for i in range(len(cells[s[6]])):
                read = cells[s[6]][i]
                s1 = read[3] + read[4]
                s2 = s[3] + s[4]
                if s1==s2:
                    cells[s[6]][i][-1] += 1
                    array.append(0)
                elif edlib.align(s1,s2)['editDistance']/1./len(s2)<0.1: 
                # edlib.align only return edit distance, may introduce errors if there are indels. Need to fix this later
                    cells[s[6]][i][-1] += 1
                    array.append(0)
                else:
                    array.append(1)
            if np.all(array==1):
                cells[s[6]].append(s)
        else:
            if UMI_count[s[6]]>threshold:
                cells[s[6]] = [s]
                cells[s[6]][0][-1] = 1

print(("Filter cells with a threshold of " + str(threshold)))
final_matrix = np.array([cell[0] for cell in list(cells.values()) if cell[0][-1]>threshold])
collision_rate = len([1 for s in cells.values() if len(s)>1])/len(cells)*100
cread = sum(final_matrix[:,-1])

# Writing stats files and fastq files 
fname = outdir + options.sample

f0 = open(fname +'.readstats','w')
f0.write("\t Counts \n" \
         "Reads(total):\t" + str(count_line) +'\n' \
         "Reads(after regex):\t" + str(l_temp) + '\n' \
         "UMI count:\t" + str(len(UMI_count)) + '\n' \
         "UMI collision rate:\t" + "{0:.0000f}%".format(collision_rate) + '\n' \
         "Reads from cells:\t"   + str(cread) +'('+ "{0:.00f}%".format(cread/1./l_temp * 100)+')'+ '\n'\
         "Cells Captured:\t"     + str(len(final_matrix)))
f0.close()

f1 = open(fname + '.umiCounts','w')
f1.write("umi\t totalCount\t PASS\t sequence\n")
for key in cells:
    read = cells[key][0]
    string = read[6] + "\t" + str(read[-1]) +"\t"+ ('Pass' if read[-1] >threshold  else 'Fail') +"\t"+  read[3] + ";" + read[4] + "\n"
    f1.write(string)
f1.close()

#Write fastqfile
f2 = open(fname + '_fwd.fastq','w')
f3 = open(fname + '_rev.fastq','w')
for read in final_matrix:
    if 'None' not in [bc_fwd,bc_rev]:
        bc_string = read[6] + "+" + read[1] #if options.umi_correction == 'bc_fwd' else read[6] + "+" + read[0]
    else:
        bc_string = read[6]
    f2.write("@"+read[2]+' Fwd:' + bc_string + '\n' + read[3]+'\n' + '+\n' + "I"*len(read[3]) + '\n')
    f3.write("@"+read[2]+' Rev:' + bc_string + '\n' + read[4]+'\n' + '+\n' + "I"*len(read[4]) + '\n')
f2.close()
f3.close()