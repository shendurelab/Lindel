#!/usr/bin/env python 
# Author: Will Chen

#System tools 
import os,sys,csv,re
from re import findall
import collections
import itertools
from collections import Counter
import pickle as pkl

# Usefule modules
import numpy as np
from math import exp
from random import randint
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord

# Umitools
import pyximport
pyximport.install(build_in_temp=False)
from umi_tools._dedup_umi import edit_distance

# Load pre-defined functions
from Seq_analysis import *

workdir = sys.argv[1] 
sample  = sys.argv[2] 
ref_l   = int(sys.argv[3]) # The length of reference. will fix this later
os.chdir(workdir)
workdir += '/' if workdir[-1]  != '/' else ''

#Create table
ref_table={}
ref = open('MH_design.txt')
for line in ref:
    if 'guide' not in line:
        guide = line.rstrip('\n').split('\t')[0]
        if guide not in ref_table:
            ref_table[guide]=0


#Create table
R1Names,R1Reads = read_fasta(workdir+sample+'_needleall.fasta')
table={}
for i in range(1,len(R1Reads),2):
    read = R1Reads[i].rstrip('-').lstrip('-')
    id1 = R1Reads[i].find(read)
    ref  = R1Reads[i-1][id1:id1+len(read)]
    if len(read)>250: # need to change depending on the length of the read
        table[R1Names[i]]      = [None]*15
        table[R1Names[i]][0:2] = read,ref
        if id1 < 50:
            table[R1Names[i]][5] = 'mh1'
        elif ref_l-10 < id1 < ref_l+30:
            table[R1Names[i]][5] = 'mh2'
        elif ref_l*2-10 < id1 < ref_l*2+30:
            table[R1Names[i]][5] = 'mh3'
        elif ref_l*3 < id1:
            table[R1Names[i]][5] = 'mh4'

temp_matrix = np.array([s for s in list(table.values()) if s[5] !=None and s[5]!='mh1' and 'W' not in s[1]])

# Step 1. identify potential guides and label template switch reads
for i in range(len(temp_matrix)):
    read  = temp_matrix[i,0]
    ref   = temp_matrix[i,1]
    s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
    id3 = temp_matrix[i,1][100:].find('NN')+100
    id1,id2,id4 = temp_matrix[i,1].find('NN'),temp_matrix[i,1].find('NG'),temp_matrix[i,1].find('NC')
    temp_matrix[i,2],temp_matrix[i,3] = temp_matrix[i,0][id1:id2+1], temp_matrix[i,0][id3:id4+1]
    if '---' not in read and '---' not in ref:
        if temp_matrix[i,2]!=temp_matrix[i,3] and '-' not in temp_matrix[i,1][id3-10:id3+30]:
            s1,s2 = temp_matrix[i,2].replace('-',''),temp_matrix[i,3].replace('-','')
            alignments = pairwise2.align.globalms(s1, s2,5,-4,-5,-0.5)
            if alignments and alignments[0][-1]>20:
                temp_matrix[i,-2] = 'template switch'
    temp_matrix[i,8] = id3-id1
    temp_matrix[i,0:2] = temp_matrix[i,0][id1:-20],temp_matrix[i,1][id1:-20]

# Step 2. Correct the guide RNA to the designed gRNA sequences
guides = {}
for guide in temp_matrix[:,2]:
    if guide:
        try:
            guides[guide] += 1
        except KeyError:
            guides[guide] = 1
guide_ref={}
for g in guides:
    if g in ref_table:
        guide_ref[g] = [g,0]
    else:
        guide_ref[g] = min_edit(ref_table,g)

# Step 3. Replace degenerate bases with sequence and re-aligen
for i in range(len(temp_matrix)):
    if temp_matrix[i,2]:
        guide,dis = guide_ref[temp_matrix[i,2]]
        if dis <=3 and temp_matrix[i,-2]!= 'template switch':
            temp_matrix[i,6:8] = guide,dis
            design = temp_matrix[i,5]
            if design == 'mh2':
                temp_matrix[i,1] = temp_matrix[i,1].replace('NNNNNNNNNNNNNNNNNNNN',guide).replace('NNNNNN',guide[11:17])
            elif design == 'mh3':
                temp_matrix[i,1] = temp_matrix[i,1].replace('NNNNNNNNNNNNNNNNNNNN',guide).replace('NNNN',guide[13:17])
            elif design == 'mh4':
                temp_matrix[i,1] = temp_matrix[i,1].replace('NNNNNNNNNNNNNNNNNNNN',guide).replace('NN',guide[15:17])
            temp_matrix[i,8] = temp_matrix[i,1].replace('-','')[100:].find(guide) +100

for i in range(len(temp_matrix)):
    read = temp_matrix[i]
    if read[7] != None and temp_matrix[i,-2]!= 'template switch':
        s1, s2 = read[0].replace('-',''),read[1].replace('-','')
        alignments = pairwise2.align.globalms(s1, s2,5,-4,-13,-0.5,penalize_end_gaps=False) 
        if alignments:
            idx = alignments[-1][1].find(read[6][0:10])
            temp_matrix[i][0:2] = alignments[-1][0][idx:],alignments[-1][1][idx:] # Here are we are taking the last alignment in the array, in which the deletion is right aligned. 
            temp_matrix[i,8] = temp_matrix[i,1].replace('-','')[100:].find(read[6]) +100
            idx2 = alignments[-1][1][100:].find(read[6][0:10])+100
            temp_matrix[i][2] = alignments[-1][0][0:20]
            temp_matrix[i][3] = alignments[-1][0][idx2:idx2+20]
            temp_matrix[i][4] = edit_distance(temp_matrix[i][2].encode('UTF-8'),temp_matrix[i][3].encode('UTF-8'))
            if temp_matrix[i,8]==99:
                temp_matrix[i,7] = None

fname = workdir + sample
pkl.dump(temp_matrix,open(fname+'_temp_matrix.pkl', 'wb'))

# editing table
# Step 4. Assign indels 
for i in range(len(temp_matrix)):
<<<<<<< HEAD
    if temp_matrix[i,7] != None and temp_matrix[i,-2]!= 'template switch':
        read  = temp_matrix[i,0]
        ref   = temp_matrix[i,1]
        ss,cs = temp_matrix[i,8], temp_matrix[i,8] +17 # Define start site of the target
        s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
        if '-' in read[cs-3:cs+4]:
            if len(s1) == 2:
            # The data is kind of messy. some bad alignment will show up as deletion in read and insertion in ref.
                if len(s2)<6:
                    locs=[]
                    for pattern in s2[1:-1]:
                        if len(locs)==0:
                            temp = list(re.search(r"\b"+pattern+r"\b", read).span())
                            locs.append(temp)
                        else:
                            temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", read[locs[-1][1]:]).span())]
                            locs.append(temp)
                    y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))
                    for loc in locs:
                        x=set(range(loc[0],loc[1]))
                        if len(x) ==1:
                            if x.intersection(y2):
                                temp_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
                        else:
                            if x.intersection(y1):
                                temp_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
                else:
                    temp_matrix[i,9] = 'complex'
            else:
                temp_matrix[i,9] = 'complex'
        # Label insertion
        elif '-' in ref[cs-3:cs+4]:
            if max([len(x) for x in s2[0:-1]]) <= 1:
        # May have indels. only count for the read with 1 deletion.
=======
    read  = temp_matrix[i,0]
    ref   = temp_matrix[i,1]
    ss,cs = temp_matrix[i,8], temp_matrix[i,8] +17 # Define start site of the target
    s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
    if '-' in read[cs-3:cs+2]:
        if len(s1) == 2:
        # The data is kind of messy. some bad alignment will show up as deletion in read and insertion in ref.
            if len(s2)<6: # the array may have some in dels but still good editing.
>>>>>>> parent of a3227e6... Update
                locs=[]
                for pattern in s2[1:-1]:
                    if len(locs)==0:
                        temp = list(re.search(r"\b"+pattern+r"\b", read).span())
                        locs.append(temp)
                    else:
                        temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", read[locs[-1][1]:]).span())]
                        locs.append(temp)
                y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))
                for loc in locs:
                    x=set(range(loc[0],loc[1]))
                    if len(x) ==1:
                        if x.intersection(y2):
                            temp_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
                    else:
                        if x.intersection(y1):
<<<<<<< HEAD
                            temp_matrix[i,9:13] = 'ins',loc[0]-ss,loc[1]-loc[0],read[loc[0]:loc[1]]
            else:
                temp_matrix[i,9] = 'complex'
        else:
            if '-' not in read and '-' not in ref and temp_matrix[i,-1]!='1bp indels':
                temp_matrix[i,9] = 'WT' 
            else:
                temp_matrix[i,9] = 'Trash' 
                
=======
                            temp_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
    # Label insertion
    elif '-' in ref[cs-3:cs+2]:
        if max([len(x) for x in s2[0:-1]]) <= 1:
        # May have indels. only count for the read with 1 deletion.
            locs=[]
            for pattern in s1[1:-1]:
                if len(locs)==0:
                    temp = list(re.search(r"\b"+pattern+r"\b", ref).span())
                    locs.append(temp)
                else:
                    temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", ref[locs[-1][1]:]).span())]
                    locs.append(temp)
            y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))  
            for loc in locs:
                x=set(range(loc[0],loc[1]))
                if len(x) ==1:
                    if x.intersection(y2):
                        temp_matrix[i,9:13] = 'ins',loc[0]-ss,loc[1]-loc[0],read[loc[0]:loc[1]]
                else:
                    if x.intersection(y1):
                        temp_matrix[i,9:13] = 'ins',loc[0]-ss,loc[1]-loc[0],read[loc[0]:loc[1]]
    else:
        if '--' not in read and '-' not in ref:
            idx = temp_matrix[i,0].find(temp_matrix[i,2]) # check if there is any deletions in tracRNA
            temp_matrix[i,9] = 'WT' if '-' not in temp_matrix[i,0][idx+1:idx+105] else 'Trash' # idx may
        else:
            temp_matrix[i,9] = 'Trash'

>>>>>>> parent of a3227e6... Update
# Step 5 take good reads
final_matrix = np.array([read for read in temp_matrix if read[7]!=None and read[9]!=None and read[9]!='Trash'\
                        and read[-2]!='template switch'])

# Step 6 For 1bp deletion, alignment index changed
for i in range(len(final_matrix)):
    if final_matrix[i,9]=='del' and final_matrix[i,11]==1:
        read  = final_matrix[i,0]
        ref   = final_matrix[i,1]
        s1, s2 = read.replace('-',''),ref.replace('-','')
        ali = pairwise2.align.globalms(s1, s2,5,-4,-13,-0.3)
        if len(ali)>1:
            final_matrix[i,0:2] = ali[0][0],ali[0][1]
            read  = final_matrix[i,0]
            ref   = final_matrix[i,1]
            ss = final_matrix[i,1].find(final_matrix[i,6][0:10])
            cs = ss + 17 # Define start site of the target
            s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
            locs=[]
            for pattern in s2[1:-1]:
                if len(locs)==0:
                    temp = list(re.search(r"\b"+pattern+r"\b", read).span())
                    locs.append(temp)
                else:
                    temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", read[locs[-1][1]:]).span())]
                    locs.append(temp)
            y = set(range(cs-2,cs+2))
            for loc in locs:
                x=set(range(loc[0],loc[1]))
                if x.intersection(y):
                    final_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]

# Step 7 Identify reads with both insertion and deletions, will show different alignment with different parameters 
for i in range(len(final_matrix)):
    read = final_matrix[i]
    if read[7] != None:
        s1, s2 = read[0][100:-15].replace('-',''),read[1][100:-15].replace('-','')
        ali1 = pairwise2.align.globalms(s1, s2,5,-4,-20,-0.5)
        ali2 = pairwise2.align.globalms(s1, s2,5,-4,-10,-0.4)
        if ali1[-1][0] != ali2[-1][0] or ali1[-1][1] != ali2[-1][1]:
<<<<<<< HEAD
            final_matrix[i][9] = 'complex'

for i in range(len(final_matrix)):
    if final_matrix[i,9]!='WT' and final_matrix[i,9]!='complex' and final_matrix[i,10]>17 and final_matrix[i,11]>1:
        read  = final_matrix[i,0]
        ref   = final_matrix[i,1]
        s1, s2 = read[30:].replace('-',''),ref[30:].replace('-','')
        ali = pairwise2.align.globalms(s1, s2,5,-4,-13,-0.3)
        if len(ali)>1:
            shift = final_matrix[i,10]-17
            if len(ali)>shift:
                final_matrix[i,0:2] = read[0:30]+ali[-1-shift][0],ref[0:30]+ali[-1-shift][1]
            elif len(ali)<=shift:
                final_matrix[i,0:2] = read[0:30]+ali[0][0],ref[0:30]+ali[0][1]
            read  = final_matrix[i,0]
            ref   = final_matrix[i,1]
            ss = final_matrix[i,1].find(final_matrix[i,6][0:10])
            cs = ss +17 # Define start site of the target
            s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
            if '-' in read[cs-3:cs+2]:
                locs=[]
                for pattern in s2[1:-1]:
                    if len(locs)==0:
                        temp = list(re.search(r"\b"+pattern+r"\b", read).span())
                        locs.append(temp)
                    else:
                        temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", read[locs[-1][1]:]).span())]
                        locs.append(temp)
                y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))
                for loc in locs:
                    x=set(range(loc[0],loc[1]))
                    if len(x) ==1:
                        if x.intersection(y2):
                            final_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
                    else:
                        if x.intersection(y1):
                            final_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
    # Label insertion
            elif '-' in ref[cs-3:cs+2]:
                locs=[]
                for pattern in s1[1:-1]:
                    if len(locs)==0:
                        temp = list(re.search(r"\b"+pattern+r"\b", ref).span())
                        locs.append(temp)
                    else:
                        temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", ref[locs[-1][1]:]).span())]
                        locs.append(temp)
                y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))  
                for loc in locs:
                    x=set(range(loc[0],loc[1]))
                    if len(x) ==1:
                        if x.intersection(y2):
                            final_matrix[i,9:14] = 'ins',loc[0]-ss,loc[1]-loc[0],read[loc[0]:loc[1]],read[loc[0]-1]
                    else:
                        if x.intersection(y1):
                            final_matrix[i,9:14] = 'ins',loc[0]-ss,loc[1]-loc[0],read[loc[0]:loc[1]],read[loc[0]-1]
            else:
                final_matrix[i,9]='complex'
=======
            final_matrix[i][-1] = 'realigned'


# For 1bp deletion, alignment index changed
for i in range(len(final_matrix)):
    if final_matrix[i,9]=='del' and final_matrix[i,11]==1:
        read  = final_matrix[i,0]
        ref   = final_matrix[i,1]
        s1, s2 = read.replace('-',''),ref.replace('-','')
        ali = pairwise2.align.globalms(s1, s2,5,-4,-9.5,-0.3)
        final_matrix[i,0:2] = ali[0][0],ali[0][1]
        read  = final_matrix[i,0]
        ref   = final_matrix[i,1]
        ss = final_matrix[i,1].find(final_matrix[i,6][0:10])
        cs = ss + 17 # Define start site of the target
        s1,s2 = re.split('[AGCTN]+',ref), re.split('[AGCTN]+',read)
        if '-' in read[cs-3:cs+2]:
            if len(s1) == 2:
        # The data is kind of messy. some bad alignment will show up as deletion in read and insertion in ref.
                if len(s2)<6: # the array may have some in dels but still good editing.
                    locs=[]
                    for pattern in s2[1:-1]:
                        if len(locs)==0:
                            temp = list(re.search(r"\b"+pattern+r"\b", read).span())
                            locs.append(temp)
                        else:
                            temp = [s+locs[-1][1] for s in list(re.search(r"\b"+pattern+r"\b", read[locs[-1][1]:]).span())]
                            locs.append(temp)
                    y1,y2  = set(range(cs-5,cs+5)),set(range(cs-2,cs+2))
                    for loc in locs:
                        x=set(range(loc[0],loc[1]))
                        if len(x) ==1:
                            if x.intersection(y2):
                                final_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
                        else:
                            if x.intersection(y1):
                                final_matrix[i,9:12] = 'del',loc[0]-ss,loc[1]-loc[0]
>>>>>>> parent of a3227e6... Update

# Stats about the percentage 
ts = len([s for s in temp_matrix if s[-2] =='template switch'])


editing_matrix = np.array([s for s in final_matrix if s[9] == 'del' or s[9]=='ins'])
edit = len(editing_matrix)/1./len(final_matrix) 


fname = workdir + sample
f0 = open(fname +'.editstat','w')
f0.write("\t Counts \n" \
         "Cells(total):\t" + str(len(table)) +'\n' \
         "Cells(designed):\t" + str(len(temp_matrix)) + '\n' \
         "Template switch ratio:\t"  +"{:.2f}%".format(100.*ts/len(temp_matrix)) + '\n'\
         "Cells for analysis:\t" + str(len(final_matrix)) + '('+ "{:.2f}%".format(100.*len(final_matrix)/len(temp_matrix)) + ')' + '\n'\
         "Editing efficiancy:\t" + str(edit)
        )
f0.close()
pkl.dump(final_matrix,open(fname+'_final_matrix.pkl', 'wb'))
