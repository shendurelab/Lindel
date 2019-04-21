
"""Functions used in this pipeline"""

def reverse_complement(seq):
    """This function returns the reverse_complement sequence of the input sequence
    from 3' to 5' """
    complement = {'A':'T', 'C':'G','G':'C', 'T':'A', 'N':'N'}
    rcomp_seq = ''
    for base in seq:
        rcomp_seq = complement[base] + rcomp_seq   
    return rcomp_seq

def reverse_seq(seq):
    """This function returns the reverse sequence of the input sequence from 3' to 5' """
    r_seq = ''
    for base in seq:
        r_seq = base + r_seq   
    return r_seq

def read_fastq(filename):
    rnames, reads = [], []
    entries=0
    for line in open(filename, 'r'):
        if entries%4 == 0: 
            rnames.append(line.rstrip().split(' ')[0][1:])
        elif entries%4 == 1:
            reads.append(line.rstrip())
        entries += 1
    return rnames, reads

def read_fasta(filename):
    """Read the aligned fasta file"""
    rnames, reads = [], []
    with open(filename,'r') as FASTA:
        for line in FASTA:
            if line[0] == '>':
                rnames.append(line.rstrip().split(' ')[0][1:])
                try:
                    reads.append(seq)
                except NameError:
                    pass
                seq=''
            else:
                seq += line.rstrip()
        reads.append(seq)
    FASTA.closed
    return rnames, reads

def min_edit(array,seq):
    m = 20
    ref=''
    for s in array.keys():
        dis = edit_distance(s.encode('UTF-8'),seq.encode('UTF-8'))
        if dis <= m:
            m = dis
            ref = s
    return [ref,m]

def seq_com(seq,ref,length):
    idx1,idx2 = 0,0
    array = np.zeros(length)
    while idx2 < length:
        try:
            if ref[idx1] == seq[idx1]:
                idx1 += 1
                idx2 += 1
            else:
                if ref[idx1] != '-' and seq[idx1] != '-': # mismatch
                    array[idx2] = 1
                    idx1 += 1
                    idx2 += 1
                elif ref[idx1] != '-' and seq[idx1] == '-': # Deletion
                    array[idx2] = 2
                    idx1 += 1
                    idx2 += 1                
                elif ref[idx1] == '-' and seq[idx1] != '-': # insertion
                    array[idx2-1] = 3
                    idx1 += 1
        except IndexError:
            break
    return array

def seq_count(sample):
    count = {}
    for read in sample:
        seq = read[6]
        if seq in count:
            count[seq] +=1
        else:
            count[seq] =1
    return count

def indel_count(sample):
    count = {}
    for read in sample:
        seq = read[6]
        if read[9]=='del':
            if seq in count:
                count[seq] += 1
            else:
                count[seq] = 1
        elif read[9] =='ins':
            if seq in count:
                count[seq] += 1
            else:
                count[seq] = 1
    return count

def parse_patterns(pattern):
    p = re.split('N',pattern)
    umi_length = pattern.count('N')
    return [p[0],p[-1]],umi_length


def ep_counter(data):
    ep ={}
    for sample in data:
        for read in sample:
            seq = read[6]
            if seq in ep:
                if read[9] == 'WT':
                    ep[seq][read[9]]+=1
                else:
                    key = str(read[10]) + '+' + str(read[11]) if read[9] == 'del' else str(read[10]) + '+' + str(read[11]) + '+' + str(read[12])
                    try:
                        ep[seq][read[9]][key] +=1
                    except KeyError:
                        ep[seq][read[9]][key] = 1
            else:
                ep[seq] = {'del':{},'ins':{},'WT':0}
                if read[9] != 'WT' :
                    key = str(read[10]) + '+' + str(read[11]) if read[9] == 'del' else str(read[10]) + '+' + str(read[11]) + '+' + str(read[12])
                    ep[seq][read[9]][key] = 1
                elif read[9] == 'WT':
                    ep[seq]['WT'] += 1
    return ep


def ep_counter_v2(data):
    ep ={}
    for read in data:
        seq = read[6]
        if seq in ep:
            if read[9] == 'WT':
                ep[seq][read[9]]+= read[-1]
            else:
                key = str(read[10]) + '+' + str(read[11]) if read[9] == 'del' else str(read[10]) + '+' + str(read[11]) + '+' + str(read[12])
                try:
                    ep[seq][read[9]][key] += read[-1]
                except KeyError:
                    ep[seq][read[9]][key] = read[-1]
        else:
            ep[seq] = {'del':{},'ins':{},'WT':0}
            if read[9] != 'WT' :
                key = str(read[10]) + '+' + str(read[11]) if read[9] == 'del' else str(read[10]) + '+' + str(read[11]) + '+' + str(read[12])
                ep[seq][read[9]][key] = read[-1]
            elif read[9] == 'WT':
                ep[seq]['WT'] += read[-1]
    return ep

def cal_entropy(array):
    t = sum(array) + 0.
    p = [x/t for x in array if x!=0]
    entro = -sum([np.log2(x)*x for x in p if x!=0])
    return entro

def read_editing_file(filename):
    data = []
    for line in open(filename, 'r'):
        if 'Indel' not in line:
            read = [None]*15
            row = line.rstrip('\n').split('\t')
            read[0:2] = row[-2:]
            read[6] = row[0]
            read[-1] = int(row[1])
            read[9:13] = row[2:6]
            data.append(read)
    return np.array(data)

def entro_cal(sample,cfilter):
    entropy = {}
    for seq in sample:
        if seq in cfilter:
            cwt, cdel,cins = sample[seq]['WT'],list(sample[seq]['del'].values()),list(sample[seq]['ins'].values())
            entropy[seq] = [cal_entropy([cwt] + cdel+cins),
                             cal_entropy(cdel+cins)]
    return entropy