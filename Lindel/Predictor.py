import numpy as np
import scipy.sparse as sparse
import re
import json
import itertools

def gen_indel(sequence: str, cut_site: int) -> list:
    """
    Generates all possible unique indels and lists the redundant classes
    which will be combined after.
    """
    nt = ['A', 'T', 'C', 'G']
    up = sequence[:cut_site]
    down = sequence[cut_site:]
    dmax = min(len(up), len(down))
    unique_seq = {}

    # Deletions
    for dstart in range(1, cut_site + 3):
        for dlen in range(1, dmax):
            if cut_site - 2 < dstart + dlen < len(sequence):
                seq = sequence[:dstart] + sequence[dstart + dlen:]
                indel = sequence[:dstart] + '-' * dlen + sequence[dstart + dlen:]
                array = [indel, sequence, 13, 'del', dstart - 30, dlen, None, None, None]
                if seq not in unique_seq or dstart - 30 < 1:
                    unique_seq[seq] = array

    # Insertions
    for length in range(1, 3):
        for bases in itertools.product(nt, repeat=length):
            base_str = "".join(bases)
            seq = sequence[:cut_site] + base_str + sequence[cut_site:]
            indel = sequence[:cut_site] + '-' * length + sequence[cut_site:]
            array = [sequence, indel, 13, 'ins', 0, length, base_str, None, None]
            unique_seq[seq] = array

    uniq_align = label_mh(list(unique_seq.values()), 4)
    for read in uniq_align:
        if read[-2] == 'mh':
            read[-3] = [(read[4] - i, read[5]) for i in range(read[-1] + 1)]
    return uniq_align

def label_mh(sample,mh_len):
    '''Function to label microhomology in deletion events'''
    for k in range(len(sample)):
        read = sample[k]
        if read[3] == 'del':
            idx = read[2] + read[4] +17
            idx2 = idx + read[5]
            x = mh_len if read[5] > mh_len else read[5]
            for i in range(x,0,-1):
                if read[1][idx-i:idx] == read[1][idx2-i:idx2] and i <= read[5]:
                    sample[k][-2] = 'mh'
                    sample[k][-1] = i
                    break
            if sample[k][-2]!='mh':
                sample[k][-1]=0
    return sample


def create_feature_array(ft,uniq_indels):
    '''Used to create microhomology feature array 
       require the features and label 
    '''
    ft_array = np.zeros(len(ft))
    for read in uniq_indels:
        if read[-2] == 'mh':
            mh = str(read[4]) + '+' + str(read[5]) + '+' + str(read[-1])
            try:
                ft_array[ft[mh]] = 1
            except KeyError:
                pass
        else:
            pt = str(read[4]) + '+' + str(read[5]) + '+' + str(0)
            try:
                ft_array[ft[pt]]=1
            except KeyError:
                pass
    return ft_array


def onehotencoder(seq: str) -> np.ndarray:
    """Converts sequence to single and di-nucleotide one-hot encoding."""
    nt = ['A', 'T', 'C', 'G']
    l = len(seq)
    
    # Pre-calculate sizes for efficiency
    single_nt_size = len(nt) * l
    di_nt_size = len(nt) * len(nt) * (l - 1)
    total_size = single_nt_size + di_nt_size

    # Create a mapping from nucleotide/dinucleotide to index
    head_idx = {char + str(i): i * len(nt) + j for i in range(l) for j, char in enumerate(nt)}
    offset = single_nt_size
    for i in range(l - 1):
        for j, char1 in enumerate(nt):
            for k, char2 in enumerate(nt):
                head_idx[char1 + char2 + str(i)] = offset + i * (len(nt) * len(nt)) + j * len(nt) + k

    # Create the one-hot encoded array
    encode = np.zeros(total_size)
    for j, char in enumerate(seq):
        encode[head_idx[char + str(j)]] = 1.0
    for k in range(l - 1):
        encode[head_idx[seq[k:k+2] + str(k)]] = 1.0
        
    return encode

def create_label_array(lb,ep_freq,seq):
    lb_array = np.zeros(len(lb))
    for pt in ep_freq[seq]['del']:
        lb_array[lb[pt]] = ep_freq[seq]['del'][pt]
    for pt in ep_freq[seq]['ins']:
        lb_array[lb[pt]] = ep_freq[seq]['ins'][pt]
    return lb_array


def gen_prediction(seq,wb,prereq):
    '''generate the prediction for all classes, redundant classes will be combined'''
    pam = {'AGG':0,'TGG':0,'CGG':0,'GGG':0}
    guide = seq[13:33]
    if seq[33:36] not in pam:
        return ('Error: No PAM sequence is identified.')
    w1,b1,w2,b2,w3,b3 = wb
    label,rev_index,features,frame_shift = prereq
    indels = gen_indel(seq,30) 
    input_indel = onehotencoder(guide)
    input_ins   = onehotencoder(guide[-6:])
    input_del   = np.concatenate((create_feature_array(features,indels),input_indel),axis=None)
    cmax = gen_cmatrix(indels,label) # combine redundant classes
    dratio, insratio = softmax(np.dot(input_indel,w1)+b1)
    ds  = softmax(np.dot(input_del,w2)+b2)
    ins = softmax(np.dot(input_ins,w3)+b3)
    y_hat = np.concatenate((ds*dratio,ins*insratio),axis=None) * cmax
    return (y_hat,np.dot(y_hat,frame_shift))

def softmax(weights):
    return (np.exp(weights)/sum(np.exp(weights)))

def gen_cmatrix(indels,label): 
    ''' Combine redundant classes based on microhomology, matrix operation'''
    combine = []
    for s in indels:
        if s[-2] == 'mh':
            tmp = []
            for k in s[-3]:
                try:
                    tmp.append(label['+'.join(list(map(str,k)))])
                except KeyError:
                    pass
            if len(tmp)>1:
                combine.append(tmp)
    temp = np.diag(np.ones(557), 0)
    for key in combine:
        for i in key[1:]:
            temp[i,key[0]] = 1
            temp[i,i]=0    
    return (sparse.csr_matrix(temp))

def format_predictions(seq: str, pred_sorted: list, pred_freq: dict, output_type: str = 'json', fname: str = None):
    """Formats predictions into JSON or a file."""
    output_data = []
    ss = 13
    cs = ss + 17
    # Add the original sequence information first
    output_data.append({
        "Sequence": f"{seq[:30]} | {seq[30:60]}",
        "Frequency": "0",
        "Indels": ""
    })

    for pt, freq in pred_sorted:
        try:
            idx1, dl = map(int, pt.split('+'))
            indel_type = f"D{dl}  {idx1 - 30}"
            idx1 += cs
            idx2 = idx1 + dl
            if idx1 < cs:
                if idx2 >= cs:
                    s = f"{seq[:idx1]}{'-' * (cs - idx1)} | {'-' * (idx2 - cs)}{seq[idx2:]}"
                else:
                    s = f"{seq[:idx1]}{'-' * (idx2 - idx1)}{seq[idx2:cs]} | {seq[cs:]}"
            elif idx1 > cs:
                s = f"{seq[:cs]} | {seq[cs:idx1]}{'-' * dl}{seq[idx2:]}"
            else:
                s = f"{seq[:idx1]} | {'-' * dl}{seq[idx2:]}"
        except ValueError:
            idx1 = int(pt.split('+')[0])
            if pt != '3':
                bp = pt.split('+')[1]
                il = str(idx1)
                indel_type = f"I{il}+{bp}"
            else:
                bp = 'X'  # label any insertion >= 3bp as X
                il = '>=3'
                indel_type = f"I3+{bp}"
            s = f"{seq[:cs]} {bp}{' ' * (2 - len(bp))}{seq[cs:]}"
        
        output_data.append({
            "Sequence": s,
            "Frequency": f"{freq * 100:.2f}",
            "Indels": indel_type
        })

    if output_type == 'json':
        return json.dumps(output_data, indent=1)
    elif output_type == 'file' and fname:
        with open(fname, 'w') as f:
            for item in output_data:
                f.write(f"{item['Sequence']}\t{item['Frequency']}\t{item['Indels']}\n")
    else:
        return None

