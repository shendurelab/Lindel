#!/usr/bin/env python 
# Author: Will Chen
''' 
1. All functions are tested under python3.5 and python 3.6
2. Add Lindel folder to your python path.
3. y_hat is the prediction of all ~450 classes of indels <30bp.
4. fs is the frameshift ratio for this sequence.
5. Input should be 65bp and the clevage site should be at 30
'''

from Lindel.Predictor import * 
import pickle as pkl

weights = pkl.load(open("Model_weights.pkl",'rb'))
prerequesites = pkl.load(open('model_prereq.pkl','rb'))

seq = sys.argv[1].upper() #input your sequence here
filename = sys.argv[2]
try:
    y_hat, fs = gen_prediction(seq,weights,prerequesites)
    filename += '_frame_shift_ratio=' + str(fs)
    print (fs)
    rev_index = prerequesites[1]
	pred_freq = {}
	for i in range(len(y_hat)):
	    if y_hat[i]!=0:
	        pred_freq[rev_index[i]] = y_hat[i]
	pred_sorted = sorted(pred_freq.items(), key=lambda kv: kv[1],reverse=True)
	werite_file(seq,pred_sorted,pred_freq,filename)
except ValueError:
    return ('Error: No PAM sequence is identified.')
