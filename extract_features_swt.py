import wfdb
import os
import glob
import numpy as np
import matplotlib.pyplot as plt 
import random
import pickle
from biosppy.signals.ecg import christov_segmenter,engzee_segmenter,gamboa_segmenter,hamilton_segmenter,extract_heartbeats
from scipy.signal import resample,filtfilt,butter,freqz, medfilt,savgol_filter
import pywt
import sampen
import time
import sys
import multiprocessing as mp

def preprocess(sig):
    sig_resampled,_=wfdb.processing.resample_sig(sig,fs,fs_resampled)
    sig_resampled_mdf_stg1=medfilt(volume=sig_resampled,kernel_size=fs_resampled//2)
    sig_resampled_mdf_stg2=medfilt(volume=sig_resampled_mdf_stg1,kernel_size=fs_resampled-1)
    sig_resampled_mdf=sig_resampled-sig_resampled_mdf_stg2
    sig_resampled_mdf_sg=savgol_filter(x=sig_resampled_mdf,window_length=15,polyorder=3)
    return sig_resampled_mdf_sg    
    
def get_data(key):
    sig,fields=wfdb.srdsamp(recordname=key,channels=[1,2,5])
    sig_processed=[preprocess(record*2000) for record in sig.T]
    sig_processed=np.array(sig_processed)
    n=int(3.072*fs_resampled)
    n_segment=sig_processed.shape[1]//n
    segments_ch2=np.reshape(sig_processed[0,:n_segment*n],[n_segment,-1])
    segments_ch3=np.reshape(sig_processed[1,:n_segment*n],[n_segment,-1])
    segments_ch_avf=np.reshape(sig_processed[2,:n_segment*n],[n_segment,-1])
    features=[]
    for i in range(n_segment):
        coeff_ch2=pywt.swt(data=segments_ch2[i],wavelet='db5',level=6)
        coeff_ch3=pywt.swt(data=segments_ch3[i],wavelet='db5',level=6)
        coeff_ch_avf=pywt.swt(data=segments_ch_avf[i],wavelet='db5',level=6)
        coeff_2_ch2=pywt.swt(data=segments_ch2[i]**2,wavelet='db5',level=6)
        coeff_2_ch3=pywt.swt(data=segments_ch3[i]**2,wavelet='db5',level=6)
        coeff_2_ch_avf=pywt.swt(data=segments_ch_avf[i]**2,wavelet='db5',level=6)
        
        sen_d_1_3_3=sampen.sampen2(data=list(coeff_ch3[-3][1]),mm=2,r=0.2,normalize=True)[2][1]
        nse_a_2_3_3=compute_en(coeff_2_ch3[-3][0])/sum([compute_en(level[0]) for level in coeff_2_ch3])
        #sen_d_1_avf_2=sampen.sampen2(data=list(coeff_ch_avf[-2][1]),mm=2,r=0.2,normalize=True)[2][1]
        lee_d_2_2_2=sum([np.log2((x)**2) if x!=0 else 0 for x in coeff_2_ch2[-2][1]])        
        nse_a_2_3_1=compute_en(coeff_2_ch3[-1][0])/sum([compute_en(level[0]) for level in coeff_2_ch3])
        sen_d_1_2_2=sampen.sampen2(data=list(coeff_ch2[-2][1]),mm=2,r=0.2,normalize=True)[2][1]
        mds_d_1_2_2=np.median(np.abs(coeff_ch2[-2][1][1:]-coeff_ch2[-2][1][:-1]))*fs_resampled
        sen_d_1_avf_2=sampen.sampen2(data=list(coeff_ch_avf[-2][1]),mm=2,r=0.2,normalize=True)[2][1]
        lee_d_1_2_1=sum([np.log2((x)**2) if x!=0 else 0 for x in coeff_ch2[-1][1]])        
        mds_d_1_avf_1=np.median(np.abs(coeff_ch_avf[-1][1][1:]-coeff_ch_avf[-1][1][:-1]))*fs_resampled    
        
        feature=[sen_d_1_3_3,nse_a_2_3_3,sen_d_1_avf_2,lee_d_2_2_2,nse_a_2_3_1,sen_d_1_2_2,mds_d_1_2_2,
                 lee_d_1_2_1,mds_d_1_avf_1]
        features.append(np.array(feature))
        
    if 'Myocardial infarction' in fields['comments'][4]:
        label=1
    if 'Healthy control' in fields['comments'][4]:
        label=0
    return (key,np.array(features), label)
def compute_en(vec):
    return sum(vec**2)/len(vec)

path='D:\\Project\\PTB\\code' 
os.chdir(path)
fs=1e3
fs_resampled=250 
    
if __name__ == "__main__":
    data_dir=os.path.join('..','ptbdb')
    filepaths=list(set([os.path.splitext(fl)[0] for fl in glob.glob(os.path.join(data_dir,'*','*'))]))
    
    print('Getting filepaths...')
    keys_imi=[]
    keys_hc=[]
    for i,path in enumerate(filepaths):
        _,fields=wfdb.srdsamp(path)
        if 'Healthy control' in fields['comments'][4]:
            keys_hc.append(path)
        else:
            if 'Myocardial infarction' in fields['comments'][4]:
                if 'inferior' in fields['comments'][5]:
                    keys_imi.append(path)
        print('processed {}/{}'.format(i+1,len(filepaths)),end='\r')
    print('\n')

    pool = mp.Pool(processes=4)
    print('Computing features...')
    results=[]
    keys=keys_hc+keys_imi
    for i, item in enumerate(pool.imap_unordered(get_data,keys), 1):
        results.append(item)
        print('processed {}/{}'.format(i,len(keys)),end='\r')
    #print('\n')
    data={}
    for result in results:
        data[result[0]]=(result[1],result[2])

    with open(os.path.join('..','data','feature_imi_hc_9.bin'),'wb') as pfile:
        pickle.dump(file=pfile,obj=data,protocol=4)