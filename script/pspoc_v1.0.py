#!/usr/bin/python
# -*- coding: utf-8 -*-

import subprocess
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from LINKS.L101 import L101
import PSPOC_web_v3 as pspoc
from sklearn.preprocessing import StandardScaler


def Predict(Data,fpname,modelname,traindafafilename):

    #参数表
    layer=20
    close_layer=1
    close_dis=4 	#单位A
    close_num=6
    PDBPATH='./pdb/'


    #读取标准化器
    with open('./mods/'+traindafafilename+'.pkl','rb') as f:
        st=pickle.load(f)
    #读取模型
    with open('./mods/'+modelname+'.pkl','rb') as f:
        xgb0=pickle.load(f)

    fp=pspoc.RUN(Data,PDBPATH,layer,close_layer,close_dis,close_num,'lower')
    fp.to_csv('./PredictPSPOC.csv')

    fp.replace(np.inf,0,inplace=True)

    standard_fp=pd.read_csv('./mods/'+fpname+'.csv',low_memory=False,index_col=0)

    FP=standard_fp._append(fp,ignore_index=True)
    FP=FP.fillna(0)

    Pre_X=FP.iloc[1:,:len(standard_fp.keys())].values
    Pre_X=st.transform(Pre_X)
    Pre_y=xgb0.predict(Pre_X)

    return Pre_y


if __name__ == '__main__':    

    PATH='./Pred/'
    datafilename='your data file'

    fpname='fp_pspoc_wv3.0b_20L1CL4D6DN_PKADpartmixb_standard'
    modelname='xgboost_AllTrain_fp_pspoc_wv3.0b_20L1CL4D6DN_PKADpartmixb_database'
    traindafafilename='StandardScalerfp_pspoc_wv3.0b_20L1CL4D6DN_PKADpartmixb_database'

    Data=pd.read_csv(PATH+datafilename+'.csv')
    Pre_pka=Predict(Data,fpname,modelname,traindafafilename)

    Data['Pre_pKa']=Pre_pka
    Data.to_csv(PATH+'After'+modelname+'Pred_'+datafilename+'.csv',index=0)
