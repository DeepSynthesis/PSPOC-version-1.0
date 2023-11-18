#!/usr/bin/env python
# -*- coding: utf-8 -*-

#XGBoost�㷨���ع�����

#����ѵ��������X_train
#����ѵ����ǩy_train

#���ѵ���õ�ģ��

import pandas as pd
import numpy as np
import xgboost as xgb

def L701(X_train, y_train):

	xgb0 = xgb.XGBRegressor(
		learning_rate= 0.04, 
		n_estimators= 1200, 
		max_depth= 10, 							
		min_child_weight= 7, 
		subsample= 0.8, 
		colsample_bytree= 0.6, 
		gamma= 0.1, 
		reg_alpha= 2, 
		reg_lambda= 0.5,
		verbosity= 2,
		n_jobs=-1,
		seed=42
		)
													
	xgb0.fit(X_train, y_train)

	return xgb0
	
	
if __name__ == '__main__':
	
	Data_X=pd.read_csv('fp_hspoc_v1.0_31_1.9_pKa_infp_20220407.csv',index_col=0)
	Data_Y=pd.read_csv('pKa_infp_20220407.csv')
	
	X=Data_X
	Y=Data_Y['pKa']
	
	