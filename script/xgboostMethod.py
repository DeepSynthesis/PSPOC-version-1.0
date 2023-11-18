#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split,KFold
from sklearn.metrics import r2_score, classification_report, mean_absolute_error, mean_squared_error

import xgboost
import pickle


def load_data(data_name, fp_name, split_ratio, random_seed=42):

	data = pd.read_csv(data_name+'.csv') 
	fingerprint = pd.read_csv(fp_name+'.csv',low_memory=False,index_col=0)
	fingerprint = fingerprint.replace(np.inf,0) 

	y = data['pKa'].values
	X = fingerprint.values
	
	st = StandardScaler()

	X = st.fit_transform(X)
	X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=split_ratio, random_state=random_seed)
	return X_train, X_test, y_train, y_test, X, y, st
	
def xgb(X_train, y_train):
	title = r'Extra Tree Regressor'
	xgb1 = xgboost.XGBRegressor(learning_rate= 0.04, 
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
							seed=42)
	xgb1.fit(X_train, y_train)

	return xgb1
		

def plotResult(y, y_train, y_pred_train, split_ratio, fp, figure_title=''):	

	rms_train = (np.mean((y_train - y_pred_train)**2))**0.5
	rms_test = (np.mean((y_test - y_pred_test)**2))**0.5
	print('rms_train: %s; r^2_train: %s' % (round(rms_train,3), round(r2_score(y_train,y_pred_train),3)))
	print('rms_test: %s; r^2_test: %s' % (round(rms_test,3), round(r2_score(y_test,y_pred_test),3)))

	# plot
	plt.figure(figsize=(5,3))
	plt.scatter(y_train,y_pred_train, label = 'Train', c='deepskyblue',alpha = 0.5,s=20)
	plt.title(figure_title)
	plt.xlabel('Experimence')
	plt.ylabel('Prediction')
	x_start = min(y)-(max(y)-min(y))*0.1
	x_end = min(y)+(max(y)-min(y))*1.1
	plt.xlim((x_start,x_end))
	plt.ylim((x_start,x_end))
	plt.scatter(y_test,y_pred_test,c='navy', label='Test',alpha = 0.5,s=20)
	x1 = [x_start,x_end]
	y1 = x1
	plt.plot(x1,y1,c='lightcoral', alpha = 0.8)
	plt.legend(loc=4)
	#
	point_x1 = min(y)-(max(y)-min(y))*0.05
	#
	test1 = 'RMSE = '+str(round(rms_test,3))
	test2 = 'R$^2$ = '+str(round(r2_score(y_test,y_pred_test),3))
	point_y3 = (max(y)-min(y))*1+min(y)
	point_y4 = (max(y)-min(y))*0.9+min(y)
	plt.text(point_x1, point_y3,test1, weight = "light")
	plt.text(point_x1, point_y4,test2, weight = "light")
	#
	plt.savefig('xgboost_'+str(split_ratio)+'train_'+fp+'.png', dpi=300, bbox_inches='tight')
	return



def modelGeneration(xgb0, st, info):

	with open('xgboost_AllTrain_'+info+'.pkl','wb') as f:
		pickle.dumps(xgb0,f)
	with open('StandardScaler'+info+'.pkl','wb') as f:
		pickle.dumps(st,f)

	return



def kFold(X,y,k,fp_name,info='',random_seed=42):

	kf=KFold(n_splits=k,shuffle=True,random_state=random_seed)
	R2=[]
	RMSE=[]
	MAE=[]
	for train_index, test_index in kf.split(X):
		xgb_model=xgb(X[train_index],y[train_index])
		prediction=xgb_model.predict(X[test_index])
		actuals=y[test_index]
		R2.append(r2_score(prediction,actuals))
		MAE.append(mean_absolute_error(prediction,actuals))
		RMSE.append(np.sqrt(mean_squared_error(prediction,actuals)))

	print('R2_avg='+str(np.mean(R2))+'\n'+'MAE_avg='+str(np.mean(MAE))+'\n'+'RMSE_avg='+str(np.mean(RMSE)))
	with open('xgboost_'+info+str(k)+'FoldResult_'+fp_name+'.txt','w') as f:
		f.write('R2_avg='+str(np.mean(R2))+'\n'+str(R2)+'\n')
		f.write('MAE_avg='+str(np.mean(MAE))+'\n'+str(MAE)+'\n')
		f.write('RMSE_avg='+str(np.mean(RMSE))+'\n'+str(RMSE)+'\n')

	return np.mean(R2),np.mean(MAE),np.mean(RMSE)
		
def Contribution_Analyst(fp,xgb_model,info=''):

	fingerprint = pd.read_csv(fp+'.csv',low_memory=False,index_col=0)
	im=pd.DataFrame({'parameter':fingerprint.columns,'importance':xgb_model.feature_importances_})
	im=im.sort_values(by='importance',ascending=False)
	im.to_csv('xgboost_'+info+'contribution_'+fp+'.csv')

	return im

if __name__ == '__main__':

###############################################################################	
	data_name = 'PKAD_database'
	fp_name = ['fp_pspoc_wv3.0b_20L1CL4D6DN_PKAD_database',]
	figure_title=''
	split_ratio=0.8
###############################################################################

	for fp in fp_name:	

		X_train, X_test, y_train, y_test, X, y, st = load_data(data_name, fp, split_ratio)

		#输出模型
		xgb0 = xgboost(X, y)
		modelGeneration(xgb0,st,fp+'_'+data_name)
		

		xgb1 = xgb(X_train, y_train)
		y_pred_train = xgb1.predict(X_train)
		y_pred_test = xgb1.predict(X_test)
		plotResult(y, y_train, y_pred_train, split_ratio, fp, figure_title)
		kFold(X,y,5,fp)
		Contribution_Analyst(fp,xgb1)



