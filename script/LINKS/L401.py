#!/usr/bin/python
# -*- coding: utf-8 -*-

#����rdkit���ɷ��Ӹ���ԭ�ӵ�gasteiger���

#�����������ӵ�SMILES

#���pandas DataFrame��ʽ
##ÿ������indexΪSMILES_GasteigerCharge_AtomIndex
##������ΪSMILES,Atom_Index,Gasteiger_Charge

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def L401(SMILES,mol):	
	
	df=pd.DataFrame({
		'SMILES':[],
		'Atom_Index':[],
		'Gasteiger_Charge':[]
		})
	
	atoms=mol.GetAtoms()
	
	AllChem.ComputeGasteigerCharges(mol)
	
	for atom in atoms:
		
		Atom_Index=atom.GetIdx()
		charge=atom.GetDoubleProp('_GasteigerCharge')
		
		df.loc[str(SMILES)+'_GasteigerCharge_'+str(Atom_Index),['SMILES','Atom_Index','Gasteiger_Charge']]=[SMILES,Atom_Index,charge]
	
	return df
	
	
if __name__ == '__main__':	
	
	Data=pd.read_csv('D:\LSY\pKa_infp_20220407.csv')
	
	df_list=[]
	for smi in set(Data['SMILES'].values.tolist()):
		
		print('\nGenerate:\t\t'+smi)
		df=L401(smi)
		df_list.append(df)
		print('Done:\t\t'+smi+'\n')
		
	DF=pd.concat(df_list,axis=0,join='outer')
	DF.to_csv('D:\LSY\charge_infp_20220407.csv')