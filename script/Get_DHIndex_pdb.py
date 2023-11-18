#!/usr/bin/python
# -*- coding: utf-8 -*-

from rdkit import Chem
from LINKS.L104 import L104
from LINKS.L404m import L404m
import numpy as np
import pandas as pd



def Get_DH(mol,df,DataSeries):
	'''
	get index of the atom which directly connected with H atom.
	args:
		mol: mol object of protein molecular
		df: dataframe of origin data 
		DataSeries: series of origin dataframe
	output:
		DH_index: an int number of index 
	
	
	'''
	Residue_Data=df
	Residue_Data=Residue_Data[Residue_Data['PDB_ID']==DataSeries['PDB_ID']] 
	Residue_Data=Residue_Data[Residue_Data['Chain_ID']==DataSeries['Chain ID']] 
	Residue_Data=Residue_Data[Residue_Data['Residue_ID']==DataSeries['Residue ID']]
	Residue_Data=Residue_Data[Residue_Data['Residue_Name']==DataSeries['Residue']] 
							 
	DH_index=-1
	if DataSeries['Residue']=='ASP':
		OD1_index=int(Residue_Data[Residue_Data['Atom_Name']==' OD1']['Atom_Index'].values.tolist()[0])
		OD2_index=int(Residue_Data[Residue_Data['Atom_Name']==' OD2']['Atom_Index'].values.tolist()[0])
		CG_index=int(Residue_Data[Residue_Data['Atom_Name']==' CG ']['Atom_Index'].values.tolist()[0])
		if L404m(mol,OD1_index,CG_index)==1:
			DH_index=OD1_index
		elif L404m(mol,OD2_index,CG_index)==1:
			DH_index=OD2_index	
			
	elif DataSeries['Residue']=='CYS':
		DH_index=int(Residue_Data[Residue_Data['Atom_Name']==' SG ']['Atom_Index'].values.tolist()[0])
	
	elif DataSeries['Residue']=='GLU':
		OE1_index=int(Residue_Data[Residue_Data['Atom_Name']==' OE1']['Atom_Index'].values.tolist()[0])
		OE2_index=int(Residue_Data[Residue_Data['Atom_Name']==' OE2']['Atom_Index'].values.tolist()[0])
		CD_index=int(Residue_Data[Residue_Data['Atom_Name']==' CD ']['Atom_Index'].values.tolist()[0])
		if L404m(mol,OE1_index,CD_index)==1:
			DH_index=OE1_index
		elif L404m(mol,OE2_index,CD_index)==1:
			DH_index=OE2_index
	
	elif DataSeries['Residue']=='HIS':
		DH_index=int(Residue_Data[Residue_Data['Atom_Name']==' ND1']['Atom_Index'].values.tolist()[0])
	
	elif DataSeries['Residue']=='LYS':
		DH_index=int(Residue_Data[Residue_Data['Atom_Name']==' NZ ']['Atom_Index'].values.tolist()[0])
	
	elif DataSeries['Residue']=='TYR':
		DH_index=int(Residue_Data[Residue_Data['Atom_Name']==' OH ']['Atom_Index'].values.tolist()[0])
	
	return DH_index

if __name__ == '__main__':
	
	MolPATH='./pdb_data/'
	PDBPATH='./pdb_data/'
	datafile='test_set_WT+MT+aSN.csv'
	
	Data=pd.read_csv(datafile)
	DH_List=[]
	for i in range(len(Data)):
		PDB_ID=Data.loc[i]['PDB_ID']	
		Chain_ID=Data.loc[i]['Chain ID']
		Residue_ID=Data.loc[i]['Residue ID']
		Residue_Name=Data.loc[i]['Residue']
		mol=Chem.MolFromPDBFile(MolPATH+PDB_ID.lower()+'.pdb')
		print('\nGenerate:\t\t'+PDB_ID+'_'+Chain_ID+str(Residue_ID))
		atoms=mol.GetAtoms()
		df=L104(PDBPATH,PDB_ID)
		DH_List.append( Get_DH(mol,df,Data.loc[i]) )
		print('Done:\t\t\t'+PDB_ID+'_'+Chain_ID+str(Residue_ID)+'\n')
	DH_data=pd.DataFrame({'DH_Index':DH_List})
	DH_data.to_csv('DH_Index_test_set_WT+MT+aSN.csv')
		
		