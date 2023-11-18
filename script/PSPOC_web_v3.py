#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from Get_DHIndex_pdb import Get_DH
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors

from LINKS.L401m import L401m
from LINKS.L403m import L403m
from LINKS.L405m import L405m
from LINKS.L104 import L104

#########################################################################################################
#(pandas Series)
def Get_PSPOC(DataSeries,mol,layer,close_layer,close_dis,close_num,PDBPATH):
	
	atoms=mol.GetAtoms()
	PDB_ID=DataSeries['PDB_ID']
	Chain_ID=DataSeries['Chain ID']
	Residue_ID=DataSeries['Residue ID']
	Residue_Name=DataSeries['Residue']

	pdb_info=L104(PDBPATH,PDB_ID)
	DH_index=Get_DH(mol,pdb_info,DataSeries)

	atoms=mol.GetAtoms()
	AllChem.ComputeGasteigerCharges(mol)
	charge_dict=dict(   zip( [atom.GetIdx() for atom in atoms],[atom.GetDoubleProp('_GasteigerCharge') for atom in atoms] )   )
	
	df0=pd.DataFrame({'entry':['<'+Residue_Name+'>'+'<'+PDB_ID+'>'+'<'+str(Chain_ID)+'>'+'<'+str(Residue_ID)+'>''<'+str(DH_index)+'>']})


	

	H_Charge_Dict={'ASP':0.296269336,'HIS':0.311248158,'GLU':0.296267722,'LYS':0.344006605,'CYS':0.102304044,'TYR':0.29309207}
	df_H_Charge=pd.DataFrame({'H_Charge':[ H_Charge_Dict[Residue_Name] ]})

	Descriptor_AtomSets_List=L405m(mol,DH_index,layer)
	Descriptor={}
	
	bonds=mol.GetBonds()
	rank_matrix=np.zeros([len(atoms),len(atoms)])
	for bond in bonds:
		rank_matrix[bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]=float(bond.GetBondTypeAsDouble())
		rank_matrix[bond.GetEndAtomIdx(),bond.GetBeginAtomIdx()]=float(bond.GetBondTypeAsDouble())

	distance_matrix=AllChem.Get3DDistanceMatrix(mol)

	for i in range(layer):
		Atom_Index_List=list(Descriptor_AtomSets_List[i+1])
		for j in range(len(Atom_Index_List)):
			atom_index=int(Atom_Index_List[j])
			AtomicNum1=atoms[atom_index].GetAtomicNum()
			if AtomicNum1 != 1:
				Descriptor['Radius_'+str(i+1)+'_AtomicNum_'+str(j+1)]=[AtomicNum1]
				Descriptor['Radius_'+str(i+1)+'_Atom_Charge_'+str(j+1)]=[charge_dict[atom_index]]
				Descriptor['Radius_'+str(i+1)+'_BondOrder_'+str(j+1)]=max([rank_matrix[atom_index,x] for x in Descriptor_AtomSets_List[i]])
			
	Close_List=L403m(mol,DH_index)
	Close_list=[ x for x in Close_List if (atoms[int(x)].GetAtomicNum() in [7,8,9,15,16,17]) ]
	close_list=[ x for x in Close_List if ((distance_matrix[DH_index,x]) <= close_dis) and (int(x)!=int(DH_index)) and (x not in [atom.GetIdx() for atom in atoms[int(DH_index)].GetNeighbors()]) and (pdb_info[pdb_info['Atom_Index']==x]['Atom_Name'].values.tolist()[0] not in [' N  ',' O  ']) ]
	
	Order_AtomSets_List=L405m(mol,DH_index,6)
			
	for i in range(len(close_list)):
		if i+1 > close_num:
			break
		if atoms[int(close_list[i])].GetAtomicNum() != 1:
			Descriptor['Close_'+str(i+1)+'_AtomicNum']=[atoms[int(close_list[i])].GetAtomicNum()]
			Descriptor['Close_'+str(i+1)+'_Atom_Charge']=[charge_dict[int(close_list[i])]]
			#Descriptor['Close_'+str(i+1)+'_Residue_'+pdb_info[pdb_info['Atom_Index']==int(close_list[i])]['Residue_Name'].values.tolist()[0]]=[1]
			
			distance=distance_matrix[DH_index,close_list[i]]
			if close_list[i] in Order_AtomSets_List[2]:
				Descriptor['Close_'+str(i+1)+'_1/dis_DisNum3']=[1/distance]
			elif close_list[i] in Order_AtomSets_List[3]:
				Descriptor['Close_'+str(i+1)+'_1/dis_DisNum4']=[1/distance]
			elif close_list[i] in Order_AtomSets_List[4]:
				Descriptor['Close_'+str(i+1)+'_1/dis_DisNum5']=[1/distance]
			elif close_list[i] in Order_AtomSets_List[5]:
				Descriptor['Close_'+str(i+1)+'_1/dis_DisNum6']=[1/distance]
			else: 
				Descriptor['Close_'+str(i+1)+'_1/dis_DisNumMore']=[1/distance]
		
			Descriptor_Hbond_AtomSets_List=L405m(mol,close_list[i],close_layer)
			for j in range(close_layer):
				Hbond_Atom_Index_List=list(Descriptor_Hbond_AtomSets_List[j+1])
				for k in range(len(Hbond_Atom_Index_List)):
					hbond_atom_index=int(Hbond_Atom_Index_List[k])
					Descriptor['Close_'+str(i+1)+'Radius_'+str(j+1)+'_AtomicNum_'+str(k+1)]=[atoms[hbond_atom_index].GetAtomicNum()]
					Descriptor['Close_'+str(i+1)+'Radius_'+str(j+1)+'_Atom_Charge_'+str(k+1)]=[charge_dict[hbond_atom_index]]
					Descriptor['Close_'+str(i+1)+'Radius_'+str(j+1)+'_BondOrder_'+str(k+1)]=max([rank_matrix[hbond_atom_index,x] for x in Descriptor_Hbond_AtomSets_List[j]])
	
	df_Descriptor=pd.DataFrame(Descriptor)
	
	df=pd.concat([df0,df_H_Charge,df_Descriptor],axis=1,join='outer')
	
	return df

#RUN
# PDB_ID Chain_ID Residue_ID Residue
def RUN(Data,PDBPATH,layer,close_layer,close_dis,close_num,ID_mode):

	df_list=[]
	for i in range(len(Data)):
		PDB_ID=Data.loc[i]['PDB_ID']
		print('Do:'+PDB_ID+'\t'+Data.loc[i]['Chain ID']+str(Data.loc[i]['Residue ID']))
		if ID_mode=='lower':
			mol=Chem.MolFromPDBFile(PDBPATH+PDB_ID.lower()+'.pdb',sanitize=False)
		elif ID_mode=='':
			mol=Chem.MolFromPDBFile(PDBPATH+PDB_ID+'.pdb',sanitize=False)
		if 1 in [atom.GetAtomicNum() for atom in mol.GetAtoms()]:
			mol=Chem.RemoveHs(mol)			
		df=Get_PSPOC(Data.loc[i],mol,layer,close_layer,close_dis,close_num,PDBPATH)
		df_list.append(df)

	DF=pd.concat(df_list,axis=0,join='outer')
	DF=DF.fillna(0)

	return DF 

#########################################################################################################
if __name__ == '__main__':
	
#parameters									
	Data_filename='PKAD_database.csv'								
	PATH='./'
	MolPATH='./pdb/'	


	layer=20
	close_layer=1
	close_dis=4	
	close_num=6
	
	Out_PATH='./'
	info=str(layer)+'L'+str(close_layer)+'CL'+str(close_dis)+'D'+str(close_num)+'DN'
	Output_filename='fp_pspoc_wv3.0b'+'_'+info+'_'+Data_filename

#########################################################################################################
	Data=pd.read_csv(PATH+Data_filename)
	DF=RUN(Data,MolPATH,layer,close_layer,close_dis,close_num,'lower')
	DF.to_csv(Out_PATH+Output_filename,index=0)	
			
