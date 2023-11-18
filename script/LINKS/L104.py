#!/usr/bin/python
# -*- coding: utf-8 -*-

#由PDB文件.pdb获得分子信息

#输入.pdb文件路径
#输入pdbID

#输出pandas DataFrame格式
##每行数据index为pdbID_xyz_AtomIndex
##列索引为PDB_ID，Chain_ID，Residue_ID，Residue_Name，Atom_Index，Atom_Type，Atom_Name，From，x，y，z

import numpy as np
import pandas as pd
from rdkit import Chem

def L104(pdb_filePATH,pdbID):
	
	df=pd.DataFrame({'PDB_ID':[],
									 'Chain_ID':[],
									 'Residue_ID':[],
									 'Residue_Name':[],
									 'Atom_Index':[],
									 'Atom_Type':[],
									 'Atom_Name':[]
									})	
	
	TER_num=0
	H_num=0
	Atom_Time=0
	with open(pdb_filePATH+pdbID.lower()+'.pdb') as f:								
		lines=f.readlines()
		for line in lines:
			if line[:6]=='ENDMDL':
				break
				
			if line[:4] == 'ATOM' or line[:6] == 'HETATM':
				
				if line[76:78]==' H':
					H_num=H_num+1
					continue
								
				if line[16] in ['B','C','D','E','F']:
					Atom_Time=Atom_Time+1
					continue
					
				Atom_Index=int(line[7:11])
				df.loc[pdbID+'_'+line[21]+'_'+str(int(line[23:26]))+'_'+line[17:20]+'_'+str(Atom_Index-TER_num-H_num-Atom_Time-1),
								['PDB_ID','Chain_ID','Residue_ID','Residue_Name','Atom_Index','Atom_Type','Atom_Name']
							]=[pdbID,line[21],int(line[23:26]),line[17:20],(Atom_Index-TER_num-H_num-Atom_Time-1),line[76:78],line[12:16]]


			if line[:3] == 'TER':
				TER_num=TER_num+1
			
	
	return df
	
	
