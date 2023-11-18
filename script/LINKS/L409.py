#!/usr/bin/python
# -*- coding: utf-8 -*-

#依据L101中的xyz表输出距离矩阵的部分值

import numpy as np
import pandas as pd
from rdkit import Chem

def L409(mol,df):
	
	atoms=mol.GetAtoms()	
	total_num=len(atoms)
	distance_matrix=np.zeros([total_num,total_num])
	
	for i in range(total_num):
		for j in range(i):
		
		    xyz1=np.array(df[df['Atom_Index']==atoms[i].GetIdx()][['x','y','z']])
		    xyz2=np.array(df[df['Atom_Index']==atoms[j].GetIdx()][['x','y','z']])
		
		    distance_matrix[i,j]=np.linalg.norm(xyz1-xyz2)
		    distance_matrix[j,i]=np.linalg.norm(xyz2-xyz1)
	
	return distance_matrix
