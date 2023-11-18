#!/usr/bin/python
# -*- coding: utf-8 -*-

#基于原子xyz坐标，将所有原子与指定原子按距离从小到大排序
##指定原子与自身距离为0

#输入待计算分子的SMILES
#目标的原子基于rdkit的索引Index
#输入pandas DataFrame,包含'SMILES','Atom_Index',‘x','y','z'索引

#输出原子索引的python list，顺序从小到大

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def L403m(mol,Index):
	
	atoms=mol.GetAtoms()
	distance_matrix=AllChem.Get3DDistanceMatrix(mol)
	
	Distance_List=[]
	for atom in atoms:
		Distance_List.append(distance_matrix[Index,atom.GetIdx()])
		
	Order_List=np.array(Distance_List).argsort()
	
	return Order_List
	
