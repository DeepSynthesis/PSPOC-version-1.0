#!/usr/bin/python
# -*- coding: utf-8 -*-

#����ԭ��xyz���꣬������ԭ����ָ��ԭ�Ӱ������С��������
##ָ��ԭ�����������Ϊ0

#�����������ӵ�SMILES
#Ŀ���ԭ�ӻ���rdkit������Index
#����pandas DataFrame,����'SMILES','Atom_Index',��x','y','z'����

#���ԭ��������python list��˳���С����

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
	
