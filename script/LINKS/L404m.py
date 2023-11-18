#!/usr/bin/python
# -*- coding: utf-8 -*-

#���ڷ��ӵ�SMILES������ָ��ԭ�Ӽ�ļ���
##��ԭ�Ӽ������ӣ������Ϊ0

#����SMILES
#���������ԭ��1�Ļ���rdkit������Index1
#���������ԭ��1�Ļ���rdkit������Index2

#�������float
##��������Index��ͬ��������-1

from rdkit import Chem

def L404m(mol,Index1,Index2):
	
	atoms=mol.GetAtoms()
	
	bonds=mol.GetBonds()
	for bond in bonds:
		if Index1 in [bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()] and Index2 in [bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]:
			rank=float(bond.GetBondTypeAsDouble())
			break
		else:
			rank=float(0)
	if Index1 == Index2:
		rank=-1
		
	return rank