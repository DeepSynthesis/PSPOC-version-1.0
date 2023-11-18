#!/usr/bin/python
# -*- coding: utf-8 -*-

#基于分子的SMILES，计算指定原子间的键级
##若原子间无连接，则键级为0

#输入SMILES
#输入待计算原子1的基于rdkit的索引Index1
#输入待计算原子1的基于rdkit的索引Index2

#输出键级float
##如果输入的Index相同，则会输出-1

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