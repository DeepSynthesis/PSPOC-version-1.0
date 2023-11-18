#!/usr/bin/python
# -*- coding: utf-8 -*-

#计算MACCSKeys分子指纹

#输入分子SMILES

#输出分子指纹list

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys

def L407(SMILES):
	
	mol=Chem.MolFromSmiles(SMILES)
	mol=Chem.AddHs(mol)
	
	fp=[x for x in MACCSkeys.GenMACCSKeys(mol)]
	
	return fp