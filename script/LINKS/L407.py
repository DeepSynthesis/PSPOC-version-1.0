#!/usr/bin/python
# -*- coding: utf-8 -*-

#����MACCSKeys����ָ��

#�������SMILES

#�������ָ��list

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys

def L407(SMILES):
	
	mol=Chem.MolFromSmiles(SMILES)
	mol=Chem.AddHs(mol)
	
	fp=[x for x in MACCSkeys.GenMACCSKeys(mol)]
	
	return fp