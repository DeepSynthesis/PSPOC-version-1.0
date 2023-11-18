#!/usr/bin/python
# -*- coding: utf-8 -*-

#由mol文件输出距离矩阵

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def L409m(mol):
	
	distance_matrix=AllChem.Get3DDistanceMatrix(mol)
	
	return distance_matrix
