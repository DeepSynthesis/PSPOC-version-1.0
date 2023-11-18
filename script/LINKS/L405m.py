#!/usr/bin/python
# -*- coding: utf-8 -*-

#������ԭ��Ϊ��0�㣬�ع��ۼ���������չ��������ÿһ���ԭ������
##����Hԭ��

#����rdkit mol
#��������ԭ��Index
#����չ������layer

#�����setΪԪ�ص�list
##list��ʽΪ[set(�����ԭ������(������ԭ��)),set([��һ��ԭ������]),set([�ڶ���ԭ������]),...]
##�б�������Ϊ����

from rdkit import Chem

def L405m(mol,Index,layer):
	
	atoms=mol.GetAtoms()
	
	ECSet_List=[set([Index])]
	
	for i in range(layer):
		Connect_Set=set([])
		for x in ECSet_List[i]:	
			neighbors=atoms[int(x)].GetNeighbors()	
			for atom in neighbors:
				Connect_Set.add(atom.GetIdx())
		ECSet_List.append(Connect_Set)
		
	for i in range(layer+1):
		for x in ECSet_List[:i]:
			ECSet_List[i]=ECSet_List[i]-x
		
	return ECSet_List