#!/usr/bin/env python
#-*- coding: utf-8 -*-

# this python script is to search through one single .pdb file with multi models and save them in different files


# syntax: ./pdb_extract.py pdbfile.pdb
# output: pdbfile_1.pdb pdbfile_2.pdb ... pdbfile_n.pdb

from sys import argv,exit


#script name
script, filename = argv

#get filename
name = filename.split('.')[0]

#read lines
opener = open(filename, 'r')
lines = opener.readlines()
opener.close()

modelIndexNum = []
endIndexNum = []

for i in range(len(lines)):
    if lines[i].split()[0] == 'MODEL':
        modelIndexNum.append(i)
    elif lines[i].split()[0] == 'ENDMDL':
        endIndexNum.append(i)
    else:
        continue

#double check length of both lists
if len(modelIndexNum) != len(endIndexNum):
    print 'FILE ERROR, PLEASE CHECK THE KEYWORD MODEL AND ENDMDL IN FILE'
    exit()
else:
    for idx in range(len(modelIndexNum)):
        fileToWrite = open(name+'_'+str(idx+1)+'.pdb','w')
        fileToWrite.writelines(lines[modelIndexNum[idx]:endIndexNum[idx]+1])
        fileToWrite.close()



