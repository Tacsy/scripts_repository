#!/usr/bin/env python 
#-*- coding: utf-8 -*-

# this python script is to calculate electronic coupling between two specific frontier orbital using GMH method.

# steps
#
#
#
#

# syntax:


from sys import argv,exit
import numpy as np
import numpy.linalg as la
from math import sqrt
import argparse

#script name
script = argv[0]

#parser setup
#for help information, use ./calHda_GMH.py -h [--help]
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest = 'dipole', help = 'Dipole matrix, XYZ lower triangular form', type=str)
parser.add_argument('-m', dest = 'MOcoef', help = 'MO coefficient', type=str)
parser.add_argument()
parser.add_argument()

#parse the arguments and store them in the defined dictionary
options = vars(parser.parse_args())

#########################
# 1) build whole matrix
#########################

#######################################
def build_matrix(filename):
    #open file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    #read all lower triangle matrix elements in the list
    numlist = []
    for line in lines:
        linelist = line.split()
        for num in linelist:
            numlist.append(float(num))

    #get the length of the list and get the dimension
    tria_num = len(numlist)
    for i in range(int(sqrt(tria_num*2)),1,-1):
        if i*(i+1)/2.0 == tria_num:
            dim = i
            break
        else:
            print "The input parsed is not a lower triangle matrix, please look into it and try again"
            exit()
    
    #build full matrix
    full_matrix = np.matrix(np.zeros((dim,dim)))
    #set the original point (0,0) as the first element
    full_matrix[0,0] = numlist[0]
    idx = 1

    #build lower triangular matrix element
    for i in range(1,dim):
        for j in range(0,i+1):
            full_matrix[i,j] = numlist[idx]
            idx = idx + 1

    #build upper triangular matrix element
    for i in range(dim):
        for j in range(i+1,dim):
            full_matrix[i,j] = np.conjugate(full_matrix[j,i])

    return full_matrix, dim
########################################

########################################
def build_dipo_matrix(filename):
   
    #open file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    #read all lower triangle matrix elements in the list
    numlist = []
    for line in lines:
        linelist = line.split()
        for num in linelist:
            numlist.append(float(num))

    #get the length of the list and get the dimension
    tria_num = len(numlist)
    
    #test wheather it can be divided by 3
    if tria_num % 3 == 0:
        num = tria_num/3
    else:
        print "The input parsed is not a 3D dipole matrix, please look into it and try again"
        exit()
    
    for i in range(int(sqrt(num*2)),1,-1):
        if i*(i+1)/2.0 == num:
            dim = i
            break
        else:
            print "The input parsed is not a lower triangle matrix, please look into it and try again"
            exit()
   
    #build X,Y,Z dipole matrix
    dipX = np.matrix(np.zeros((dim,dim)))
    dipY = np.matrix(np.zeros((dim,dim)))
    dipZ = np.matrix(np.zeros((dim,dim)))

    idx = 0

    #build lower triangular matrix element
    for i in range(dim):
        for j in range(i+1):
            dipX[i,j] = numlist[num*0+idx]
            dipY[i,j] = numlist[num*1+idx]
            dipZ[i,j] = numlist[num*2+idx]
            idx = idx + 1

    #build upper triangular matrix element
    for i in range(dim):
        for j in range(i+1,dim):
            dipX[i,j] = np.conjugate(dipX[j,i])
            dipY[i,j] = np.conjugate(dipY[j,i])
            dipZ[i,j] = np.conjugate(dipZ[j,i])

    return dipX, dipY, dipZ, dim
########################################

########################################
def build_mo_matrix(filename):
    
    #open file
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    #read all lower triangle matrix elements in the list
    numlist = []
    for line in lines:
        linelist = line.split()
        for num in linelist:
            numlist.append(float(num))

    #get the length of the list and get the dimension
    num = len(numlist)
    dim = int(sqrt(num))

    #build MO matrix
    MO = np.matrix(np.zeros((dim,dim)))

    idx = 0
    for i in range(dim):
        for j in range(dim):
            MO[i,j] = numlist[idx]
            idx = idx + 1

    return MO
########################################

