#!/usr/bin/env python 
#-*- coding: utf-8 -*-

# this python script is to calculate electronic coupling between two specific frontier orbital using GMH method.

# steps
# 1) build matrix: Dipole matrix in XYZ, MO matrix
# 2) use GMH method to calculate the coupling
# 3) output streaming

# syntax: calHda_GMH.py -d dipomatrix -m MOcoefmatrix 


from sys import argv,exit
import numpy as np
import numpy.linalg as la
from scipy.linalg import sqrtm
from math import sqrt
import argparse

#script name
script = argv[0]

#parser setup
#for help information, use ./calHda_GMH.py -h [--help]
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest = 'dipole', help = 'Dipole matrix, XYZ lower triangular form', type=str)
parser.add_argument('-m', dest = 'MOcoef', help = 'MO coefficient', type=str)
parser.add_argument('-f', dest = 'fock', help = 'Fock matrix', type=str)
parser.add_argument('-s', dest = 'overlap', help = 'Overlap matrix', type=str)
parser.add_argument('-no', dest = 'orbital', nargs = '+', help = 'Number of donor and acceptor orbitals (not index)',  type=int)
parser.add_argument('-o', dest = 'outfile', help = 'Output prefix, gmh as default', type=str)

#parse the arguments and store them in the defined dictionary
options = vars(parser.parse_args())

#input validation
if options['dipole'] == None:
    print 'No dipole matrix file found, please try again'
    exit()
else:
    dipofile = options['dipole']
if options['MOcoef'] == None:
    print 'No MO coefficient matrix file found, please try again'
    exit()
else:
    mofile = options['MOcoef']
if options['fock'] == None:
    print 'No fock matrix file found, please try again'
    exit()
else:
    fockfile = options['fock']
if options['overlap'] == None:
    print 'No overlap matrix file found, please try again'
    exit()
else:
    overfile = options['overlap']
if options['outfile'] == None:
    outprefix = 'gmh'
else:
    outprefix = options['outfile']
if options['orbital'] == None:
    print 'No orbital specified, please try again'
    exit()
else:
    dOrb,aOrb = options['orbital']
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

########################################
def gmh(dipX,dipY,dipZ,MO,energy,dOrb,aOrb):

    #AO basis dipole to MO basis dipole
    MOdipX = MO * dipX * MO.transpose()
    MOdipY = MO * dipY * MO.transpose()
    MOdipZ = MO * dipZ * MO.transpose()

    #get the overall dipole
    #X = np.trace(MOdipX)
    #Y = np.trace(MOdipY)
    #Z = np.trace(MOdipZ)

    #since both donor and acceptor orbitals are their corresponding number, not index, minus 1 when in use.
    H11 = energy[dOrb-1]
    H22 = energy[aOrb-1]
    delE = np.abs(H11-H22)
    mu11 = np.array([MOdipX[dOrb-1,dOrb-1], MOdipY[dOrb-1,dOrb-1], MOdipZ[dOrb-1,dOrb-1]])
    mu22 = np.array([MOdipX[aOrb-1,aOrb-1], MOdipY[aOrb-1,aOrb-1], MOdipZ[aOrb-1,aOrb-1]])
    mu12 = np.array([MOdipX[dOrb-1,aOrb-1], MOdipY[dOrb-1,aOrb-1], MOdipZ[dOrb-1,aOrb-1]])
    delmu = mu11 - mu22
    Hda = (la.norm(mu12) * delE * 27.2114)/ np.sqrt(la.norm(delmu)**2 + 4.0 * la.norm(mu12)**2) 

    return Hda, X ,Y, Z
########################################

#main function
#build matrix
fock, fock_dim = build_matrix(fockfile)
over, over_dim = build_matrix(overfile)

if fock_dim != over_dim:
    print 'The dimension of both matrixs don\'t match. Please look into it and try again'
    exit()

U_eval, U_evec = la.eig(over)
over_prime = np.matrix(np.diag(U_eval))
over_T = U_evec * la.inv(sqrtm(over_prime)) * U_evec.transpose()
fock_O = over_T.transpose() * fock * over_T

MO_val, MO_vec = la.eig(fock_O)
MO_index = MO_val.argsort()
energy = MO_val[MO_index]

dipX, dipY, dipZ, dip_dim = build_dipo_matrix(dipofile)
MO = build_mo_matrix(mofile)

#use gmh mathod to calculate the coupling and overall dipole
Hda = gmh(dipX, dipY, dipZ, MO, energy, dOrb, aOrb)

#################################################
###############                   ###############  
###############  OUTPUT STREAMING ############### 
###############                   ###############
#################################################

#gmh derived coupling and overall dipole output
gmh_out = open(outprefix+'_coupling.dat', 'w')

gmh_out.write('%e\n'%Hda)
#gmh_out.write('%e\n'%X)
#gmh_out.write('%e\n'%Y)
#gmh_out.write('%e\n'%Z)

gmh_out.close()
