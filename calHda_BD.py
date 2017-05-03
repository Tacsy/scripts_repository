#!/usr/bin/env python
#-*- coding: utf-8 -*-

# this python script is to calculate electronic coupling between two specific frontier orbital using block diagonalization method.

# steps
# 1) build whole matrix from the raw lower triangle matrix
# 2) block-diagonalize the matrix 
# 3) get the wanted orbital and the related site energy and coupling

# syntax: ./calHda_BD.py -f fockmatrix -s overlapmatrix -no number_of_orbital -ne number_of_electron -o output

from sys import argv,exit
import numpy as np
import numpy.linalg as la
from scipy.linalg import block_diag,sqrtm
from math import sqrt
import argparse

#script name
script = argv[0]

#parser setup
#for help information, use ./calHda_BD.py -h [--help]
parser = argparse.ArgumentParser()
parser.add_argument('-f', dest = 'fock', help = 'Fock matrix, lower triangle form', type = str)
parser.add_argument('-s', dest = 'overlap', help = 'Overlap matrix, lower triangle form', type = str)
parser.add_argument('-nb', dest = 'numofbasis', nargs = '+', help = 'Number of basis in each block', type = int)
parser.add_argument('-ne', dest = 'numofelectron', nargs = '+', help = 'Number of electron in each block', type = int)
parser.add_argument('-o', dest = 'output', help = 'Output file prefix, block_diag as the default', type = str)

#parse the arguments and store them in the defined directory
options = vars(parser.parse_args())

#test
#print options

#input validation
if options['output'] == None:
    out_prefix = 'block_diag'
else:
    out_prefix = options['output']

basislist = options['numofbasis']
electronlist = options['numofelectron']

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
    tria_num = len(fock_list)
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
            full_matrix[i.j] = np.conjugate(full_matrix[j,i])

    return full_matrix, dim
#######################################

#get the full matrix and examine the dimension
fockfile = options['fock'] 
overfile = options['overlap']

fock_full, fock_dim = build_matrix(fockfile)
over_full, over_dim = build_matrix(overfile)

if fock_dim != over_dim:
    print "The dimension of both files don't match. Please look into it and try again."
    exit()

#########################
# 2) block-diagonlize 
#########################

# need to know 
# basis set of this calculation is required for calculating the electronic coupling correctly

#############################################
def build_BD(fock,over,basis):
    #get the dimension of the matrix
    dim = np.shape(fock)[0]
    
    #test the basis input
    if sum(basis) != dim:
        print "Input in -nb section is wrong! Please look into it and try again."
        exit()

    #get np.matrix format
    fock = np.matrix(fock)
    over = np.matrix(over)

    #transform the fock matrix in orthogonal basis
    U_eval, U_evec = la.eig(over)
    over_prime = np.matrix(np.diag(U_eval))
    over_T = U_evec * la.inv(sqrtm(over_prime)) * U_evec.transpose()
    fock_O = over_T.transpose() * fock * over_T

    #block diagonalize submatries
    submatrix = []
    num_basis = 0
    #store each submatrix
    for i in range(len(basis)):
        submatrix.append(fock_O[num_basis:num_basis+basis[i],num_basis:num_basis+basis[i]])
        num_basis = num_basis + basis[i]
    
    #build both block eval and evec
    eval_sub = []
    evec_sub = []
    for matrix in submatrix:
        eval_tmp, evec_tmp = la.eig(matrix)
        index = eval_tmp.argsort()
        eval_sub.append(eval_tmp[index])
        evec_sub.append(evec_tmp[:,index])

    #get the block diagonalized fock matrix
    U = block_diag(*evec_sub)
    F = U.transpose() * fock_O * U
    
    #for test, get the full MO Fock matrix
    MO_val, MO_vec = la.eig(fock_O)
    MO_index = MO_val.argsort()
    MO_val = MO_val[MO_index]

    return F, MO_val
#############################################

#########################
# 3) get reduced fock  
#########################

#############################################
def get_BD(fock,basis,electron):
    #since it's a closed shell calculation, both alpha and beta electron stay in the same orbital
    #so the orbital we want for HOMO and hole transfer is divided by 2 from the total electron
    dim = len(electron)

    orb_index = []
    num_basis = 0
    for i in range(len(electron)):
        orb = (electron[i]+1)/2
        orb_index.append(num_basis+orb-1)
        num_basis = num_basis + basis[i]

    reduced_F = []
    for i in orb_index:
        for j in orb_index:
            reduced_F.append(fock[i,j])

    F = np.array(reduced_F).reshape([dim,dim])

    return F
#############################################

dim = fock_dim

fock_BD, MO_val = build_BD(fock_full,over_full,basislist)

reduced_F = get_BD(fock_BD,basislist,electronlist)

#################################################
###############                   ###############  
###############  OUTPUT STREAMING ############### 
###############                   ###############
#################################################

#block diagnolized fock matrix output
fock_BD_out = open(out_prefix+'_BDmatrix.dat', 'w')

for i in range(dim):
    for j in range(dim):
        fock_BD_out.write('%e    '%fock_BD[i,j].real)
    fock_BD_out.write('\n')

fock_BD_out.close()

#reduced-fock matrix output
fock_reduce_out = open(out_prefix+'_redF.dat', 'w')

unit = len(basislist)
for i in range(unit):
    for j in range(unit):
        fock_reduce_out.write('%e     '%reduced_F[i,j].real)
    fock_reduce_out.write('\n')

fock_reduce_out.close()

#full MO energy output
MO_out = open(out_prefix+'_MO.dat', 'w')

for i in range(len(MO_val)):
    MO_out.write('%e'%MO_val[i].real)
    MO_out.write('\n')

MO_out.close()


