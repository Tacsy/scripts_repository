#!/usr/bin/env python
#-*- coding: utf-8 -*-

#this script is to convert GROMACS topology file (.itp) to psf file for every single monomer.

#syntax: ./top2psf.py -p input -o output
#caution: this script can be processed one at a time, then use catpsf.py to concatenate them together.


from sys import argv, exit
import argparse

script = argv[0]

#parse input
#for help info, do top2psf.py -h (--help)
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest = 'top', help = 'GROMACS topology .itp file containing needed info', type = str)
parser.add_argument('-o', dest = 'output', help = 'Output name, defualt is $filename.psf', type = str)

#parsing input
options = vars(parser.parse_args())

if options['top'] == None:
    print 'No input detected, please look into it and try again'
    exit()
else:
    fullname = options['top']
    filename = fullname.split('.')[0]
#options validation
if options['output'] == None:
    outfile = filename+'.psf'
else:
    outfile = options['output']

#readline for parsing
f = open(fullname,'r')
lines = f.readlines()
f.close()

######################################################
###############                        ###############
###############     INPUT PARSING      ###############
###############                        ###############
######################################################

#define container list for storing relevant info
ATOMS = []
BONDS = []
bond_tmp = []

#get the index of specific section, here we only need atoms and bonds info
atom = lines.index('[ atoms ]')
bond = lines.index('[ bonds ]')
pair = lines.index('[ pairs ]')

#read atom info
for i in range(atom+1,bond):
    
    line = lines[i]
    if line[0] == ';' or line == '\n':
        continue
    
    #do real line parsing
    element = line.split()
    num = element[0]
    atomtype = element[1]
    resnum = element[2]
    resname = element[3]
    atom = element[4]
    charge = float(element[6])
    mass = float(element[7])
    
    #line formatting
    newline = '{:>8s}{:>5s} {:<5s}{:<5s}{:<5s}{:<5s}{:>14.10f} {:>14.10f}\n'.format(num,'MAIN',resnum,resname,atom,atomtype,charge,mass)

    ATOMS.append(newline)

#read bond info
for i in range(bond+1,pair):

    line = lines[i]
    if line[0] == ';' or line == '\n':
        continue

    #do real line parsing
    element = line.split()
    bond_tmp.append(element[0])
    bond_tmp.append(element[1])

    if len(bond_tmp) == 8:
        newline = '{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}\n'.format(*bond_tmp)
        BONDS.append(newline)
        bond_tmp[:] = []



######################################################
###############                        ###############
###############    OUTPUT STREAMING    ###############
###############                        ###############
######################################################

out = open(outfile, 'w')
#header streaming
out.write('PSF\n')
out.write('       4 !NTITLE\n')
out.write(' REMARKS TOPOLOGY CREATED BY PYTHON SCRIPT (XUYAN RU)')
out.write(' REMARKS CONTACT XUYAN.RU@DUKE.EDU')
out.write('\n')

#atom info streaming
num_of_atom = len(ATOMS)
out.write('{:>8s} !NATOM\n'.format(str(num_of_atom)))

for line in ATOMS:
    out.write(line)

out.write('\n')

#bond info streaming
num_of_bond = len(BONDS)*4 + len(bond_tmp)/2
out.write('{:>8s} !NBOND: bonds\n'.format(str(num_of_bond)))

for line in BONDS:
    out.write(line)

lastline = ''
for string in bond_tmp:
    lastline = lastline + '{:>8s}'.format(string)

lastline = lastline + '\n'
 
out.write(lastline)
out.write('\n') 
#finish writing the output file, close the file
out.close()

