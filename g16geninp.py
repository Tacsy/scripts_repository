#!/usr/bin/env python
#-*- coding: utf-8 -*-
# this script is to convert pre-processed pdb file into g16 template input .com file:
# syntax: g16geninp.py -p input -n nproc -m mem -o output [specific function]

import numpy as np
from sys import argv, exit
import argparse

#script name
script = argv[0]

#parse input
#for help info, use g16geninp.py -h (--help)
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest = 'input', help = 'Pre-processed PDB file with only QM part contained', type = str)
parser.add_argument('-n', dest = 'nproc', help = 'Number of CPU cores', type = str)
parser.add_argument('-m', dest = 'mem', help = 'Memory required for the calculation, unit in GB', type = str)
parser.add_argument('-o', dest = 'output', help = 'Output name, input name as default', type = str)
#add further specific commend for analyzing
parser.add_argument('-s', dest = 'specific', nargs = '+', help = 'Specific feature for parsing input data', type = str)
parser.add_argument('-ct', dest = 'Cterminal', help = 'C terminal residue number', type = str)
parser.add_argument('-nt', dest = 'Nterminal', help = 'N terminal residue number', type = str)

######################################################
###############                        ###############
###############   PARAMETER SETTINGS   ###############
###############                        ###############
######################################################

# cainclude: not use H to replace CA, use H to replace C and N bonded with CA instead
# addh: add hydrogen to neutralize HEME carboxylic acid group
# nocap: just use the default coordinates and do no changes to it (warning: this will discard all other features and 
#        just do nocap) 
# capback: cap backbone with hydrogen to NH and OH to C=O

#parsing input
options = vars(parser.parse_args())

#options validation
# -p argument
if options['input'] == None:
    print 'No input detected, please try again'
    exit()
else:
    inpfile = options['input']

# -n argument
if options['nproc'] == None:
    nprocshared = '4'
else:
    nprocshared = options['nprocshared']

# -m argument
if options['mem'] == None:
    mem = '4GB'
else:
    if 'GB' in options['mem'].upper():
        mem = options['mem'].upper()
    else:
        mem = options['mem']+'GB'

# -o argument
if options['output'] == None:
    name = inpfile.split('.')[0]
else:
    name = options['output']

# -s argument
if options['specific'] == None:
    specific_feature = []
else:
    specific_feature = map(lambda x:x.lower(),options['specific'])

#read pdb file
fin = open(inpfile, 'r')
rawlines = fin.readlines()
fin.close()

#remove the redundent lines (like ter or comment)
lines = []
for line in rawlines:
    if 'ATOM' in line or 'HETATM' in line:
        lines.append(line)
    else:
        continue

#start parsing PDB file
#should know PDB format in advance
#PDB ATOM and HETATM
#c 1-6 ATOM or HETATM
#c 7-11 serial number
#c 13-16 atom name
#c 18-20 residue name
#c 22  chainID
#c 23-26 residue sequence number
#c 31-38 x coordinate
#c 39-46 y coordinate
#c 47-54 z coordinate
#c 77-78 element

#######################################################
###############                         ###############
###############   FUNCTION DEFINITION   ###############
###############                         ###############
#######################################################

#convert function
#from pdb format to com format coordinates
def pdb2com(lines):
    newlines = []
    for i in range(len(lines)):
        line = lines[i]
        #parse one line
        element = line[76:78].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        #write line
        newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
        newlines.append(newline)

    return newlines

#refine function
#from pdb format to pdb format, only refine several lines
def cap_backbone(lines, nterminal, cterminal):
    newlines = []
    for i in range(len(lines)):
        line = lines[i]
        #parse one line
        if line[22:26].strip() == nt:
            if line[12:16].strip() == 'C':
                #write line
                linelist = list(line)
                #revise atom type and atom name
                linelist[12:16] = ' HR '
                linelist[76:78] = ' H'
                newline = ''.join(linelist)
                newlines.append(newline)
                continue
            else:
                newlines.append(line)
                continue
        elif line[22:26].strip() == ct:
            if line[12:16].strip() == 'N':
                #write line
                linelist = list(line)
                #revise atom type and atom name
                linelist[12:16] = ' OR '
                linelist[76:78] = ' O'
                newline = ''.join(linelist)
                newlines.append(newline)
                continue
            elif line[12:16].strip() == 'HN':
                #write line
                linelist = list(line)
                #revise atom type and atom name
                linelist[12:16] = ' HR '
                linelist[76:78] = ' H'
                newline = ''.join(linelist)
                newlines.append(newline)
                continue
            else:    
                newlines.append(line)
                continue
        else:
            newlines.append(line)
    
    return newlines

def addh_heme(lines):
    #default value for OH bond length
    O_H_bond= 0.98
    newlines = []
    for i in range(len(lines)):
        line = lines[i]
        #parse one line
        if (line[12:16].strip() == 'O1A' and line[17:20].strip() == 'HEME') or (line[12:16].strip() == 'O1D' and line[17:20].strip() == 'HEME'):
            element = line[76:78].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            #get coord of CG
            refline = lines[i-2]
            x_ref = float(refline[30:38].strip())
            y_ref = float(refline[38:46].strip())
            z_ref = float(refline[46:54].strip())
            #calculate the new coord
            CO_length = ((x-x_ref)**2+(y-y_ref)**2+(z-z_ref)**2)**0.5
            newx = x + (x - x_ref)*O_H_bond/CO_length
            newy = y + (y - y_ref)*O_H_bond/CO_length
            newz = z + (z - z_ref)*O_H_bond/CO_length
            #write two line 
            addline = line
            linelist = list(addline)
            #revise atom type and atom name
            if line[12:16].strip() == 'O1A':
                linelist[12:16] = ' HR1'
            else:
                linelist[12:16] = ' HR2'
            linelist[30:54] = '{:>8.3f}{:>8.3f}{:>8.3f}'.format(newx,newy,newz)
            #linelist[38:46] = '{:>8.3f}'.format(newy)
            #linelist[46:54] = '{:>8.3f}'.format(newz)
            linelist[76:78] = ' H'
            newline = ''.join(linelist)
            newlines.append(line)
            newlines.append(newline)
        else:
            newlines.append(line)

    return newlines

def cap_sidechain(lines):
    C_H_bond = 1.09
    newlines = []
    for i in range(len(lines)):
        line = lines[i]
        if line[12:16].strip() == 'CA':
            x_ca = float(line[30:38].strip())
            y_ca = float(line[38:46].strip())
            z_ca = float(line[46:54].strip())
            #get coord of beta-C
            refline = lines[i+1]
            x_ref = float(refline[30:38].strip())
            y_ref = float(refline[38:46].strip())
            z_ref = float(refline[46:54].strip())
            #calculate the new coord
            length = ((x_ca-x_ref)**2+(y_ca-y_ref)**2+(z_ca-z_ref)**2)**0.5
            newx = x_ref + (x_ca - x_ref)*C_H_bond/length
            newy = y_ref + (y_ca - y_ref)*C_H_bond/length
            newz = z_ref + (z_ca - z_ref)*C_H_bond/length
            #get the new coordinates replaced 
            linelist = list(line)
            linelist[12:16] = ' HR '
            linelist[30:54] = '{:>8.3f}{:>8.3f}{:>8.3f}'.format(newx,newy,newz)
            #linelist[38:46] = '{:>8.3f}'.format(newy)
            #linelist[46:54] = '{:>8.3f}'.format(newz)
            linelist[76:78] = ' H'
            newline = ''.join(linelist)
            newlines.append(newline)
        else:
            newlines.append(line)
    return newlines

def cap_sidechain_more(lines):
    C_H_bond = 1.09
    newlines = []
    for i in range(len(lines)):
        line = lines[i]
        if line[12:16].strip() == 'N':
            x_n = float(line[30:38].strip())
            y_n = float(line[38:46].strip())
            z_n = float(line[46:54].strip())
            #get coord of alpha-C
            refline = lines[i+1]
            x_ref = float(refline[30:38].strip())
            y_ref = float(refline[38:46].strip())
            z_ref = float(refline[46:54].strip())
            #calculate the new coord
            length = ((x_n-x_ref)**2+(y_n-y_ref)**2+(z_n-z_ref)**2)**0.5
            newx = x_ref + (x_n - x_ref)*C_H_bond/length
            newy = y_ref + (y_n - y_ref)*C_H_bond/length
            newz = z_ref + (z_n - z_ref)*C_H_bond/length
            #get the new coordinates replaced
            linelist = list(line)
            linelist[12:16] = ' HR1'
            linelist[30:54] = '{:>8.3f}{:>8.3f}{:>8.3f}'.format(newx,newy,newz)
            #linelist[38:46] = '{:>8.3f}'.format(newy)
            #linelist[46:54] = '{:>8.3f}'.format(newz)
            linelist[76:78] = ' H'
            newline = ''.join(linelist)
            newlines.append(newline)
        elif line[12:16].strip() == 'C':
            x_c = float(line[30:38].strip())
            y_c = float(line[38:46].strip())
            z_c = float(line[46:54].strip())
            #calculate the new coord
            length = ((x_c-x_ref)**2+(y_c-y_ref)**2+(z_c-z_ref)**2)**0.5
            newx = x_ref + (x_c - x_ref)*C_H_bond/length
            newy = y_ref + (y_c - y_ref)*C_H_bond/length
            newz = z_ref + (z_c - z_ref)*C_H_bond/length
            #get the new coordinates replaced
            linelist = list(line)
            linelist[12:16] = ' HR2'
            linelist[30:54] = '{:>8.3f}{:>8.3f}{:>8.3f}'.format(newx,newy,newz)
            #linelist[38:46] = '{:>8.3f}'.format(newy)
            #linelist[46:54] = '{:>8.3f}'.format(newz)
            linelist[76:78] = ' H'
            newline = ''.join(linelist)
            newlines.append(newline)
        else:
            newlines.append(line)
    return newlines

######################################################
###############                        ###############
###############    OUTPUT STREAMING    ###############
###############                        ###############
######################################################

#start parsing and refining the coordinates

#find specific feature in -s argument
if 'addh' in specific_feature:
    lines = addh_heme(lines)
if 'cap_back' in specific_feature:
    if options['Cterminal'] == None or options['Nterminal'] == None:
        print 'No specific terminal residue number decleared, please try again'
        exit()
    else:
        ct = options['Cterminal']
        nt = options['Nterminal']
        lines = cap_backbone(lines, nt, ct)
if 'cap_side' in specific_feature:
    lines = cap_sidechain(lines)
if 'cap_sidemore' in specific_feature:
    lines = cap_sidechain_more(lines)

#after these specific features, convert pdb to com format
coordinates = pdb2com(lines)


#start writing g16 input files
#write header
fo = open(name+'_template.com', 'w')
fo.write('%chk='+name+'_template.chk\n')
fo.write('%rwf='+name+'_template.rwf\n')
fo.write('%nprocshared='+nprocshared+'\n')
fo.write('%mem='+mem+'\n')
#inplement in this #p command will be develped further
fo.write('#p functional/basis Pop=Full keyword\n')
fo.write('\n')
fo.write('ab initio calculation for '+name+'\n')
fo.write('\n')
#default charge and multiplicity is 0 1
fo.write('0 1\n')
for line in coordinates:
    fo.write(line)

#extend several blank lines
fo.write('\n')

fo.close()


