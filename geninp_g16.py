#!/usr/bin/env python
#-*- coding: utf-8 -*-
# this script is to convert pre-processed pdb file into g16 input .com files:
# 1) for density functional theory calculation: b3lyp/6-31g
# 2) for semi-empirical calculation: zindo=(nstates=1)
# syntax: geninp_g16.py -p inputfile -n nprocshared -m memory -o outputfile

import numpy as np
from sys import argv, exit
import argparse

#script name
script = argv[0]

#parse input
#for help info, do geninp_g16.py -h (--help)
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest = 'input', help = 'Pre-processed PDB file for only QM part contained', type = str)
parser.add_argument('-n', dest = 'nprocshared', help = 'Number of CPU cores in the calculation', type = str)
parser.add_argument('-m', dest = 'memory', help = 'Memory required for the calculation, unit in GB', type = str)
parser.add_argument('-o', dest = 'output', help = 'Output name, input name as default', type = str)
#add further specific commend for analyzing
parser.add_argument('-s', dest = 'specific', nargs = '+', help = 'Specific feature for parsing input data', type = str)

######################################################
###############                        ###############
###############   PARAMETER SETTINGS   ###############
###############                        ###############
######################################################

# cainclude: not use H to replace CA, use H to replace C and N bonded with CA instead
# addh: add hydrogen to neutralize HEME carboxylic acid group
# nocap: just use the default coordinates and do no changes to it (warning: this will discard all other features and 
#        just do nocap) 

#parsing input
options = vars(parser.parse_args())

#get specific feature
specific_feature = map(lambda x:x.lower(),options['specific'])
#options validation
# -p argument
if options['input'] == None:
    print 'No input detected, please try again'
    exit()
else:
    inpfile = options['input']

# -n argument
if options['nprocshared'] == None:
    nprocshared = '4'
else:
    nprocshared = options['nprocshared']

# -m argument
if options['memory'] == None:
    mem = '4GB'
else:
    if 'GB' in options['memory'].upper():
        mem = options['memory'].upper()
    else:
        mem = options['memory']+'GB'

# -o argument
if options['output'] == None:
    name = inpfile.split('.')[0]
else:
    name = options['output']

#read pdb file
fin = open(inpfile, 'r')
lines = fin.readlines()
fin.close()

#write output .com file
fo_dft = open(name+'_dft.com', 'w')
fo_zindo = open(name+'_zindo.com', 'w')

#start reading PDB file
#should know the format of PDB inadvance
#PDB ATOM and HETATM format
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

######################################################
###############                        ###############
###############     INPUT PARSING      ###############
###############                        ###############
######################################################

#define coord list to store all new lines for input in g16 .com file
coord = []
#we mannually add hydrogen to replace the alpha-C as the default feature, with specific C-H distance
#pre-defined C-H bond distance = 1.09 A
C_H_bond= 1.09
#we add hydrogen to neutralize COO-, with specific O-H distance
O_H_bond= 0.98

#case 1: replace CA and not addh
if 'cainclude' not in specific_feature and 'addh' not in specific_feature:
    #since the output pdb file 1st line and -1st line are not coord, discard them in the for loop
    for i in range(1,len(lines)-1):
        line = lines[i] 
        if 'CA' == line[12:16].strip():
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
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        else:
            element = line[76:78].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())

            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)

#case 2: replace CA and addh
if 'cainclude' not in specific_feature and 'addh' in specific_feature:
    #since the output pdb file 1st line and -1st line are not coord, discard them in the for loop
    for i in range(1,len(lines)-1):
        line = lines[i] 
        if 'CA' == line[12:16].strip():
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
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        elif 'O1A' == line[12:16].strip() or 'O1D' == line[12:16].strip():
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
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)
            addline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(addline)
        else:
            element = line[76:78].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            #write line
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)

#case 3: replace N and C and not addh
if 'cainclude' in specific_feature and 'addh' not in specific_feature:
    #since the output pdb file 1st line and -1st line are not coord, discard them in the for loop
    for i in range(1,len(lines)-1):
        line = lines[i] 
        if 'N' == line[12:16].strip():
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
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        elif 'C' == line[12:16].strip():
            x_c = float(line[30:38].strip())
            y_c = float(line[38:46].strip())
            z_c = float(line[46:54].strip())
            #calculate the new coord
            length = ((x_c-x_ref)**2+(y_c-y_ref)**2+(z_c-z_ref)**2)**0.5
            newx = x_ref + (x_c - x_ref)*C_H_bond/length
            newy = y_ref + (y_c - y_ref)*C_H_bond/length
            newz = z_ref + (z_c - z_ref)*C_H_bond/length
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        else:
            element = line[76:78].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            #write line
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)

#case 4: replace N and C and addh
if 'cainclude' in specific_feature and 'addh' in specific_feature:
    #since the output pdb file 1st line and -1st line are not coord, discard them in the for loop
    for i in range(1,len(lines)-1):
        line = lines[i] 
        if 'N' == line[12:16].strip():
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
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        elif 'C' == line[12:16].strip():
            x_c = float(line[30:38].strip())
            y_c = float(line[38:46].strip())
            z_c = float(line[46:54].strip())
            #calculate the new coord
            length = ((x_c-x_ref)**2+(y_c-y_ref)**2+(z_c-z_ref)**2)**0.5
            newx = x_ref + (x_c - x_ref)*C_H_bond/length
            newy = y_ref + (y_c - y_ref)*C_H_bond/length
            newz = z_ref + (z_c - z_ref)*C_H_bond/length
            #get the new line added
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(newline)
        elif 'O1A' == line[12:16].strip() or 'O1D' == line[12:16].strip():
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
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)
            addline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format('H',newx,newy,newz)
            coord.append(addline)
        else:
            element = line[76:78].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            #write line
            newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
            coord.append(newline)

#case 5: no cap
if 'nocap' in specific_feature:
    #since the output pdb file 1st line and -1st line are not coord, discard them in the for loop
    for i in range(1,len(lines)-1):
        line = lines[i]
        #parse one line
        element = line[76:78].strip()
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        #write line
        newline = '{:>2}             {:>14.8f}{:>14.8f}{:>14.8f}\n'.format(element,x,y,z)
        coord.append(newline)


######################################################
###############                        ###############
###############    OUTPUT STREAMING    ###############
###############                        ###############
######################################################

#right now we just use it for general single point calculation, not for further parameter setting
#start writing g16 input files: dft version
#write header
fo_dft.write('%chk='+name+'_dft.chk\n')
fo_dft.write('%rwf='+name+'_dft.rwf\n')
fo_dft.write('%nprocshared='+nprocshared+'\n')
fo_dft.write('%mem='+mem+'\n')
#what to do in the calculation, default is single point 
#inplement in this #p command will be develped further
fo_dft.write('#p b3lyp/cc-pVDZ Pop=Full\n')
fo_dft.write('\n')
fo_dft.write('dft calculation for '+name+'\n')
fo_dft.write('\n')
#default charge and multiplicity is 0 1
fo_dft.write('0 1\n')
for line in coord:
    fo_dft.write(line)

#extend several blank lines
fo_dft.write('\n')
fo_dft.write('\n')

fo_dft.close()


#start writing g16 input files: zindo version
#write header
fo_zindo.write('%chk='+name+'_zindo.chk\n')
fo_zindo.write('%rwf='+name+'_zindo.rwf\n')
fo_zindo.write('%nprocshared='+nprocshared+'\n')
fo_zindo.write('%mem='+mem+'\n')
#what to do in the calculation, default is single point 
#inplement in this #p command will be develped further
fo_zindo.write('#p zindo=(nstates=1) Pop=Full\n')
fo_zindo.write('\n')
fo_zindo.write('zindo calculation for '+name+'\n')
fo_zindo.write('\n')
#default charge and multiplicity is 0 1
fo_zindo.write('0 1\n')
for line in coord:
    fo_zindo.write(line)

#extend several blank lines
fo_zindo.write('\n')
fo_zindo.write('\n')

fo_zindo.close()


