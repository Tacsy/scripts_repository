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

 
