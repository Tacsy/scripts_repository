#!/usr/bin/env python
#-*- coding: utf-8 -*-

# this python script is to calculate outer-sphere reorganization energy between two sites with Marcus two-sphere model.

# syntax: ./calReorg_marcus.py -es epsilon static index value defualt is 4 in general protein environment
# then provide D diameter, A diameter, DA distance in Angstrom
# when want exit, type "exit"

from sys import argv, exit
import numpy as np
import argparse

#script name
script = argv[0]

#parser setup
#for help information, use ./calReorg_marcus.py -h [--help]
parser = argparse.ArgumentParser()
parser.add_argument('-es', dest = 'epsilonS', help = 'static dielectric constants, if not defined, use 4 as the typical value in protein environment', type = float)

#parse the arguments
options = vars(parser.parse_args())

if options['epsilonS'] == None:
    es = 4.0
else:
    es = options['epsilonS']

while True:
    inputstr = raw_input('type in D diameter, A diameter, DA distance in angstrom\n')
    if inputstr.lower() == 'exit':
        exit()
    else:
        if len(inputstr.split()) != 3:
            print "Wrong input format, please try again!"
            continue
        else:
            Ddiameter = float(inputstr.split()[0])
            Adiameter = float(inputstr.split()[1])
            DAdistance = float(inputstr.split()[2])
            
            epsilon = 8.854*10**(-12) #unit: F/m
            eo = 2.0 #optic dielectric constant
            E = 1.602*10**(-19) #unit: C

            reorg = E*(1.0/eo-1.0/es)*(1.0/(4*np.pi*epsilon))*(10**10)*(1.0/Ddiameter + 1.0/Adiameter - 1.0/DAdistance)

            print "The outer-sphere reorganization energy is:  "+str(reorg)+"  eV"

            
