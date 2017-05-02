#!/usr/bin/env python
#-*- coding: utf-8 -*-

#this script is to convert GROMACS topology file (.itp) to psf file for every single monomer.

#syntax: ./top2psf.py -p input -o output
#caution: this script can be processed one at a time, then use catpsf.py to concatenate them together.


from sys import argv, exit
import argparse

script = argv[0]


