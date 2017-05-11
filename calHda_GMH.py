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

