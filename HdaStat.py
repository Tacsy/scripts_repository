#!/usr/bin/env python
#-*- coding: utf-8 -*-
#this python script is to parse a file with only electronic coupling and calculate the following terms:
#1) mean of squared coupling <H^2>
#2) mean of coupling <H>
#3) variance of coupling <H^2>-<H>^2

#syntax: ./HdaStat.py filename

from sys import argv
import numpy as np

#get script and filename
script, filename = argv

fin = open(filename, 'r')
lines = fin.readlines()
fin.close()

length = len(lines)

Hda = []
Hda_rev = []
for line in lines:
    Hda.append(float(line.split()[0]))
    Hda_rev.append(float(line.split()[1]))

#doing statistics

mean_Hda = np.mean(Hda)
var_Hda = np.var(Hda)
msq_Hda = var_Hda + mean_Hda**2.0

mean_Hda_rev = np.mean(Hda_rev)
var_Hda_rev = np.var(Hda_rev)
msq_Hda_rev = var_Hda_rev + mean_Hda_rev**2.0

#doing test: calculating by hand
#test_mean_Hda = np.sum(Hda)/len(Hda)
#test_msq_Hda = np.sum(map(lambda x:x*x, Hda))/len(Hda)
#test_var_Hda = test_msq_Hda - test_mean_Hda**2.0

#output
print 'Statistical Analysis of Electronic Coupling (eV):'
print 'Highest coupling'
print 'The mean squared coupling is:  {}'.format(msq_Hda)
print 'The mean coupling is:  {}'.format(mean_Hda)
print 'The varience of coupling is:  {}'.format(var_Hda)

print ''
print 'Statistical averaged coupling'
print 'The mean squared coupling is:  {}'.format(msq_Hda_rev)
print 'The mean coupling is:  {}'.format(mean_Hda_rev)
print 'The varience of coupling is:  {}'.format(var_Hda_rev)
#test output
#print 'test output:'
#print 'The mean squared couplnig is: {}'.format(test_msq_Hda)
#print 'The mean coupling is: {}'.format(test_mean_Hda)
#print 'The varience of coupling is: {}'.format(test_var_Hda)

