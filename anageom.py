#!/usr/bin/env python
#-*- coding: utf-8 -*-

#this script is to analysis geomatry/coordinates of specific molecule or system

from sys import argv, exit
import numpy as np



def centeroid(coordinate,atlist=None):
    #coordinate is list of lists of coordinates in x,y,z format
    #atlist is list of atom selection for calculating centeroid

    #since number in atom list is number, not list index, so we should -1 when use
    X = 0.0
    Y = 0.0
    Z = 0.0
    if atlist == None:
        atlist = range(1,len(coordinate)+1)

    for i in atlist:
        X += coordinate[i-1][0]
        Y += coordinate[i-1][1]
        Z += coordinate[i-1][2]
    #calculating centeroid
    alen = len(atlist)
    x = X/alen
    y = Y/alen
    z = Z/alen

    return [x, y, z]

def radius(coordinate,atlist=None):
    #coordinate is list of lists of coordinates in x,y,z format
    #atlist is list of atom selection for calculating radius

    #since number in atom list is number, not list index, -1 when in use
    if atlist == None:
        atlist = range(1,len(coordinate)+1)

    for i in atlist:
    X = [coordinate[i-1][0] for i in atlist]
    Y = [coordinate[i-1][1] for i in atlist]
    Z = [coordinate[i-1][2] for i in atlist]

    xspan = np.abs(X.min()-X.max())
    yspan = np.abs(Y.min()-Y.max())
    zspan = np.abs(Z.min()-Z.max())
    #arithmetic mean for radius
    radius_arithmetic = (xspan+yspan+zspan)/3
    #geometric mean for radisu
    xyz = xspan*yspan*zspan
    if xyz > 0.0:
        radius_geometric = xyz**(1./3)
    else:
        radius_geometric = -(-xyz)**(1./3)

    return [raidus_arithmetic, radius_geometric]

def distance(a,b):
    #calculate the distance between two points
    #both a and b are list of coordinates
    return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5 

def edge2edgedis(alist,blist):
    #calculate shortest distance between two lists of coordinates a and b
    pairs = [(a,b) for a in alist for b in blist]
    return np.argmin([distance(a,b) for (a,b) in pairs])


