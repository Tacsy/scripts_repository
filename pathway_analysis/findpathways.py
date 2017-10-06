#!/usr/bin/env python
#This python script is to extract the pathway in a set of raw
#pathway files and summarize in a single file

#syntax: ./findpathways.py filenum #ofpathways
from sys import argv

#take the first argument as the number of top-n pathways
script, filenum, num = argv

#outputfile initialize
pathwayout = open('sum_of_pathway_'+filenum+'.log','w')

#read inputfile one by one
for i in range(int(num)):
    #open the file
    opener = open('pathway_'+filenum+'_'+str(i+1)+'.log','r')
    #read and store the lines
    path = opener.readlines()
    #close the file
    opener.close()

    #put the first line of pathway and relative coupling in the output file
    line = '#'+path[0]
    pathwayout.write(line)

    #analyze the pathway and then write it in the output file
    pathway = []
    for idx in range(1,len(path)):
        keys = path[idx].split()
        #get the residue, with resname+resnum
        res = keys[2]+keys[3]
        if res not in pathway:
            pathway.append(res)
        else:
            continue

    #to make sure acceptor is not the only atom of that residue associated with
    #the pathway so that it can not be displayed in this script, double check the
    #last atom
    lastkeys = path[-1].split()
    lastres = lastkeys[8]+lastkeys[9]
    if lastres not in pathway:
        pathway.append(lastres)
    
    pathwayout.write('>'.join(pathway)+'\n')


#after analyzing all the files, close the output file
pathwayout.close()

