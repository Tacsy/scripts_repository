#!/usr/bin/env python
#This python script is to analyze a bunch of pathways files and summerize them
#in a new output file

#syntax: ./analyzepathways.py number_of_frames

from sys import argv

script, num = argv

#we define a dictionary to store all the pathway informations
#specifically, every key represents a pathway, and every value
#represents the number of this pathway been seen in all defined
#pathways.

#define pathway dictionary
pathway_dict = {}

for i in range(int(num)):
    
    #open file
    opener = open('sum_of_pathway_'+str(i+1)+'.log','r')
    #read all lines and store in a list
    path = opener.readlines()
    #close the file
    opener.close()

    for line in path:
        #discard the comment line
        if '#' in line:
            continue
        else:
            #define the key
            key = line[:-1]
            if key not in pathway_dict.keys():
                #initialize the key and value
                pathway_dict[key] = 1
            else:
                pathway_dict[key] += 1
    

#write output file
outputfile = open('SUMMERY_PATHWAY.log','w')

#write lines, every line with key:value pair
for key in pathway_dict.keys():
    #generate the lines
    line = key+':'+str(pathway_dict[key])+'\n'
    
    outputfile.write(line)

#after summerizing all the pathways, close the outputfile
outputfile.close()
