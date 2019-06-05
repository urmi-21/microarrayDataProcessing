# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:47:27 2019

@author: mrbai
[1] input file (from MOG)
[2] output file

Example input file
Tx  Tx  Correlation
tx1 tx1 1.0000
tx1 tx2 0.1322
tx1 tx3 -0.0072
tx1 tx4 -0.2158
tx1 tx5 -0.1029
tx1 tx6 0.1522
tx1 tx7 0.1834
tx1 tx8 0.0444
tx1 tx9 -0.1598
tx1 tx10    -0.0731
"""


import markov_clustering as mc
import networkx as nx
import random
import sys
import numpy as np
import scipy.sparse as sparse

#read input file to create adj matrix

print("reading data")
with open(sys.argv[1],"r") as f:
        data=f.read().splitlines()
#remove header
data.pop(0)
#create a dict for gene name to id
nameDict={}
gid=0
gidToName={}
for l in data:
        temp=l.split("\t")
        if(not temp[0] in nameDict):
                nameDict[temp[0]]=gid
                gidToName[gid]=temp[0]
                gid=gid+1
        if(not temp[1] in nameDict):
                nameDict[temp[1]]=gid
                gidToName[gid]=temp[1]
                gid=gid+1

#print (nameDict)


dList=[]
#init d list with empty lists
for i in range(len(nameDict)):
        dList.append([0 for i in range(len(nameDict))])

#print (dList)

print("Converting to sparse mat")
thresh=0.6
for l in data:
        #print (l)
        thisList=[]
        temp=l.split("\t")
        thisRow=nameDict[temp[0]]
        thisCol=nameDict[temp[1]]
        if(float(temp[2])>thresh):
                edge=1
        else:
                edge=0
        thisList.append(edge)
        dList[thisRow][thisCol]=edge
        dList[thisCol][thisRow]=edge

#print (dList)

spArr=sparse.csr_matrix(dList)
print("Dim:"+str(spArr.shape))
print("Running MCL")
#choose inflation
maxMod=-1
bestInf=2
bestCluster=None
for inflation in [i / 10 for i in range(15, 26)]:
        result = mc.run_mcl(spArr, inflation=inflation)
        clusters = mc.get_clusters(result)
        Q = mc.modularity(matrix=result, clusters=clusters)
        print("inflation:", inflation, "modularity:", Q)
        if(Q>maxMod):
                maxMod=Q
                bestInf=inflation
                bestCluster=clusters
clusters=bestCluster
print("best inflation:"+str(bestInf))
print("best modularity:"+str(maxMod))

'''
result = mc.run_mcl(spArr)           # run MCL with default parameters
clusters = mc.get_clusters(result)
'''

#print(clusters)

print("Parsing results")

clusterwithNames={}
clid=0
for c in clusters:
        thisNames=[]
        for x in c:
                key=gidToName[x]
                thisNames.append(key)
        clusterwithNames[clid]=thisNames
        clid=clid+1

#print (clusterwithNames)

#write to file
print("writing to file...")
f=open(sys.argv[2]+"_inf_"+str(bestInf),"w")
f.write("ClusterID\tGenes")
for clid in clusterwithNames:
        thisList=clusterwithNames[clid]
        if len(thisList)>1:
                f.write("\n"+str(clid)+"\t"+";".join(thisList))
f.close()

print("writing to file...Done!!!")
