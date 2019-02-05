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
thresh=0.5
for l in data:
        #print (l)
        thisList=[]
        temp=l.split("\t")
        thisRow=nameDict[temp[0]]
        thisCol=nameDict[temp[1]]
        if(float(temp[2])>0.5):
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
'''maxMod=-1
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
print("writeing to file...")
f=open(sys.argv[2]+"_inf_2","w")
f.write("ClusterID\tGenes")
for clid in clusterwithNames:
        thisList=clusterwithNames[clid]
        if len(thisList)>1:
                f.write("\n"+str(clid)+"\t"+";".join(thisList))
f.close()

print("writing to file...Done!!!")

