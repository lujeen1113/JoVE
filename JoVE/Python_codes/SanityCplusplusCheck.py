#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 17:46:39 2021

@author: min
"""

from networkx.generators.random_graphs import erdos_renyi_graph
import networkx as nx
import copy
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import dimod,numpy,random

from dwave.embedding import chain_strength
import sys
flowm=[0, 0, 1, 0]
pos_matrix=[[10, 10],[10, 30],[30, 10],[30, 30]]
adjacency_matrix=[[0, 1, 1, 1],[1, 0, 1, 1],[1, 1, 0, 1],[1, 1, 1, 0]]

singlePathBuff=[]

def trimHJmatrix(H, J, fixed_dict):
    dictfixed=fixed_dict.copy()
    fixed_inds=list(dictfixed)
    print(fixed_dict)
    fixednum=len(fixed_inds)
    orig_size=len(H)
    for i in range(fixednum):
        fixed_inds[i]=fixed_inds[i]-1
    print(fixed_inds)
    orig_size = len(H)
    unfixed_inds = [i for i in range(orig_size) if not i in fixed_inds]
    num_fixed = len(fixed_inds)
    num_unfixed = orig_size - num_fixed
    if num_fixed == 0:
        return H.copy(), J.copy()
    Hcur = numpy.zeros(num_unfixed)
    Jcur = numpy.zeros((num_unfixed, num_unfixed))
    
    for j_unfixed in range(num_unfixed):
        j_unfixed_ind = unfixed_inds[j_unfixed]
        Hcur[j_unfixed] += H[j_unfixed_ind]
        for k_fixed in range(num_fixed):
            k_fixed_ind = fixed_inds[k_fixed]
            Hcur[j_unfixed] += (J[j_unfixed_ind, k_fixed_ind] + J[k_fixed_ind, j_unfixed_ind])*fixed_dict[k_fixed_ind+1]
        for k_unfixed in range(num_unfixed):
            k_unfixed_ind = unfixed_inds[k_unfixed]
            Jcur[j_unfixed, k_unfixed] = J[j_unfixed_ind, k_unfixed_ind]
    
    return Hcur, Jcur      

def cost(s, C, p):
    sT = numpy.transpose(s)
    #transpose for a column vector rather than for a row vector
    val = numpy.dot(s, numpy.transpose(numpy.dot(C,s))) + numpy.dot(p, sT)
    val = float(val)
    return val

def PathCollector(dest,flowm):
    global singlePathBuff
    pathCollector=[]
    sNodeID=[]
    for i in range(vertices):
        nodekey=i+1
        if nodekey!=dest and flowm[i]!=0:
            printAllPaths(nodekey,dest)
            #print(singlePathBuff)
            pathCollector.append(copy.copy(singlePathBuff))
            singlePathBuff=[]
            sNodeID.append(copy.copy(nodekey))
    return pathCollector,sNodeID
path_index=10
flag=0
srcID=0
pcollector=[0]*vertices
def printAllPaths(src,dest):
    global path_index,flag,srcID,count,pcollector
    srcID=src
    visited=[0]*vertices
    #pcollector=[0]*vertices
    path_index=1
    flag=0
    count=0
    printAllPathsUtil(src,visited,dest)
count=0   
def printAllPathsUtil(src,visited,dest):
    global path_index,flag,srcID,count,pcollector
    visited[src-1]=1
    #print(pcollector)
    pcollector[path_index-1]=src

    path_index=path_index+1
    if src==dest and count<maxflownum:
        #print(pcollector)
        singlePathBuff.append(copy.copy(pcollector))
        count=count+1
        flag=1
    elif count<maxflownum:
        for i in range(vertices):
            nodekey=i+1
            adjIdx=adjacency_matrix[src-1][nodekey-1]
            #print(adjacency_matrix)
            vIdx=visited[nodekey-1]
            #print(adjIdx)
            #print(vIdx)
            if adjIdx==1 and vIdx==0:
                
                printAllPathsUtil(nodekey, visited,dest)
            
    path_index=path_index-1
    visited[src-1]=0
    
   # if src==srcID:
   #     pcollector=[0]*vertices
    
#A=PathCollector(1)
#print(A)
import numpy
def getEdgeNum():
    count=0
    edge_count=numpy.zeros((vertices,vertices))
    for i in range(vertices):
        for j in range(i+1,vertices):
            if adjacency_matrix[i][j]==1:
               count=count+1
               edge_count[i][j]=count 
               edge_count[j][i]=count
    return edge_count,count
edge_count,edgenum=getEdgeNum() 


import copy
def getFlowAountperNode(PathCollector):
    pathnum=len(PathCollector)
    
    flowAmountBuff=[]
    
    for i in range(pathnum):
        B=PathCollector[i]
        num=len(B)
        flowAmountBuff.append(copy.copy(num))
    
    return flowAmountBuff
#flowm=[0,1,1,0]
PA,sNodeID=PathCollector(dest,flowm)
def getEdgeM(PA):
    PathArray=numpy.array(PA)
    #print(type(PathArray))
    #print(PathArray)
    pathnum=len(PathArray)
    #print(pathnum)
    flowPathBuff=getFlowAountperNode(PA)
    pathtotal=sum(flowPathBuff)
    vertices=pathnum/maxflownum+1
    edgem=numpy.zeros((vertices*vertices,pathtotal))
    flag=0
    count=0
    pathID=0
    for i in range(pathnum):
        pathIdx=i+1
        count=count+1
        flowamount=flowPathBuff[i]
        amount=0
        for j in range(i):
            amount=amount+flowPathBuff[j]
        pathID=amount
        for n in range(flowamount):
            pathID=pathID+1   
            flowIdx=n+1
            for j in range(vertices-1):
                nodeIdx=j+1
                k=PathArray[pathIdx-1][flowIdx-1][nodeIdx-1]
                p=PathArray[pathIdx-1][flowIdx-1][nodeIdx]
                #print(p)
                if p==1:   
                    #print('pathtotal'+str(pathtotal))
                    for o in range(pathtotal):
                        Idx2=o
                        # print(type(k))
                        kval=k-1
                        pval=p-1
                        Idx1=kval*vertices+pval
                        #print(edgem)
                        #print(Idx1)
                        #print(Idx2)
                        if edgem[Idx1][Idx2]==0:
                            edgem[Idx1][Idx2]=pathID
                            break
                    break
                else:
                    for o in range(pathtotal):
                        Idx2=o
                        #print(type(k))
                        kval=k-1
                        pval=p-1
                        Idx1=kval*vertices+pval
                        if edgem[Idx1][Idx2]==0:
                            edgem[Idx1][Idx2]=pathID
                            break
        #pathID=pathID+flowamount        
    return edgem                
        
#edgem=getEdgeM(PA)       
#EDGEMATRIX=edgem 
linkrate=5
d0=87.7085
lambdap=20
E=50
epson_fs=10
epson_mp=0.0013
t=5

dist_matrix=[[0,20, 20, 28.28],[20 ,0, 28.28, 20],[20, 28.28, 0, 20],[28.28, 20, 20, 0]]
def DistFun(n_i,n_j):
    dist=dist_matrix[n_i][n_j]
    return dist
def makeEffArray(edgem):
    flownum=0
    Idx=[]
    for i in range(vertices):
        if flowm[i] !=0:
            flownum=flownum+1
            Idx.append(copy.copy(i))
    edge_count,edgenum=getEdgeNum()        
    raterange=linkrate+1
    varnum=edgenum+flownum
    #print(varnum)
    var_size=numpy.zeros(varnum)
    flowPathBuff=getFlowAountperNode(PA)
    #print(flowPathBuff)
    pathnum=sum(flowPathBuff)
    msize1=0
    for i in range(flownum):
        NodeIdx=Idx[i]
        var_size[i]=flowPathBuff[i]
        msize1=msize1+flowPathBuff[i]
    msize2=0    
    for i in range(flownum,varnum,1):
        var_size[i]=raterange
        msize2=msize2+raterange
    
    size=msize1+msize2
    penalty=numpy.zeros(size,int)
    interaction=numpy.zeros((size,size),int)
    for i in range(vertices):
        NodeIdx1=i+1
        if NodeIdx1 !=dest:
            for j in range(vertices):
                NodeIdx2=j+1
                virtualedgenum=i*vertices+j
                #print('current edgem is:'+str(edgem))
                if edgem[virtualedgenum][0]!=0:
                    edgecount=edge_count[i][j]
                    for k in range(raterange):
                        l=k+1
                        m=msize1+(edgecount-1)*raterange+l
                        m=int(m)
                        #penalty[m-1]=(1-2*linkrate)*lambdap
                        penalty[m-1]=lambdap*((l-1)*(l-1)-2*linkrate*(l-1))
                    pathnum=len(edgem[0])	
                    #print("this penalty is :"+str(penalty))
                    for k in range(pathnum):
                        #print(edgem)
                        print("pathnum is :"+str(pathnum))
                        temp=edgem[virtualedgenum][k]
                        if temp==0:
                            break
                        else:
                            markflow,y,pathind,srnodeid=checkvalidity(temp,flowm,flowPathBuff)
                            if y!=0:
                                #print('the right path')
                                temp1=0
                                for p in range(markflow-1):
                                    temp1=temp1+flowPathBuff[p]
                                m=temp1+pathind
                                flowrate=flowm[srnodeid-1]
                                dist=DistFun(i,j)
                                val=0
                            
                                if dist<d0:
                                    val=E+epson_fs*dist*dist
                                else:
                                    val=E+epson_mp*dist*dist*dist*dist
                                m=int(m)
                                #print("val is :"+str(val))
                                #print("current m is :"+str(m))
                                #print("interval is :"+str(t))
                                #print("flowrate is :"+str(flowrate))
                                #print("linkrate is :"+str(linkrate))
                                penalty[m-1]=penalty[m-1]+t*flowrate*val+lambdap*flowrate*flowrate-2*linkrate*flowrate*lambdap
                                #print("current pval is :"+str(penalty[m-1]))
                                for o in range(raterange):
                                    l=o+1
                                    n1=pathnum+(edgecount-1)*raterange+l
                                    print("edgecount is:"+str(edgecount))
                                    print("l is:"+str(l))
                                    print("m is:"+str(m))
                                    print("n is:"+str(n1))
                                    n1=int(n1)
                                    print("int n is :"+str(n1))
                                    #print(*interaction)
                                    interaction[m-1][n1-1]=interaction[m-1][n1-1]+2*flowrate*(l-1)*lambdap
 
                                for p in range(k,pathnum):
                                    q=p+1
                                    temp2=edgem[virtualedgenum][q-1]
                                    if temp2==0:
                                       # print("the right path1")
                                        break
                                    else:
                                        markflow,y,pathind,srnodeid=checkvalidity(temp2,flowm,flowPathBuff)
                                        if y!=0:
                                            #print('the right path2')
                                            temp3=0
                                            for p1 in range(markflow-1):
                                                temp3=temp3+flowPathBuff[p1]
                                            n1=temp3+pathind
                                            flowrate2=flowm[srnodeid-1]
                                            n1=int(n1)
                                            #print("m is:"+str(m))
                                            #print("n is:"+str(n))
                                            interaction[m-1][n1-1]=interaction[m-1][n1-1]+2*flowrate*flowrate2*lambdap
                                            
    return penalty, interaction,var_size 

def checkvalidity(routeid,flowm,flowPathBuff):
    for i in range(vertices-1):
        amount=flowPathBuff[i]
        if routeid<=amount:
            NodeID=sNodeID[i]
            break
        else:
            routeid=routeid-amount
    k=NodeID
    pathid=routeid%amount
    if pathid==0:
        pathid=amount
    y=0
    nodepos=0
    flownum=0
    nodeid=k
    for i in range(vertices):
        if flowm[i]!=0:
            flownum=flownum+1
            if i==(k-1):
                nodepos=flownum
                y=1
    markflow=nodepos
    y=y
    pathind=pathid
    srnodeid=nodeid
    return markflow,y,pathind,srnodeid
#print('flowarray is:'+str(flowm))
print('PA is:'+str(PA))


PAnum=len(PA)
PAnullFlag=True
if PAnum>1:
    for i in range(PAnum):
        PAnullFlag=numpy.any(PA[i])
else:
    PAnullFlag=numpy.any(PA)
    
print(PAnullFlag)
edgem=getEdgeM(PA)       
EDGEMATRIX=edgem 
print(edgem)
msize=len(edgem)
from ez_domain_wall import make_domain_wall_encoding

def CombineTwoDict(fixed_dict,result_dict):
    fixed_num=len(fixed_dict)
    fixed_ind=list(fixed_dict)
    fixed_val=list(fixed_dict.values())
    fixed_ind = numpy.array(fixed_ind)

    fixed_val = numpy.array(fixed_val)

    unfixed_num=len(result_dict)
    unfixed_ind=list(result_dict)
    unfixed_val=list(result_dict.values())
    
    unfixed_ind = numpy.array(unfixed_ind)
    unfixed_val = numpy.array(unfixed_val)
    
    orig_num=fixed_num+unfixed_num
    
    orig_dict={}
    
    count=0
    for i in range(orig_num):
    	flag=False
    	for j in range(fixed_num):
    	    if i+1 == fixed_ind[j]:
    	    	orig_dict[i+1]=fixed_val[j]
    	    	del fixed_dict[fixed_ind[j]]
    	    	fixed_num=len(fixed_dict)
    	    	fixed_ind=list(fixed_dict)
    	    	fixed_ind = numpy.array(fixed_ind)
    	    	fixed_val=list(fixed_dict.values())
    	    	fixed_val=numpy.array(fixed_val)
    	    	flag=True
    	    	break
    	if flag == False:
    	   orig_dict[i+1]=unfixed_val[0]
    	   count=count+1
    	   del result_dict[count]
    	   unfixed_ind=list(result_dict)
    	   unfixed_ind=numpy.array(unfixed_ind)
    	   unfixed_val=list(result_dict.values())
    	   unfixed_val=numpy.array(unfixed_val)
    
    return orig_dict  
    
if PAnullFlag :
   print('PAFlagTrueLoop entered')
   penalty, interaction,var_size=makeEffArray(EDGEMATRIX)
   print(*penalty)
   print(*interaction)
   #numpy.savetxt('Documents/EE-Routing/penalty.txt',penalty)
   #numpy.savetxt('Documents/EE-Routing/interaction.txt',interaction)
   varnum=len(var_size)
   SingularRoute=0
   var_sizes=numpy.zeros(varnum,int)
   for i in range(varnum):
   	var_sizes[i]=int(var_size[i])
   	if var_size[i]==1:
            SingularRoute=1
            
   if SingularRoute ==0:
       J_core,H_core,J_prob,H_prob=make_domain_wall_encoding(var_sizes,penalty, interaction)
       #print("J_core is :"+str(J_core))
       #print("H_core is :"+str(H_core))
       #print("J_prob is :"+str(J_prob))
       #print("H_prob is :"+str(H_prob))
       constVal=100
       H=constVal*H_core+H_prob
       J=constVal*J_core+J_prob
       
       Hdict=dict(enumerate(H,1))
       size=len(H)
       Jdict={}
       for m in range(1,size+1):
           for n in range(1,size+1):
               if m != n:
                   Jdict[(m,n)]=J[m-1][n-1]
       #print("succeed before BQM")
       
       #print("Hdict is :"+str(Hdict))
       #print("Jdict is :" +str(Jdict))
       
       bqm=dimod.BinaryQuadraticModel.from_ising(Hdict,Jdict,offset=0.0)
       # fixed_dict=dimod.fix_variables(bqm)
       minene,fixed_dict=roof_duality(bqm)
       	
       size1=len(fixed_dict)
       
       if size1<size:
          Hcur,Jcur=trimHJmatrix(H,J,fixed_dict)
          Hdict=dict(enumerate(Hcur,1))
          size=len(Hcur)
          Jdict={}
          for m in range(1,size+1):
              for n in range(1,size+1):
                  if m != n:
                     Jdict[(m,n)]=Jcur[m-1][n-1]
          sampler = EmbeddingComposite(DWaveSampler(token='DEV-34ebdd33f77f73904ed58bac612191ccdce89841'))
          bqm=dimod.BinaryQuadraticModel.from_ising(Hdict,Jdict,offset=0.0)
          ChainVal=chain_strength.uniform_torque_compensation(bqm)
          token=100
          response=sampler.sample(bqm,num_reads=token,chain_strength=ChainVal,auto_scale=True)
          sol2=next(response.data(['sample']))[0]
          print(response)
          sample=list(sol2.values())
          SOL=CombineTwoDict(fixed_dict,sol2)
          SOL=list(SOL.values())
       else:
          SOL=list(fixed_dict.values())
       
       print(SOL)
    

