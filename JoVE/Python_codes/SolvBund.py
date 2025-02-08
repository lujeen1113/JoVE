from networkx.generators.random_graphs import erdos_renyi_graph
import networkx as nx
import copy
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import dimod,numpy,random

from dwave.embedding import chain_strength
import sys
from SolveBundleCalls import getEdgeNum,PathCollector,getEdgeM,makeEffArray,trimHJmatrix,CombineTwoDict
def SolverBundle(h):
    maxflownum=3
    path_index=10
    flag=0
    srcID=0
    vertices=4
    pcollector=[0]*vertices
    singlePathBuff=[]
    count=0 
    dest=1
    linkrate=5
    flowm=[0, 1, 1, 1]
    pos_matrix=[[10, 10],[10, 30],[30, 10],[30, 30]]
    adjacency_matrix=[[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 1],[1, 1, 1, 0]]
    edge_count,edgenum=getEdgeNum(adjacency_matrix)
    
    PA,sNodeID=PathCollector(dest,flowm,adjacency_matrix)
    
    print("Path Array is:"+str(PA))
    print("sNodeId is :"+str(sNodeID))
    
    linkrate=5
    d0=87.7085
    lambdap=20
    E=50
    epson_fs=10
    epson_mp=0.0013
    t=5

    dist_matrix=[[0,20, 20, 28.28],[20 ,0, 28.28, 20],[20, 28.28, 0, 20],[28.28, 20, 20, 0]]
    
    PAnum=len(PA)
    PAnullFlag=True
    if PAnum>1:
       for i in range(PAnum):
           PAnullFlag=numpy.any(PA[i])
    else:
       PAnullFlag=numpy.any(PA)
    
    print(PAnullFlag)
    edgem=getEdgeM(PA,vertices)       
    EDGEMATRIX=edgem 
    print(edgem)
    msize=len(edgem)
    
    from ez_domain_wall import make_domain_wall_encoding
    
    if PAnullFlag :
       print('PAFlagTrueLoop entered')
       penalty, interaction,var_size=makeEffArray(EDGEMATRIX,PA,vertices,flowm,adjacency_matrix,linkrate,sNodeID,dist_matrix)
       #print(*penalty)
       #print(*interaction)
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
       
          constVal=500
          H=constVal*H_core+H_prob
          J=constVal*J_core+J_prob
       
          Hdict=dict(enumerate(H,1))
          size=len(H)
          Jdict={}
          for m in range(1,size+1):
              for n in range(1,size+1):
                  if m != n:
                     Jdict[(m,n)]=J[m-1][n-1]
       
       
          bqm=dimod.BinaryQuadraticModel.from_ising(Hdict,Jdict,offset=0.0)
       
          minene,fixed_dict=roof_duality(bqm)
       	
          size1=len(fixed_dict)
          print("fixed_size is:"+str(size1))
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
             token=1000
             response=sampler.sample(bqm,num_reads=token,chain_strength=ChainVal,auto_scale=True)
             sol2=next(response.data(['sample']))[0]
             print(response)
             sample=list(sol2.values())
             SOL=CombineTwoDict(fixed_dict,sol2)
             SOL=list(SOL.values())
          else:
             SOL=list(fixed_dict.values())
       
          print(SOL)
          return SOL
          
#SolverBundle()    
