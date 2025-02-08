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
import math
import random
def segGraph(segNo,segAmount,theta_core_pos,ang):
    max_graph_size=segAmount
    radius=2000
    theta=0
    ang_seg=360/segAmount
    	
    pos_matrix=numpy.zeros((2,max_graph_size))
    adjacency_matrix=numpy.zeros((max_graph_size,max_graph_size))
    
    pos_matrix[0,0]=theta_core_pos[0,segNo]
    pos_matrix[1,0]=theta_core_pos[1,segNo]
    
    for i in range(1,max_graph_size):
    	r_radius=random.random()
    	temp_radius=r_radius*radius
    	
    	r_angular=random.random()
    	temp_angular=r_angular*ang_seg+ang-ang_seg
    	
    	x=temp_radius*math.cos(temp_angular)
    	y=temp_radius*math.sin(temp_angular)
    	
    	pos_matrix[0,i]=x
    	pos_matrix[1,i]=y
    	
    	
    
    for i in range(max_graph_size-1):
    	for j in range(i+1,max_graph_size,1):
    	    dice=random.random()
    	    if dice<0.7:
    	    	adjacency_matrix[i,j]=1
    	    	adjacency_matrix[j,i]=1
    
    
    return pos_matrix, adjacency_matrix
    
def ClusterComputing(opt):
    max_graph_size=6
    theta_core=[30,90,150,210,270,330]
    theta_bound=[60,120,180,240,300,360]
    theta_core_radius=1000
    theta_core_pos=numpy.zeros((2,max_graph_size))
    pi=3.1415926
    
    for i in range(max_graph_size):
    	temp=theta_core[i]*2*pi/360
    	theta_core[i]=temp
    
    for i in range(max_graph_size):
    	x=theta_core_radius*math.cos(theta_core[i])
    	y=theta_core_radius*math.sin(theta_core[i])
    	theta_core_pos[0,i]=x
    	theta_core_pos[1,i]=y
    
    segAmount = 6
    
    print("theta_core_pos is: "+str(theta_core_pos))
    flowm=numpy.zeros(segAmount)
    
    routing_SOL={}
    for i in range(segAmount):
        ang=theta_bound[i]
        pos_matrix,adjacency_matrix=segGraph(i,segAmount,theta_core_pos,ang)
        
        for j in range(1,segAmount):
            dice=random.random()
            if dice<0.5:
            	flowm[j]=1
        
        #print("pos_matrix is:"+str(pos_matrix))
        #print("adjacency_matrix is:"+str(adjacency_matrix))
        #pos_matrix=pos_matrix.transpose()
        routing_sol=SolveBundle(flowm,pos_matrix,adjacency_matrix,opt)
        print(routing_sol)


def SolveBundle(flowm,pos_matrix,adjacency_matrix,opt):
    maxflownum=3
    path_index=10
    flag=0
    srcID=0
    vertices=len(flowm)
    pcollector=[0]*vertices
    singlePathBuff=[]
    count=0 
    dest=1
    #dest=0
    linkrate=5
    #flowm=[0, 1, 1, 1]
    #pos_matrix=[[10, 10],[10, 30],[30, 10],[30, 30]]
    #adjacency_matrix=[[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 1],[1, 1, 1, 0]]
    edge_count,edgenum=getEdgeNum(adjacency_matrix)
    
    PA,sNodeID=PathCollector(dest,flowm,adjacency_matrix)
    
    #print("Path Array is:"+str(PA))
    #print("sNodeId is :"+str(sNodeID))
    
    linkrate=5
    d0=87.7085
    lambdap=20
    E=50
    epson_fs=10
    epson_mp=0.0013
    t=5

    #dist_matrix=[[0,20, 20, 28.28],[20 ,0, 28.28, 20],[20, 28.28, 0, 20],[28.28, 20, 20, 0]]
    
    PAnum=len(PA)
    PAnullFlag=True
    if PAnum>1:
       for i in range(PAnum):
           PAnullFlag=numpy.any(PA[i])
    else:
       PAnullFlag=numpy.any(PA)
    
    edgem=getEdgeM(PA,vertices)       
    EDGEMATRIX=edgem 
    msize=len(edgem)
    
    dist_matrix=numpy.zeros((vertices,vertices))
    
    for i in range(vertices):
    	for j in range(vertices):
    	    dist_x=pos_matrix[0,i]-pos_matrix[0,j]
    	    dist_y=pos_matrix[1,i]-pos_matrix[1,j]
    	    temp=math.sqrt(dist_x*dist_x+dist_y*dist_y)
    	    dist_matrix[i,j]=temp
    
    
    from ez_domain_wall import make_domain_wall_encoding
    
    if PAnullFlag :
       #print('PAFlagTrueLoop entered')
       #print("pos_matrix is :"+str(pos_matrix))
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
          #print("fixed_size is:"+str(size1))
          if size1<size:
             Hcur,Jcur=trimHJmatrix(H,J,fixed_dict)
             Hdict=dict(enumerate(Hcur,1))
             size=len(Hcur)
             Jdict={}
             for m in range(1,size+1):
                 for n in range(1,size+1):
                     if m != n:
                        Jdict[(m,n)]=Jcur[m-1][n-1]
             bqm=dimod.BinaryQuadraticModel.from_ising(Hdict,Jdict,offset=0.0)
             if opt==1:
             	sampler = EmbeddingComposite(DWaveSampler(token='DEV-34ebdd33f77f73904ed58bac612191ccdce89841'))
             	ChainVal=chain_strength.uniform_torque_compensation(bqm)
             	token=1000
             	response=sampler.sample(bqm,num_reads=token,chain_strength=ChainVal,auto_scale=True)
             	sol2=next(response.data(['sample']))[0]
             	sample=list(sol2.values())
             	SOL=CombineTwoDict(fixed_dict,sol2)
             	SOL=list(SOL.values())
             elif opt==2:
             	import cplex
             	bqm1=bqm
             	QUBO=bqm1.to_qubo()
             	QUBO_dw=numpy.zeros((size,size))
             	for m in range(1,size+1):
                   for n in range(1,size+1):
                    	if m>=n:
                          QUBO_dw[m-1][n-1]=QUBO[0][(m,n)]
		num=len(QUBO_dw)
		H=numpy.zeros(num)
		J=numpy.zeros((num,num))
            	for i in range(num):
                    H[i]=QUBO_dw[i][i]
                    for j in range(i+1,num):
                        J[i][j]=0.5*QUBO_dw[i][j]
                        J[j][i]=0.5*QUBO_dw[i][j]
                          
                p=cplex.Cplex()
                p.objective.set_sense(p.objective.sense.minimize)
                obj=H
                ub=[1]*num
            
            	for i in range(num):
                    if i==0:
                       typestr="I"
                    else:
                        typestr=typestr+"I"


            	p.variables.add(obj=obj,ub=ub,types=typestr)

                qmat=[1]*num

                for i in range(num):
                    rowarray=[1,1]
                    indexarray=[]
                    valarray=[]
                    for j in range(num):
                        if J[i][j]!=0:
                           indexarray.append(copy.copy(j))
                           valarray.append(copy.copy(J[i][j]))
                           rowarray[0]=indexarray
                           rowarray[1]=valarray
                           qmat[i]=rowarray

                p.objective.set_quadratic(qmat)
                beforetime=time.time()
                p.solve()
                aftertime=time.time()
                #print('Cplex time is :'+str(aftertime-beforetime))
                sol=p.solution
                SOL4=numpy.zeros(num)
  
                for i in range(num):
                    SOL4[i]=sol.get_values(i)

                d=SOL4.dot(QUBO_dw)
                ans4=d.dot(numpy.transpose(SOL4))
                SOL=SOL4
             
          else:
             SOL=list(fixed_dict.values())
       
          #print(SOL)
          return SOL
          
#SolverBundle() 
#ClusterComputing()
import time
import threading   
def SpeedSlideWindow(weight):
    rate=80
    #meter per second
    t_start=time.time()
    count=2
    for i in range(count):
    	ClusterComputing(2)
    t_stop=time.time()
    
    t_lap=t_stop-t_start
    print("elapsed time is :"+str(t_lap))
    
SpeedSlideWindow(1)
    
    
