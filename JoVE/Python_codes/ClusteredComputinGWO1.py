import threading
import time
from networkx.generators.random_graphs import erdos_renyi_graph
import networkx as nx
import copy
from dwave.preprocessing.lower_bounds import roof_duality
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import dimod
import numpy
import random

from dwave.embedding import chain_strength
import sys
from SolveBundleCalls import getEdgeNum, PathCollector, getEdgeM, makeEffArray, trimHJmatrix, CombineTwoDict
import math
import random


def segGraph(segNo, segAmount, theta_core_pos, ang):
    max_graph_size = 6
    radius = 100
    theta = 0
    ang_seg = 360/segAmount

    pos_matrix = numpy.zeros((2, max_graph_size))
    adjacency_matrix = numpy.zeros((max_graph_size, max_graph_size))

    #pos_matrix[0, 0] = theta_core_pos[0, segNo]
    #pos_matrix[1, 0] = theta_core_pos[1, segNo]
    
    centr_x=theta_core_pos[0,segNo]
    centr_y=theta_core_pos[1,segNo]

    dist=radius
    min_indx=0
    temp_dist=radius
    for i in range(max_graph_size):
        r_radius = random.random()
        temp_radius = r_radius*radius

        r_angular = random.random()
        temp_angular = r_angular*ang_seg+ang-ang_seg

        x = temp_radius*math.cos(temp_angular)
        y = temp_radius*math.sin(temp_angular)
        
        temp_dist=math.sqrt((x-centr_x)*(x-centr_x)+(y-centr_y)*(y-centr_y))
        
        if temp_dist<=dist:
        	dist=temp_dist
        	min_indx=i

        pos_matrix[0, i] = x
        pos_matrix[1, i] = y
        
       
    temp_x=pos_matrix[0,0]
    temp_y=pos_matrix[1,0]
    
    pos_matrix[0,0]=pos_matrix[0,min_indx]
    pos_matrix[1,0]=pos_matrix[1,min_indx]
    
    pos_matrix[0,min_indx]=temp_x
    pos_matrix[1,min_indx]=temp_y
        
    for i in range(max_graph_size-1):
        for j in range(i+1, max_graph_size, 1):
            dice = random.random()
            if dice < 1:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1

    return pos_matrix, adjacency_matrix
    
    
def segGraph_update(segNo, segAmount, theta_core_pos, ang,whole_buff):
    max_graph_size = 6
    radius = 50
    theta = 0
    ang_seg = 360/max_graph_size

    pos_matrix = numpy.zeros((2, max_graph_size))
    adjacency_matrix = numpy.zeros((max_graph_size, max_graph_size))

    #pos_matrix[0, 0] = theta_core_pos[0, segNo]
    #pos_matrix[1, 0] = theta_core_pos[1, segNo]
    
    centr_x=theta_core_pos[0,segNo]
    centr_y=theta_core_pos[1,segNo]

    dist=radius
    min_indx=0
    temp_dist=radius
    for i in range(max_graph_size):
        r_radius = random.random()
        temp_radius = r_radius*radius

        r_angular = random.random()
        temp_angular = r_angular*ang_seg+ang-ang_seg

        x = temp_radius*math.cos(temp_angular)
        y = temp_radius*math.sin(temp_angular)
        
        temp_dist=math.sqrt((x-centr_x)*(x-centr_x)+(y-centr_y)*(y-centr_y))
        
        if temp_dist<=dist:
        	dist=temp_dist
        	min_indx=i

        pos_matrix[0, i] = x
        pos_matrix[1, i] = y
        
       
    temp_x=pos_matrix[0,0]
    temp_y=pos_matrix[1,0]
    
    pos_matrix[0,0]=pos_matrix[0,min_indx]
    pos_matrix[1,0]=pos_matrix[1,min_indx]
    
    pos_matrix[0,min_indx]=temp_x
    pos_matrix[1,min_indx]=temp_y
        
    for i in range(max_graph_size-1):
        for j in range(i+1, max_graph_size, 1):
            dice = random.random()
            if dice < 1 and whole_buff[segNo*max_graph_size+i]>0 and whole_buff[segNo*max_graph_size+j]>0:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1

    return pos_matrix, adjacency_matrix


def ClusterComputing(opt,whole_buff,pos_matrix_buff,adjacency_matrix_buff):
    max_graph_size = 6
    theta_core = [30, 90, 150, 210, 270, 330]
    theta_bound = [60, 120, 180, 240, 300, 360]
    theta_core_radius = 25
    theta_core_pos = numpy.zeros((2, max_graph_size))
    pi = 3.1415926
    
    d0 = 87.7085
    lambdap = 20
    
    E=0.000000005
    epson_fs=0.0000000010
    epson_mp=0.0013*0.000000000001
    bitsnum=4000

    for i in range(max_graph_size):
        temp = theta_core[i]*2*pi/360
        theta_core[i] = temp

    for i in range(max_graph_size):
        x = theta_core_radius*math.cos(theta_core[i])
        y = theta_core_radius*math.sin(theta_core[i])
        theta_core_pos[0, i] = x
        theta_core_pos[1, i] = y

    segAmount = 6

    print("theta_core_pos is: "+str(theta_core_pos))
    flowm = numpy.zeros(segAmount)
    
    routing_PATH = []
    totaltime=0
    
    for i in range(segAmount):
        for j in range(1, segAmount):
            dice = random.random()
            if dice < 0.5 and whole_buff[i*segAmount+j]>0:
                flowm[j] = 1
                #1000bits/second

        print("pos_matrix is:"+str(pos_matrix[i]))
        print("adjacency_matrix is:"+str(adjacency_matrix[i]))
        print("flowmatrix is:"+str(flowm))
        # pos_matrix=pos_matrix.transpose()
       
        path=SolveBundle(flowm,pos_matrix[i],adjacency_matrix[i],opt)
        #totaltime=totaltime+t_val
        
        num=len(path)
        
        for j in range(num):
            step_count=0
            for k in range(segAmount):
            	if path[j][k]!=0:
            	   step_count=step_count+1
            	   if step_count==1:
            	   	node1=path[j][k]
            	   	node2=path[j][k+1]
            	   	dist=math.sqrt((pos_matrix[i][k][1]-pos_matrix[i][k+1][1])*(pos_matrix[i][k][1]-pos_matrix[i][k+1][1])-(pos_matrix[i][k][2]-pos_matrix[i][k+1][2])*(pos_matrix[i][k][2]-pos_matrix[i][k+1][2]))
            	   	if dist<=87.7085:
            	   		ene=bitsnum*E+bitsnum*epson_fs*dist*dist
            	   	else:
            	   		ene=bitsnum*E+bitsnum*epson_mf*dist*dist*dist*dist

            	   	node_idx=i*segAmount+node1
            	        whole_buff[node_idx]=whole_buff[node_idx]-ene

            	        if whole_buff[node_idx]<0:
            	        	whole_buff[node_idx]=-1
            	        
            	   	else:
            	        	node1=path[j][k]
            	   		node2=path[j][k+1]
            	   		dist=math.sqrt((pos_matrix[i][k][1]-pos_matrix[i][k+1][1])*(pos_matrix[i][k][1]-pos_matrix[i][k+1][1])-(pos_matrix[i][k][2]-pos_matrix[i][k+1][2])*(pos_matrix[i][k][2]-pos_matrix[i][k+1][2]))
            	   		if dist<=87.7085:
            	   			ene=bitsnum*E*2+bitsnum*epson_fs*dist*dist
            	        	else:
            	        		ene=bitsnum*E*2+bitsnum*epson_mf*dist*dist*dist*dist
            	        
            	        	node_idx=i*segAmount+node1
            	        	whole_buff[node_idx]=whole_buff[node_idx]-ene
            	        
            	        	if whole_buff[node_idx]<0:
            	        		whole_buff[node_idx]=-1
            	        
             whole_buff[i*segAmount]=whole_buff[i*segAmount]-bitsnum*E
             if whole_buff[node_idx]<0:
            	        	whole_buff[node_idx]=-1    
        routing_PATH.append(copy.copy(path))
        print(routing_Path)
        
        
    print("routing sol bundle is:"+str(routing_Path))
  
    #print("total computational time spent is:"+str(totaltime))

    return whole_buff

def SolveBundle(flowm, pos_matrix, adjacency_matrix, opt):
    maxflownum = 3
    path_index = 10
    flag = 0
    srcID = 0
    vertices = len(flowm)
    pcollector = [0]*vertices
    singlePathBuff = []
    count = 0
    dest = 1
    # dest=0
    linkrate = 5
    #flowm=[0, 1, 1, 1]
    #pos_matrix=[[10, 10],[10, 30],[30, 10],[30, 30]]
    #adjacency_matrix=[[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 1],[1, 1, 1, 0]]
    edge_count, edgenum = getEdgeNum(adjacency_matrix)

    PA, sNodeID = PathCollector(dest, flowm, adjacency_matrix)

    print("Path Array is:"+str(PA))
    print("sNodeId is :"+str(sNodeID))

    linkrate = 5
    d0 = 87.7085
    lambdap = 20
    
    E=0.000000005
    epson_fs=0.0000000010
    epson_mp=0.0013*0.000000000001
    
    t = 4

    #dist_matrix=[[0,20, 20, 28.28],[20 ,0, 28.28, 20],[20, 28.28, 0, 20],[28.28, 20, 20, 0]]

    PAnum = len(PA)
    PAnullFlag = True
    if PAnum > 1:
        for i in range(PAnum):
            PAnullFlag = numpy.any(PA[i])
    else:
        PAnullFlag = numpy.any(PA)

    edgem = getEdgeM(PA, vertices)
    EDGEMATRIX = edgem
    msize = len(edgem)

    dist_matrix = numpy.zeros((vertices, vertices))

    for i in range(vertices):
        for j in range(vertices):
            dist_x = pos_matrix[0, i]-pos_matrix[0, j]
            dist_y = pos_matrix[1, i]-pos_matrix[1, j]
            temp = math.sqrt(dist_x*dist_x+dist_y*dist_y)
            dist_matrix[i, j] = temp

    from ez_domain_wall import make_domain_wall_encoding

    if PAnullFlag:
        #print('PAFlagTrueLoop entered')
        #print("pos_matrix is :"+str(pos_matrix))
        penalty, interaction, var_size = makeEffArray(
            EDGEMATRIX, PA, vertices, flowm, adjacency_matrix, linkrate, sNodeID, dist_matrix)
        print("penality is :" +str(penalty))
        # print(*interaction)
    # numpy.savetxt('Documents/EE-Routing/penalty.txt',penalty)
    # numpy.savetxt('Documents/EE-Routing/interaction.txt',interaction)
        varnum = len(var_size)
        SingularRoute = 0
        var_sizes = numpy.zeros(varnum, int)
        for i in range(varnum):
            var_sizes[i] = int(var_size[i])
            if var_size[i] == 1:
                SingularRoute = 1

        if SingularRoute == 0:
            print("not a singular route")
            J_core, H_core, J_prob, H_prob = make_domain_wall_encoding(
                var_sizes, penalty, interaction)
            print("makeDWencoding finir")
            constVal = 500
            H = constVal*H_core+H_prob
            J = constVal*J_core+J_prob

            Hdict = dict(enumerate(H, 1))
            size = len(H)
            Jdict = {}
            for m in range(1, size+1):
                for n in range(1, size+1):
                    if m != n:
                        Jdict[(m, n)] = J[m-1][n-1]

            bqm = dimod.BinaryQuadraticModel.from_ising(
                Hdict, Jdict, offset=0.0)
            
            minene, fixed_dict = roof_duality(bqm)

            size1 = len(fixed_dict)
            print("fixed_size is:"+str(size1))
            if size1 < size:
                Hcur, Jcur = trimHJmatrix(H, J, fixed_dict)
                Hdict = dict(enumerate(Hcur, 1))
                size = len(Hcur)
                Jdict = {}
                for m in range(1, size+1):
                    for n in range(1, size+1):
                        if m != n:
                            Jdict[(m, n)] = Jcur[m-1][n-1]
                bqm = dimod.BinaryQuadraticModel.from_ising(
                    Hdict, Jdict, offset=0.0)
                if opt == 1:
                    sampler = EmbeddingComposite(DWaveSampler(
                        token='DEV-34ebdd33f77f73904ed58bac612191ccdce89841'))
                    ChainVal = chain_strength.uniform_torque_compensation(bqm)
                    token = 10
                    response = sampler.sample(
                        bqm, num_reads=token, chain_strength=ChainVal, auto_scale=True)
                    timeval=response.info["timing"]['qpu_sampling_time']/1000000
                    sol2 = next(response.data(['sample']))[0]
                    sample = list(sol2.values())
                    SOL = CombineTwoDict(fixed_dict, sol2)
                    SOL = list(SOL.values())
                elif opt == 2:
                    import cplex
                    QUBO = bqm.to_qubo()
                    QUBO_dw = numpy.zeros((size, size))
                    for m in range(1, size+1):
                        for n in range(1, size+1):
                            if m >= n:
                                QUBO_dw[m-1][n-1] = QUBO[0][(m, n)]
                    H = numpy.zeros(size)
                    J = numpy.zeros((size, size))
                    for i in range(size):
                        H[i] = QUBO_dw[i][i]
                        for j in range(i+1, size):
                            J[i][j] = 0.5*QUBO_dw[i][j]
                            J[j][i] = 0.5*QUBO_dw[i][j]

                    p = cplex.Cplex()
                    p.objective.set_sense(p.objective.sense.minimize)
                    obj = H
                    ub = [1]*size

                    for i in range(size):
                        if i == 0:
                            typestr = "I"
                        else:
                            typestr = typestr+"I"

                    p.variables.add(obj=obj, ub=ub, types=typestr)

                    qmat = [1]*size

                    for i in range(size):
                        rowarray = [1, 1]
                        indexarray = []
                        valarray = []
                        for j in range(size):
                            if J[i][j] != 0:
                                indexarray.append(copy.copy(j))
                                valarray.append(copy.copy(J[i][j]))
                                rowarray[0] = indexarray
                                rowarray[1] = valarray
                                qmat[i] = rowarray

                    p.objective.set_quadratic(qmat)
                    beforetime = time.time()
                    p.solve()
                    aftertime = time.time()
                    timeval=aftertime-beforetime
                    sol = p.solution
                    SOL4 = numpy.zeros(size)

                    for i in range(size):
                        SOL4[i] = sol.get_values(i)

                    d = SOL4.dot(QUBO_dw)
                    ans4 = d.dot(numpy.transpose(SOL4))
                    SOL = SOL4
                
                   

            else:
                SOL = list(fixed_dict.values())
                timeval=0.0
		
        print("sol is:"+str(SOL))
        print("time is:"+str(timeval))
        
        num_path=len(sNodeID)
        pathIdx_buff=numpy.zeros(num_path)
        path_buff=[]
        
        for i in range(num_path):
            if opt==1:
            	bit1=SOL[0+i*1]
            	bit2=SOL[1+i*1]
            
            if opt==2:
            	bit1=SOL[0+i*1]*2-1
            	bit2=SOL[1+i*1]*2-1
            print(type(bit1))
            print(str(bit1)+str(bit2))
            
            bit1=int(bit1)
            bit2=int(bit2)
            if bit1==1 and bit2==1:
               pathIndx=0
            if bit1==-1 and bit2==1:
               pathIdx=1
            if bit1==-1 and bit2==-1:
               pathIdx=2
            else:
               pathIdx=3
               print("invalid solution")
            	
            pathIdx_buff[i]=pathIdx
             
            print("path Idx is:"+str(pathIdx))  
            path_temp=PA[i][pathIdx]		
            path_buff.append(copy.copy(path_temp))
            
        
        
        return path_buff



# SolverBundle()
#ClusterComputing(2,1)


def SpeedSlideWindow(opt):
    rate = 80
    # meter per second
    t_start = time.time()
    aus_num=6
    whole_buff=[]
    for i in range(aus_num):
    	ene_buff=numpy.ones(aus_num)*0.5
    	whole_buff.append(copy.copy(ene_buff))
    run_flag=True
    count=0
    FND_num=0
    HND_num=0
    LND_num=0
    half_count=0
    all_count=0
    
    max_graph_size = 6
    theta_core = [30, 90, 150, 210, 270, 330]
    theta_bound = [60, 120, 180, 240, 300, 360]
    theta_core_radius = 50
    theta_core_pos = numpy.zeros((2, max_graph_size))
    pi = 3.1415926

    for i in range(max_graph_size):
        temp = theta_core[i]*2*pi/360
        theta_core[i] = temp

    for i in range(max_graph_size):
        x = theta_core_radius*math.cos(theta_core[i])
        y = theta_core_radius*math.sin(theta_core[i])
        theta_core_pos[0, i] = x
        theta_core_pos[1, i] = y

    segAmount = 6

    print("theta_core_pos is: "+str(theta_core_pos))
    flowm = numpy.zeros(segAmount)
    
    routing_SOL = []
    totaltime=0
    pos_matrix_buff=[]
    adjacency_matrix_buff=[]
    
    for i in range(segAmount):
        ang = theta_bound[i]
        pos_matrix, adjacency_matrix = segGraph(
            i, segAmount, theta_core_pos, ang)
    	pos_matrix_buff.append(copy.copy(pos_matrix))
    	adjacency_matrix_buff.append(copy.copy(adjacency_matrix))
    
    
    while run_flag==True:
    	whole_buff_update = ClusterComputing(opt,whole_buff,pos_matrix_buff,adjacency_matrix_buff)
    	count=count+1
    	for i in range(aus_num):
    		for j in range(aus_num):
    			if whole_buff_update[i][j] <=0:
    				FND_num=count
    				half_count=half_count+1	
    				all_count=all_count+1
    				
    				if half_count==aus_num*aus_num*0.5:
    					HND_num=count
    				if all_count==aus_num*aus_num:
    					LND_num=count
    					
    	pos_matrix_buff=[]
        adjacency_matrix_buff=[]
    
    	for i in range(segAmount):
            ang = theta_bound[i]
            pos_matrix, adjacency_matrix = segGraph_update(
            i, segAmount, theta_core_pos, ang,whole_buff_update)
    	    pos_matrix_buff.append(copy.copy(pos_matrix))
    	    adjacency_matrix_buff.append(copy.copy(adjacency_matrix))
    					
    	whole_buff=whole_buff_update				
        	
    #ClusterComputing(1)
    t_stop = time.time()

    t_lap = t_stop-t_start
    print("elapsed time is :"+str(t_lap))

SpeedSlideWindow(opt)
