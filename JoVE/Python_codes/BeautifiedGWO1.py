#!/usr/bin/python
# -*- coding: utf-8 -*-

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
from SolveBundleCalls import getEdgeNum, PathCollector, getEdgeM, \
    makeEffArray, trimHJmatrix, CombineTwoDict
import math
import random


def segGraph(
    segNo,
    segAmount,
    theta_core_pos,
    ang,
    ):

    max_graph_size = 6
    radius = 100
    theta = 0
    ang_seg = 360 / segAmount

    pos_matrix = numpy.zeros((2, max_graph_size))
    adjacency_matrix = numpy.zeros((max_graph_size, max_graph_size))

    # pos_matrix[0, 0] = theta_core_pos[0, segNo]
    # pos_matrix[1, 0] = theta_core_pos[1, segNo]

    centr_x = theta_core_pos[0, segNo]
    centr_y = theta_core_pos[1, segNo]

    dist = radius
    min_indx = 0
    temp_dist = radius
    pi = 3.1415926
    for i in range(max_graph_size):
        r_radius = random.random()
        temp_radius = r_radius * radius

        r_angular = random.random()
        temp_angular = (r_angular * ang_seg + ang - ang_seg) * 2 * pi \
            / 360

        x = temp_radius * math.cos(temp_angular)
        y = temp_radius * math.sin(temp_angular)

        temp_dist = math.sqrt((x - centr_x) * (x - centr_x) + (y
                              - centr_y) * (y - centr_y))

        if temp_dist <= dist:
            dist = temp_dist
            min_indx = i

        pos_matrix[0, i] = x
        pos_matrix[1, i] = y

    temp_x = pos_matrix[0, 0]
    temp_y = pos_matrix[1, 0]

    pos_matrix[0, 0] = pos_matrix[0, min_indx]
    pos_matrix[1, 0] = pos_matrix[1, min_indx]

    pos_matrix[0, min_indx] = temp_x
    pos_matrix[1, min_indx] = temp_y

    for i in range(max_graph_size - 1):
        for j in range(i + 1, max_graph_size, 1):
            dice = random.random()
            if dice < 1:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1
    ID = 0
    return (pos_matrix, adjacency_matrix, ID)


def segGraph_update(
    segNo,
    segAmount,
    theta_core_pos,
    ang,
    whole_buff,
    pos_matrix,
    ):

    max_graph_size = 6
    radius = 100
    theta = 0

    adjacency_matrix = numpy.zeros((max_graph_size, max_graph_size))

    centr_x = theta_core_pos[0, segNo]
    centr_y = theta_core_pos[1, segNo]

    # #print('centr_x is:'+str(centr_x))
    # #print('centr_y is:'+str(centr_y))

    dist = radius
    min_indx = 0
    temp_dist = radius
    NotFoundFlag = True

    weight = 0
    dist_dic = {}
    dist_list = []

    ene_dic = {}
    ene_list = []
    dist_array=[]
    count=0
    
    ene=0
    for i in range(max_graph_size):
        
        if whole_buff[segNo][i] > 0:
            x = pos_matrix[0][i]
            y = pos_matrix[1][i]
            count = count + 1
            temp_dist = math.sqrt((x - centr_x) * (x - centr_x) + (y
                                  - centr_y) * (y - centr_y))
            dist_array.append(copy.copy(temp_dist))
            dist_list.append(copy.copy(temp_dist))
            dist_dic[count] = temp_dist
            ene_list.append(copy.copy(whole_buff[segNo][i]))
            ene_dic[count] = whole_buff[segNo][i]
            temp_ene=whole_buff[segNo][i]
            
            if temp_ene >= ene:
            # if temp_weight>weight:
               ene = temp_ene
             #   weight=temp_weight
               min_indx = i
               NotFoundFlag=False
                # #print('temp_dist is:'+str(temp_dist))
               print('min_idx is:'+str(i))
        else:
        	dist_array.append(copy.copy(radius))
      
    #calculate the difference between the sorted energy and distance indx
    dist_list.sort()
    ene_list.sort()

    invalid_node_num = len(dist_list)

    decr = -invalid_node_num
    print(dist_dic)
    print(ene_dic)
    print(invalid_node_num)
    ene_val_finir=0.5
    dist_val_finir=100
    for i in range(invalid_node_num):
        dist_ID = i
        dist_val = dist_list[i]
        dist_idx=1
        for (key, value) in dist_dic.items():
            if value == dist_val:
                dist_idx = key

        ene_val = ene_dic[dist_idx]

        for j in range(invalid_node_num):
            if ene_list[j] == ene_val:
                ene_ID = j
        temp_decr = ene_ID - dist_ID

        if temp_decr >= decr:
            decr = temp_decr
            #min_indx = dist_idx
            ene_val_finir=ene_val
            dist_val_finir=dist_val
            NotFoundFlag = False

    #print('min_idx is:'+str(min_indx))
    
    for i in range(max_graph_size):
    	if whole_buff[segNo][i]==ene_val_finir and dist_array[i]==dist_val_finir:
    	   #min_indx=i
    	   print('calc1 finished')
    #####################################################################
    
    	

    if NotFoundFlag:
        print('not found!')
    
    for i in range(max_graph_size - 1):
        for j in range(i + 1, max_graph_size, 1):
            dice = random.random()
            if whole_buff[segNo][i] > 0 and whole_buff[segNo][j] > 0:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1

    return (adjacency_matrix, int(min_indx))


def ClusterComputing(
    opt,
    whole_buff,
    pos_matrix_buff,
    adjacency_matrix_buff,
    destID,
    segAmount,
    ):

    max_graph_size = 6
    theta_core = [
        30,
        90,
        150,
        210,
        270,
        330,
        ]
    theta_bound = [
        60,
        120,
        180,
        240,
        300,
        360,
        ]
    theta_core_radius = 50
    theta_core_pos = numpy.zeros((2, max_graph_size))
    pi = 3.1415926

    d0 = 87.7085
    lambdap = 20

    E = 50 * 1e-09
    epson_fs = 1e-09 * 10 * 0.001
    epson_mp = 0.0013 * 0.000000000001 * 0.001
    bitsnum = 4000

    for i in range(max_graph_size):
        temp = theta_core[i] * 2 * pi / 360
        theta_core[i] = temp

    for i in range(max_graph_size):
        x = theta_core_radius * math.cos(theta_core[i])
        y = theta_core_radius * math.sin(theta_core[i])
        theta_core_pos[0, i] = x
        theta_core_pos[1, i] = y

    max_graph_size = 6

    flowm = numpy.zeros(max_graph_size)

    routing_PATH = []
    totaltime = 0

    for i in range(segAmount):
        flowm = numpy.zeros(max_graph_size)
        for j in range(max_graph_size):
            dice = random.random()
            if dice <= 1 and whole_buff[i][j] > 0 and j != destID[i]:
                flowm[j] = 1

        if numpy.all(flowm == 0) == True:
            for k in range(max_graph_size):
                ID = k
                if ID != destID[i] and whole_buff[i][ID] > 0:
                    flowm[ID] = 1
                    break

                # 1000bits/second
        # #print('whole_buff is:'+str(whole_buff))

        part_buff = whole_buff[i]

        # ##print('pos_matrix is:' + str(pos_matrix_buff[i]))
        print('adjacency_matrix is:' + str(adjacency_matrix_buff[i]))
        print('flowmatrix is:' + str(flowm))

        # pos_matrix=pos_matrix.transpose()

        dest_id = int(destID[i])

        print('destID is:'+str(destID))
        print('segNO is:'+str(i))

        if numpy.all(adjacency_matrix_buff[i] == 0) != True:
            path = SolveBundle(flowm, pos_matrix_buff[i],
                               adjacency_matrix_buff[i], opt, dest_id)
        else:
            single = False
            count = 0
            for p in range(max_graph_size):
                if part_buff[p] == -1:
                    count = count + 1
            if count == 5:
                single = True
            ID = 0
            if single == True:
                for p in range(max_graph_size):
                    if part_buff[p] > 0:
                        ID = p

            dist = math.sqrt((pos_matrix_buff[i][0][ID] - 0)
                             * (pos_matrix_buff[i][0][ID] - 0)
                             + (pos_matrix_buff[i][1][ID] - 0)
                             * (pos_matrix_buff[i][1][ID] - 0))

            if dist <= 87.7085:
                whole_buff[i][ID] = whole_buff[i][ID] - bitsnum * E \
                    - bitsnum * epson_fs * dist * dist
            else:
                whole_buff[i][ID] = whole_buff[i][ID] - bitsnum * E \
                    - bitsnum * epson_mp * dist * dist * dist * dist

            if whole_buff[i][ID] < 0:
                whole_buff[i][ID] = -1

            path = 1

        # totaltime=totaltime+t_val
        # #print('path is:'+str(path))

        if path == 1:
            num = 0
        else:
            num = len(path)

        for j in range(num):
            step_count = 0
            for k in range(max_graph_size):
                if path[j][k] != destID[i]:
                    step_count = step_count + 1
                    if step_count == 1:
                        node1 = path[j][k] - 1
                        node2 = path[j][k + 1] - 1
                        dist = math.sqrt((pos_matrix_buff[i][1][k]
                                - pos_matrix_buff[i][1][k + 1])
                                * (pos_matrix_buff[i][1][k]
                                - pos_matrix_buff[i][1][k + 1])
                                + (pos_matrix_buff[i][0][k]
                                - pos_matrix_buff[i][0][k + 1])
                                * (pos_matrix_buff[i][0][k]
                                - pos_matrix_buff[i][0][k + 1]))
                        if dist <= 87.7085:
                            ene = bitsnum * E + bitsnum * epson_fs \
                                * dist * dist
                        else:
                            ene = bitsnum * E + bitsnum * epson_mp \
                                * dist * dist * dist * dist

                        # node_idx = i * segAmount + node1
                        # #print('i is:'+str(i))
                        # #print('node is:'+str(node1))

                        whole_buff[i][node1] = whole_buff[i][node1] \
                            - ene

                        if whole_buff[i][node1] < 0:
                            whole_buff[i][node1] = -1
                        else:

                            node1 = path[j][k] - 1
                            node2 = path[j][k + 1] - 1
                            dist = math.sqrt((pos_matrix_buff[i][1][k]
                                    - pos_matrix_buff[i][1][1 + k])
                                    * (pos_matrix_buff[i][1][k]
                                    - pos_matrix_buff[i][1][1 + k])
                                    + (pos_matrix_buff[i][0][k]
                                    - pos_matrix_buff[i][0][k + 1])
                                    * (pos_matrix_buff[i][0][k]
                                    - pos_matrix_buff[i][0][k + 1]))
                            if dist <= 87.7085:
                                ene = bitsnum * E * 2 + bitsnum \
                                    * epson_fs * dist * dist
                            else:
                                ene = bitsnum * E * 2 + bitsnum \
                                    * epson_mp * dist * dist * dist \
                                    * dist

                            # node_idx = i * segAmount + node1

                            whole_buff[i][node1] = whole_buff[i][node1] \
                                - ene

                            if whole_buff[i][node1] < 0:
                                whole_buff[i][node1] = -1

            # #print(dest_id)
            # #print(pos_matrix_buff)
            # ##print(pos_matrix_buff[i][dest_id])

            whole_buff[i][dest_id] = whole_buff[i][dest_id] -bitsnum \
                * E
            dist = math.sqrt((pos_matrix_buff[i][1][dest_id] - 0)
                             * (pos_matrix_buff[i][1][dest_id] - 0)
                             + (pos_matrix_buff[i][0][dest_id] - 0)
                             * (pos_matrix_buff[i][0][dest_id] - 0))
            if dist <= 87.7085:
                whole_buff[i][dest_id] = whole_buff[i][dest_id] \
                    - bitsnum * E - bitsnum * epson_fs * dist * dist
            else:
                whole_buff[i][dest_id] = whole_buff[i][dest_id] \
                    - bitsnum * E - bitsnum * epson_mp * dist * dist \
                    * dist * dist

            if whole_buff[i][dest_id] < 0:
                whole_buff[i][dest_id] = -1

        routing_PATH.append(copy.copy(path))

    # #print('routing sol bundle is:' + str(routing_PATH))

    # ##print("whole_buff:"+str(whole_buff))

    return (whole_buff, routing_PATH)


def SolveBundle(
    flowm,
    pos_matrix,
    adjacency_matrix,
    opt,
    ID,
    ):

    maxflownum = 3
    path_index = 10
    flag = 0
    srcID = 0
    vertices = len(flowm)
    pcollector = [0] * vertices
    singlePathBuff = []
    count = 0
    dest = ID + 1

    # dest=0

    linkrate = 5

    # flowm=[0, 1, 1, 1]
    # pos_matrix=[[10, 10],[10, 30],[30, 10],[30, 30]]
    # adjacency_matrix=[[0, 0, 1, 1],[0, 0, 1, 1],[1, 1, 0, 1],[1, 1, 1, 0]]

    (edge_count, edgenum) = getEdgeNum(adjacency_matrix)

    (PA, sNodeID) = PathCollector(dest, flowm, adjacency_matrix)

    # #print('Path Array is:' + str(PA))
    # #print('sNodeId is :' + str(sNodeID))

    linkrate = 5
    d0 = 87.7085
    lambdap = 20

    #E = 0.000000005
    #epson_fs = 1e-09
    #epson_mp = 0.0013 * 0.000000000001
    
    d0 = 87.7085
    lambdap = 20

    E = 50 * 1e-09
    epson_fs = 1e-09 * 10 * 0.001
    epson_mp = 0.0013 * 0.000000000001 * 0.001
    bitsnum = 4000

    t = 1000

    # dist_matrix=[[0,20, 20, 28.28],[20 ,0, 28.28, 20],[20, 28.28, 0, 20],[28.28, 20, 20, 0]]

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
            dist_x = pos_matrix[0, i] - pos_matrix[0, j]
            dist_y = pos_matrix[1, i] - pos_matrix[1, j]
            temp = math.sqrt(dist_x * dist_x + dist_y * dist_y)
            dist_matrix[i, j] = temp

    from ez_domain_wall import make_domain_wall_encoding

    if PAnullFlag:

        # ##print('PAFlagTrueLoop entered')
        # ##print("pos_matrix is :"+str(pos_matrix))

        (penalty, interaction, var_size) = makeEffArray(
            EDGEMATRIX,
            PA,
            vertices,
            flowm,
            adjacency_matrix,
            linkrate,
            sNodeID,
            dist_matrix,
            )

        # #print('penality is :' + str(penalty))

        # ##print(*interaction)
    # numpy.savetxt('Documents/EE-Routing/penalty.txt',penalty)
    # numpy.savetxt('Documents/EE-Routing/interaction.txt',interaction)

        varnum = len(var_size)
        SingularRoute = 0
        var_sizes = numpy.zeros(varnum, int)
        path_buff = []
        for i in range(varnum):
            var_sizes[i] = int(var_size[i])
            if var_size[i] == 1:
                SingularRoute = 1
                path = PA[0][i]
                path_buff.append(copy.copy(path))

                # #print('path_buff is:'+str(path_buff))

        if SingularRoute == 0:

            # #print('not a singular route')

            (J_core, H_core, J_prob, H_prob) = \
                make_domain_wall_encoding(var_sizes, penalty,
                    interaction)

            # #print('makeDWencoding finir')

            constVal = 500
            H = constVal * H_core + H_prob
            J = constVal * J_core + J_prob

            Hdict = dict(enumerate(H, 1))
            size = len(H)
            Jdict = {}
            for m in range(1, size + 1):
                for n in range(1, size + 1):
                    if m != n:
                        Jdict[(m, n)] = J[m - 1][n - 1]

            bqm = dimod.BinaryQuadraticModel.from_ising(Hdict, Jdict,
                    offset=0.0)

            (minene, fixed_dict) = roof_duality(bqm)

            size1 = len(fixed_dict)

            # #print('fixed_size is:' + str(size1))

            if size1 < size:
                (Hcur, Jcur) = trimHJmatrix(H, J, fixed_dict)
                Hdict = dict(enumerate(Hcur, 1))
                size = len(Hcur)
                Jdict = {}
                for m in range(1, size + 1):
                    for n in range(1, size + 1):
                        if m != n:
                            Jdict[(m, n)] = Jcur[m - 1][n - 1]
                bqm = dimod.BinaryQuadraticModel.from_ising(Hdict,
                        Jdict, offset=0.0)
                if opt == 1:
                    sampler = \
                        EmbeddingComposite(DWaveSampler(token='DEV-34ebdd33f77f73904ed58bac612191ccdce89841'
                            ))
                    ChainVal = \
                        chain_strength.uniform_torque_compensation(bqm)
                    token = 10
                    response = sampler.sample(bqm, num_reads=token,
                            chain_strength=ChainVal, auto_scale=True)
                    timeval = response.info['timing'
                            ]['qpu_sampling_time'] / 1000000
                    sol2 = next(response.data(['sample']))[0]
                    sample = list(sol2.values())
                    SOL = CombineTwoDict(fixed_dict, sol2)
                    SOL = list(SOL.values())
                elif opt == 2:
                    import cplex
                    QUBO = bqm.to_qubo()
                    QUBO_dw = numpy.zeros((size, size))
                    for m in range(1, size + 1):
                        for n in range(1, size + 1):
                            if m >= n:
                                QUBO_dw[m - 1][n - 1] = QUBO[0][(m, n)]
                    H = numpy.zeros(size)
                    J = numpy.zeros((size, size))
                    for i in range(size):
                        H[i] = QUBO_dw[i][i]
                        for j in range(i + 1, size):
                            J[i][j] = 0.5 * QUBO_dw[i][j]
                            J[j][i] = 0.5 * QUBO_dw[i][j]

                    p = cplex.Cplex()
                    p.objective.set_sense(p.objective.sense.minimize)
                    obj = H
                    ub = [1] * size

                    for i in range(size):
                        if i == 0:
                            typestr = 'I'
                        else:
                            typestr = typestr + 'I'

                    p.variables.add(obj=obj, ub=ub, types=typestr)

                    qmat = [1] * size

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
                    timeval = aftertime - beforetime
                    sol = p.solution
                    SOL4 = numpy.zeros(size)

                    for i in range(size):
                        SOL4[i] = 2 * sol.get_values(i) - 1

                    SOL = dict()
                    for (index, value) in enumerate(SOL4):
                        SOL[index + 1] = value

                    # #print('dictionary is:'+str(SOL))

                    SOL = CombineTwoDict(fixed_dict, SOL)
                    SOL = list(SOL.values())
            else:

                SOL = list(fixed_dict.values())
                timeval = 0.0

            # #print('sol is:' + str(SOL))
            # #print('time is:' + str(timeval))

            num_path = len(sNodeID)
            pathIdx_buff = numpy.zeros(num_path)
            path_buff = []

            for i in range(num_path):
                bitnum = int(var_size[i] - 1)

                if bitnum == 2:
                    bit1 = int(SOL[0 + i * bitnum])
                    bit2 = int(SOL[1 + i * bitnum])
                    if bit1 == 1 and bit2 == 1:
                        pathIdx = 0
                    elif bit1 == -1 and bit2 == 1:
                        pathIdx = 1
                    elif bit1 == -1 and bit2 == -1:
                        pathIdx = 2
                    else:
                        pathIdx = 3
                elif bitnum == 1:

                        # #print('invalid solution')

                    bit = int(SOL[0 + i * bitnum])
                    if bit == 1:
                        pathIdx = 0
                    elif bit == -1:
                        pathIdx = 1
                    else:
                        pathIdx = 2

                        # #print('invalid solution')

                pathIdx_buff[i] = pathIdx

                # #print('path Idx is:' + str(pathIdx))
                # #print('PA is:'+str(PA))
                # #print('i is:'+str(i))
                # #print('adj_matrix is:'+str(adjacency_matrix))

                path_temp = PA[i][pathIdx]
                path_buff.append(copy.copy(path_temp))

        return path_buff


# SolverBundle()
# ClusterComputing(2,1)

def SpeedSlideWindow(opt):
    rate = 80

    # meter per second

    t_start = time.time()
    aus_num = 6
    whole_buff = []
    for i in range(aus_num):
        ene_buff = numpy.ones(aus_num) * 0.5
        whole_buff.append(copy.copy(ene_buff))

    # ##print('whole_buff is'+str(whole_buff))

    run_flag = True
    count = 0
    FND_num = 0
    HND_num = 0
    LND_num = 0
    half_count = 0
    all_count = 0

    d0 = 87.7085
    lambdap = 20

    E = 50 * 1e-09
    epson_fs = 1e-09 * 10 * 0.001
    epson_mp = 0.0013 * 0.000000000001 * 0.001
    bitsnum = 4000

    max_graph_size = 6
    theta_core = [
        30,
        90,
        150,
        210,
        270,
        330,
        ]
    theta_bound = [
        60,
        120,
        180,
        240,
        300,
        360,
        ]
    theta_core_radius = 50
    theta_core_pos = numpy.zeros((2, max_graph_size))
    pi = 3.1415926

    for i in range(max_graph_size):
        temp = theta_core[i] * 2 * pi / 360
        theta_core[i] = temp

    for i in range(max_graph_size):
        x = theta_core_radius * math.cos(theta_core[i])
        y = theta_core_radius * math.sin(theta_core[i])
        theta_core_pos[0, i] = x
        theta_core_pos[1, i] = y

    segAmount = 6

    # #print('theta_core_pos is: ' + str(theta_core_pos))

    flowm = numpy.zeros(segAmount)

    routing_SOL = []
    totaltime = 0
    pos_matrix_buff = []
    adjacency_matrix_buff = []
    destID = numpy.zeros(segAmount)
    for i in range(segAmount):
        ang = theta_bound[i]
        (pos_matrix, adjacency_matrix, ID) = segGraph(i, segAmount,
                theta_core_pos, ang)
        pos_matrix_buff.append(copy.copy(pos_matrix))
        adjacency_matrix_buff.append(copy.copy(adjacency_matrix))
        destID[i] = ID

    obsolete_idx = []

    # destID=numpy.zeros(segAmount)

    tig_flag = False
   
    while run_flag==True:

        # ##print('FND_NUM is :' + str(FND_num))
        # ##print('HND_NUM is :' + str(HND_num))
        # ##print('LND_NUM is :' + str(LND_num))
        # #print('segAmount is:'+str(segAmount))

        (whole_buff_update, sol_bundle) = ClusterComputing(
            opt,
            whole_buff,
            pos_matrix_buff,
            adjacency_matrix_buff,
            destID,
            segAmount,
            )
        print('whole_buff is:'+str(whole_buff_update))
        count = count + 1
        delete_ID = []
        segInvalid_flag = False
        for i in range(segAmount):
            count1 = 0

            ndim = len(whole_buff_update)
            

            if numpy.all(whole_buff_update[i] < 0) == True and ndim > 1:
                delete_ID.append(copy.copy(i))
                segInvalid_flag = True
                count1=(max_graph_size-ndim)*6+6
                
            elif numpy.all(whole_buff_update[i]<0)==False and ndim>1:
                count1=(6-ndim)*6
                for j in range(max_graph_size):
                	if whole_buff_update[i][j]<0:
                		count1=count1+1
                
            elif ndim == 1 and numpy.all(whole_buff_update < 0) == True:
                delete_ID.append(copy.copy(i))
                segInvalid_flag = True
                count1=36
            elif ndim==1 and numpy.all(whole_buff_update<0)==False:
                 count1=30
            	 for j in range(max_graph_size):
                	if whole_buff_update[i][j]<0:
                		count1=count1+1
            
       if count1==1:
          FND_num=count
                
       if count1==18:
          HND_num=count
         
       if count1=36:
          LND_num=count
          break	

        # #print('invalid SegList is:'+str(delete_ID))
        print(whole_buff_update)
        if segInvalid_flag == True and ndim > 1:
            for o in range(len(delete_ID)):
            	idx = delete_ID[o]

            	theta_bound = numpy.delete(theta_bound, idx)

            # #print('theta_POS is:'+str(theta_core_pos))

            	theta_core_pos = numpy.delete(theta_core_pos, idx, axis=1)

            # #print('theta_POS after is:'+str(theta_core_pos))
                # ##print(pos_matrix_buff)

            	pos_matrix_buff = numpy.delete(pos_matrix_buff, idx, axis=0)

                # ##print('before'+str(whole_buff_update))

            	whole_buff_update = numpy.delete(whole_buff_update, idx,
                    axis=0)
        elif segInvalid_flag == True and ndim == 1:

                # ##print('after'+str(whole_buff_update))

            whole_buff_update = []

        segAmount = len(whole_buff_update)

        # #print('Amount is:'+str(segAmount))

        adjacency_matrix_buff = []
        #if segAmount == 0:
        #   LND_num = count

        #   print('all segments depleted')
        #   print('the final whole_buff is:'+str(whole_buff_update))

        #   run_flag = False
           #break

        destID = numpy.zeros(segAmount)
        for i in range(segAmount):
            ang = theta_bound[i]
            (adjacency_matrix, ID) = segGraph_update(
                i,
                segAmount,
                theta_core_pos,
                ang,
                whole_buff_update,
                pos_matrix_buff[i],
                )

            adjacency_matrix_buff.append(copy.copy(adjacency_matrix))
            destID[i] = int(ID)
        whole_buff = whole_buff_update
        #if numpy.all(adjacency_matrix_buff == 0) == True:
        #    run_flag = False

    # ClusterComputing(1)

    t_stop = time.time()

    t_lap = t_stop - t_start
    print('FND_NUM is :' + str(FND_num))
    print('HND_NUM is :' + str(HND_num))
    print('LND_NUM is :' + str(LND_num))
    print(obsolete_idx)


SpeedSlideWindow(2)
