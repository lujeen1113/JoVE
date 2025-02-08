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
from SolveBundleCallsv2 import getEdgeNum, PathCollector, getEdgeM, \
    makeEffArray, trimHJmatrix, CombineTwoDict
import math
import random
from operator import itemgetter

def segGraph(
    segNo,
    segAmount,
    theta_core_pos,
    ang,
    ):

    max_graph_size = 6
    radius = 50
    theta = 0
    ang_seg = 360 / segAmount

    pos_matrix = numpy.zeros((2, max_graph_size + 1))
    adjacency_matrix = numpy.zeros((max_graph_size + 1, max_graph_size
                                   + 1))

    # pos_matrix[0, 0] = theta_core_pos[0, segNo]
    # pos_matrix[1, 0] = theta_core_pos[1, segNo]

    centr_x = theta_core_pos[0, segNo]
    centr_y = theta_core_pos[1, segNo]

    dist = radius
    min_indx = 0
    temp_dist = radius
    pi = 3.1415926
    for i in range(1, max_graph_size + 1):
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

    for i in range(max_graph_size + 1):
        for j in range(i + 1, max_graph_size + 1, 1):
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
    adjacency_matrix = numpy.zeros((max_graph_size + 1, max_graph_size
                                   + 1))
    for i in range(1, max_graph_size + 1):
        adjacency_matrix[0, i] = 1
        adjacency_matrix[i, 0] = 1

    for i in range(1, max_graph_size + 1):
        for j in range(i + 2, max_graph_size + 1, 1):
            dice = random.random()
            if whole_buff[segNo][i] > 0 and whole_buff[segNo][j] > 0:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1

    return (adjacency_matrix, 0)


def ClusterComputing(
    opt,
    whole_buff,
    pos_matrix,
    adjacency_matrix,
    destID,
    CHID_buff,
    ):

    max_graph_size = 6

    pos_matrix_temp = numpy.zeros((1 + len(CHID_buff), 2))

    for i in range(len(CHID_buff)):
        ID = CHID_buff[i]
        pos_temp_x = pos_matrix[ID][0]
        pos_temp_y = pos_matrix[ID][1]

        for j in range(1, 1 + len(CHID_buff)):
            if pos_matrix_temp[j][0] == 0 and pos_matrix_temp[j][1] \
                == 0:
                indx = j
                break

        pos_matrix_temp[indx][0] = pos_temp_x
        pos_matrix_temp[indx][1] = pos_temp_y

    pi = 3.1415926

    d0 = 87.7085
    lambdap = 20

    E = 50 * 1e-09
    epson_fs = 0.000000000001 * 10
    epson_mp = 0.0013 * 0.000000000001
    bitsnum = 4000

    routing_PATH = []

    for i in range(1):
        flowm = numpy.zeros(len(CHID_buff) + 1)
        for j in range(1, len(CHID_buff) + 1):
            dice = random.random()
            ID = CHID_buff[j-1]
            if dice <= 1 and whole_buff[ID] > 0:
                flowm[j] = 1

        print('adjacency_matrix is:' + str(adjacency_matrix))
        print('flowmatrix is:' + str(flowm))

        dest_id = 0

        # #print('destID is:' + str(destID[i]))
        # #print 'segNO is:' + str(i)

        edge_count = 0
        singlepath_flag = False
        for o in range(len(CHID_buff) + 1):
            for p in range(len(CHID_buff) + 1):
                if adjacency_matrix[o][p] == 1:
                    edge_count = edge_count + 1

        if numpy.all(adjacency_matrix[i] == 0) != True \
            and edge_count > 2:
            ##print('to enter SolvBundle')
            path = SolveBundle(flowm, pos_matrix_temp,
                               adjacency_matrix, opt, dest_id)
            ##print('solved solution is:' + str(path))
        elif numpy.all(adjacency_matrix[i] == 0) != True \
            and edge_count == 2:
            for l in range(len(CHID_buff) + 1):
                if flowm[l] == 1 and l != dest_id:
                    src_id = l
            path_sub = numpy.zeros(len(CHID_buff) + 1)
            path_sub[0] = src_id
            path_sub[1] = dest_id
            path = []
            path.append(copy.copy(path_sub))
            singlepath_flag = True
        ##print(type(path))
        print('path is:'+str(path))
        if path == 1:
            num = 0
        else:
            num = len(path)

        # ##print('path is'+str(path))

        for j in range(num):
            step_count = 0
            for k in range(len(CHID_buff) + 1):
                ID = int(dest_id)
                if path[j][k] != ID:
                    step_count = step_count + 1
                    if step_count == 1:

                        # #print('ENERGYDEDUCTIONBEGINS')

                        node1 = path[j][k] - 2
                        node2 = path[j][k + 1] - 2
                        dist = math.sqrt((pos_matrix_temp[1][k]
                                - pos_matrix_temp[1][k + 1])
                                * (pos_matrix_temp[1][k]
                                - pos_matrix_temp[1][k + 1])
                                + (pos_matrix_temp[0][k]
                                - pos_matrix_temp[0][k + 1])
                                * (pos_matrix_temp[0][k]
                                - pos_matrix_temp[0][k + 1]))
                        if dist <= 87.7085:
                            ene = bitsnum * E + bitsnum * epson_fs \
                                * dist * dist
                        else:
                            ene = bitsnum * E + bitsnum * epson_mp \
                                * dist * dist * dist * dist

                        node1 = int(node1)
                        node2 = int(node2)
                        ##print(path[j])
                        ##print(node1)
                        ##print(CHID_buff)
                        ID1 = CHID_buff[node1]
                        #ID2 = CHID_buff[node2]
                        whole_buff[ID1] = whole_buff[ID1] - ene

                        if whole_buff[ID1] < 0:
                            whole_buff[ID1] = -1

        routing_PATH.append(copy.copy(path))

    return whole_buff


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

    # ####print('Path Array is:' + str(PA))
    # ####print('sNodeId is :' + str(sNodeID))

    linkrate = 5
    d0 = 87.7085
    lambdap = 20

    # E = 0.000000005
    # epson_fs = 1e-09
    # epson_mp = 0.0013 * 0.000000000001

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
    #print(pos_matrix)
    for i in range(vertices):
        for j in range(vertices):
            dist_x = pos_matrix[i, 0] - pos_matrix[j, 0]
            dist_y = pos_matrix[i, 1] - pos_matrix[j, 1]
            temp = math.sqrt(dist_x * dist_x + dist_y * dist_y)
            dist_matrix[i, j] = temp

    from ez_domain_wall import make_domain_wall_encoding
    #print('PAnullFlag is:' + str(PAnullFlag))
    if PAnullFlag:

        # #####print('PAFlagTrueLoop entered')
        # #####print("pos_matrix is :"+str(pos_matrix))

        (penalty, interaction, var_size) = makeEffArray(
            EDGEMATRIX,
            PA,
            vertices,
            flowm,
            adjacency_matrix,
            linkrate,
            sNodeID,
            dist_matrix,
            dest,
            )

        # ####print('penality is :' + str(penalty))

        # #####print(*interaction)
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

                # ####print('path_buff is:'+str(path_buff))

        if SingularRoute == 0:

            # ####print('not a singular route')

            (J_core, H_core, J_prob, H_prob) = \
                make_domain_wall_encoding(var_sizes, penalty,
                    interaction)

            # ####print('makeDWencoding finir')

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

            # ####print('fixed_size is:' + str(size1))

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

                    # #print('timing is:'+str(response.info['timing']))

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

                    # ####print('dictionary is:'+str(SOL))

                    SOL = CombineTwoDict(fixed_dict, SOL)
                    SOL = list(SOL.values())
            else:

                SOL = list(fixed_dict.values())
                timeval = 0.0

            # ####print('sol is:' + str(SOL))
            # ####print('time is:' + str(timeval))

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

                        # ####print('invalid solution')

                    bit = int(SOL[0 + i * bitnum])
                    if bit == 1:
                        pathIdx = 0
                    elif bit == -1:
                        pathIdx = 1
                    else:
                        pathIdx = 2

                        # ####print('invalid solution')

                pathIdx_buff[i] = pathIdx

                # ####print('path Idx is:' + str(pathIdx))
                # ####print('PA is:'+str(PA))
                # ####print('i is:'+str(i))
                # ####print('adj_matrix is:'+str(adjacency_matrix))

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
    max_graph_size=6
    addon=0
    whole_buff = numpy.zeros(aus_num * (max_graph_size+addon + 1))

    for i in range(aus_num * (aus_num + 1)):
        if i % (aus_num+1) == 0:
            whole_buff[i] = -1
        elif (i - 1) % 2 == 0:
            whole_buff[i] = 0.5
        else:
            whole_buff[i] = 1

    run_flag = True
    count = 0
    FND_num = 0
    HND_num = 0
    LND_num = 0
    half_count = 0
    all_count = 0

    d0 = 87.7085
    lambdap = 20

    E = 50 * 0.000000001
    epson_fs = 0.001 * 10*0.000000001
    epson_mp = 0.0013 *0.000000001*0.001
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

    # ####print('theta_core_pos is: ' + str(theta_core_pos))

    routing_SOL = []
    totaltime = 0
    pos_matrix = []
    adjacency_matrix_buff = []
    destID = numpy.zeros(segAmount)
    pos_matrix=numpy.zeros((aus_num*(max_graph_size+1),2))
    for i in range(segAmount):
        ang = theta_bound[i]
        (pos_matrix1, adjacency_matrix, ID) = segGraph(i, segAmount,
                theta_core_pos, ang)
        name = 'posdata' + str(i) + '.txt'
        pos_matrix_temp = numpy.loadtxt(name)
        col = [0, 0]
        pos_matrix_temp = numpy.insert(pos_matrix_temp, 0, col, axis=1)
        #print('pos_matrix_temp:'+str(pos_matrix_temp))
        for j in range(max_graph_size+1):
        	idx=i*(max_graph_size+1)+j
        	pos_matrix[idx][0]=pos_matrix_temp[0][j]
        	pos_matrix[idx][1]=pos_matrix_temp[1][j]
       	
        adjacency_matrix_buff.append(copy.copy(adjacency_matrix))
        destID[i] = ID

    obsolete_idx = []

    # destID=numpy.zeros(segAmount)
    #pos_matrix=pos_matrix.transpose()
    FND_flag = False
    HND_flag = False
    LND_flag=False
    DeadRatio_1500R = False
    DR_1500 = 0
    seg_count = 0

    LND_num_1 = 0

    cl_buff = numpy.ones(aus_num * (max_graph_size + 1)) * 60

    run_flag = True
    
    P = 0.173
 
    round_limit = int(1 / P)
    alpha = 1
    m = 0.5
    E_in = 0.5

    P_nrm = P / (1 + alpha * m)
    P_adv = P * (1 + alpha) / (1 + alpha * m)

    #b = random.random()
    b=0.8
    u = random.random()
    v = 1 - u

    P_bnrm = b * P_nrm
    P_badv = b * P_adv
    
    count=0

    while run_flag == True:

        print('FND_NUM is :' + str(FND_num))
        print('HND_NUM is :' + str(HND_num))
        print('LND_NUM is :' + str(LND_num))
        print('segAmount is:'+str(segAmount))
        
        count=count+1

        d_max = 0
        d_sum = 0
        gcount = 0
        e_sum = 0

        # ##print('round_buff is:'+str(round_buff))

        # #print 'whole buff is:' + str(whole_buff)

        p_time_start = time.time()

        CH_Amount = 0
        CHIdx_buff = numpy.ones(aus_num * (max_graph_size + 1))

        for i in range(aus_num * (max_graph_size + 1)):
            if i % (max_graph_size + 1) != 0:
                if whole_buff[i] > 0:
                    gcount=gcount+1
                    x_temp = pos_matrix[i][0]
                    y_temp = pos_matrix[i][1]
                    d_temp = math.sqrt(x_temp * x_temp + y_temp
                            * y_temp)
                    d_sum = d_sum + d_temp
                    e_sum = e_sum + whole_buff[i]
                    if d_temp > d_max:
                        d_max = d_temp
            else:
                CHIdx_buff[i] = 66

        d_ave = d_sum / gcount
        d_bs = d_ave
        e_ave = e_sum / gcount
        CHID_buff = []

        T_n = P / (1 - P * (count % (1 / P)))
        CH_Amount = len(CHID_buff)
        
        count_num=0
        cl_flag=True
        for i in range(aus_num*(max_graph_size+1)):
        	if i%(max_graph_size+1)!=0 and whole_buff[i]>0:
        		count_num=count_num+1
        
        if count_num<=max_graph_size:
           cl_flag=False
           CH_BUFF=[]
           ch_buff=[]
           for i in range(aus_num*(max_graph_size+1)):
           	if i%(max_graph_size+1)!=0 and whole_buff[i]>0:
           	   ch_buff.append(copy.copy(i))
           CH_BUFF.append(copy.copy(ch_buff))
        
        
        count2=0
        while CH_Amount < max_graph_size and count_num>max_graph_size:
            count2=count2+1
            CHID_buff = []
            CHIdx_buff = numpy.zeros(aus_num * (max_graph_size + 1))

            for i in range(aus_num * (max_graph_size + 1)):
            	if i % (max_graph_size + 1) != 0:
            		if whole_buff[i] > 0:
                    		gcount=gcount+1
                    		x_temp = pos_matrix[i][0]
                    		y_temp = pos_matrix[i][1]
                    		d_temp = math.sqrt(x_temp * x_temp + y_temp
                            * y_temp)
                    		d_sum = d_sum + d_temp
                    		e_sum = e_sum + whole_buff[i]
                    		if d_temp > d_max:
                       		 d_max = d_temp
            	else:
                	CHIdx_buff[i] = 66
            for i in range(aus_num * (max_graph_size + 1)):
                if i % (max_graph_size + 1) != 0 and whole_buff[i] > 0:
                    if (i-1) % 2 == 0:
                        E_start = 0.5
                    else:
                        E_start = 1

                    E_cur = whole_buff[i]

                    x_temp = pos_matrix[i][0]
                    y_temp = pos_matrix[i][1]
                    d_cur = math.sqrt(x_temp * x_temp + y_temp
                                * y_temp)

                    T_n_mnrm = P_bnrm * (u * (E_cur / E_start) + v
                                * (d_cur / d_max)) / (1 - P_nrm
                                * (count % (1 / P_nrm)))

                    T_n_madv = P_badv * (u * (E_cur / E_start) + v
                                * (d_cur / d_max)) / (1 - P_adv
                                * (count % (1 / P_adv)))

                    T_n_nrm1 = T_n_mnrm * (E_cur / E_start)
                    T_n_adv1 = T_n_madv * (E_cur / E_start)

                    T_n_nrm2 = T_n_nrm1 * e_ave
                    T_n_adv2 = T_n_adv1 * e_ave
                        ##print('d_cur:'+str(d_cur))
                        ##print('d_ave:'+str(d_ave))
                        ##print('i is:'+str(i))
                        ##print(pos_matrix)
                    T_n_fnrm = T_n_nrm2 * (1 / d_ave + 1 / d_cur)
                    T_n_fadv = T_n_adv2 * (1 / d_ave + 1 / d_cur)

                    threshold_nrm = random.random()*1.5
                    threshold_adv = random.random()*1.5

                        ##print('mnrm is:' + str(T_n_mnrm))
                        ##print('madv is:' + str(T_n_madv))

                    # T_n_mnrm=0
                    # T_n_madv=0

                    if (i-1) % 2 == 0:
                       if threshold_nrm < T_n:
                          CHIdx_buff[i] = i
                          cl_buff[i] = 0
                          ID = i
                          CHID_buff.append(copy.copy(ID))
                       
                               
                    else:
                        if threshold_adv < T_n:
                           CHIdx_buff[i] = i
                           cl_buff[i] = 0
                           ID = i
                           CHID_buff.append(copy.copy(ID))
                           
                elif whole_buff[j]<0 and j%(max_graph_size+1)!=0:
            	     CHIdx_buff[j]=-1                

            CH_Amount = len(CHID_buff)

            chid_buff = CHID_buff

         
        if cl_flag==True:
            CH_BUFF =[]
            for i in range(len(CHID_buff)):
            	temp=[]
            	CH_BUFF.append(copy.copy(temp))
        	
        #print('CHID_buff loop is:'+str(chid_buff))
        #print('CHIdx_buff loop is:'+str(CHIdx_buff))
        #print('whole_buff loop is:'+str(whole_buff)) 
            cl_size=numpy.ones(len(CHID_buff))
            for j in range(aus_num * (max_graph_size + 1)):
                if j % (max_graph_size + 1) != 0 and CHIdx_buff[j]== 0 and whole_buff[j]>0:
                   x = pos_matrix[j][0]
                   y = pos_matrix[j][1]

                   dist = 20000
                   ID = 0
                   flag = True
                   dist_array=numpy.zeros((len(CHID_buff),2))
                   for p in range(len(CHID_buff)):
                       id_ch = CHID_buff[p]

                       x_ch = pos_matrix[id_ch][0]
                       y_ch = pos_matrix[id_ch][1]
                       dist_temp = math.sqrt((x - x_ch) * (x
                                    - x_ch) + (y - y_ch) * (y - y_ch))
                       dist_array[p][0]=dist_temp
                       dist_array[p][1]=p
                   
                   dist_array=list(dist_array)
                   ##print(dist_array)
                   dist_array=list(dist_array)
                   dist_array.sort(key=itemgetter(0))
                   ##print(numpy.array(dist_array))
                   
                   for p in range(len(CHID_buff)):
                   	CH_no=int(dist_array[p][1])
                   	CH_ID=CHID_buff[CH_no]
                   	if cl_size[CH_no]<max_graph_size:
                   		if cl_size[CH_no]==1:
                   			
                   			CH_BUFF[CH_no].append(copy.copy(CH_ID))
                   			CH_BUFF[CH_no].append(copy.copy(j))
                   		else:
                   			CH_BUFF[CH_no].append(copy.copy(j))
                   		cl_size[CH_no]=cl_size[CH_no]+1
                   		CHIdx_buff[j]=CH_ID
                   		break
          	     
                   
                  
            #CH_BUFF.append(copy.copy(ch_buff))
        print('CHID_buff :'+str(CHID_buff))
        print('CHIdx_buff:'+str(CHIdx_buff))
        print('CH_BUFF:'+str(CH_BUFF))
        print('whole_buff:'+str(whole_buff))
                
        segAmount = len(CH_BUFF)
        for i in range(segAmount):
            num = len(CH_BUFF[i])
            if num!=0:
            	adjacency_matrix = numpy.zeros((num+1, num+1))
            	for j in range(num):
            		for k in range(j+1,num+1):
                		adjacency_matrix[j][k] = 1
                		adjacency_matrix[k][j] = 1

            	(whole_buff_update) = ClusterComputing(
                opt,
                whole_buff,
                pos_matrix,
                adjacency_matrix,
                0,
                CH_BUFF[i],
                )
            	whole_buff = whole_buff_update
            else:
            	ID=CHID_buff[i]
            	x=pos_matrix[ID][0]
            	y=pos_matrix[ID][1]
            	
            	
            	dist_temp=numpy.sqrt(x*x+y*y)
            	
            	print(dist_temp)
            	ene=0
            	if dist_temp<87:
            		ene=bitsnum*E+bitsnum*epson_fs*dist_temp*dist_temp
            	else:
            		ene=bitsnum*E+bitsnum*epson_mp*dist_temp*dist_temp*dist_temp*dist_temp
            	
            	print('ene is:'+str(ene))
            	print(whole_buff[ID])
            	whole_buff[ID]=whole_buff[ID]-ene
            	if whole_buff[ID]<0:
            		whole_buff[ID]=-1	

        # #print('whole_buff  is:' + str(whole_buff_update))

        count1 = 0
        for i in range(aus_num * (max_graph_size + 1)):
            if whole_buff[i] < 0:
                count1 = count1 + 1

        if count1 >= 7 and FND_flag == False:
            FND_num = count
            FND_flag = True

            # break

        if count1 >= 21 and HND_flag == False:
            HND_num = count
            HND_flag = True

            # break

        if count1 == 42 and LND_flag == False:
            LND_num = count
            break
            run_flag=False

        if count == 1500:
            DeadRatio_1500R = True
            DR_1500 = count1 / 36

    # ClusterComputing(1)

    t_stop = time.time()

    t_lap = t_stop - t_start
    print('processing time is:' + str(t_lap))
    print('FND_NUM is :' + str(FND_num))
    print('HND_NUM is :' + str(HND_num))
    print('LND_NUM is :' + str(LND_num))
    #print('LND_NUM1 is :' + str(LND_num_1))
    ##print 'seg_count is:' + str(seg_count)

    print('average time per round is:'+str(t_lap/LND_num))

    print('dead ratio at #1500 rounds' + str(DR_1500))


    # #print obsolete_idx

SpeedSlideWindow(2)
