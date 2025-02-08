#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy
import math
import copy
import random
import time
import sys

pi=3.1415926
d0=87.7085
E = 50 * 1e-09
epson_fs = 1e-12 * 10 
epson_mp = 0.0013 * 1e-12
bitsnum = 4000

def setSysVar():
	global pi,d0,E, epson_fs, epson_mp, bitsnum
	pi=3.1415926
	x=input('d0 is :')
	y=input('E is :')
	z=input('epson_fs is :')
	u=input('epson_mp is :')
	v=input('bitsnum is :')
	d0=float(x)
	E=float(y)
	epson_fs=float(z)
	epson_mp=float(u)
	bitsnum=int(v)
	print(d0)
	print(E)
	print(epson_fs)
	print(epson_mp)
	print(bitsnum)

#function to assign 2D position to all the sensor nodes within the sector i
def assignSensorPos(i, segAmount, ang,addon):
#assign radius value
    radius = 50   
   #assign sensor node amount per sector
    max_graph_size = 6+addon
    #assign accumulated angular value per sector
    ang_seg = 60 
    #create an array to save the 2D positions
    pos_matrix = numpy.zeros((2, max_graph_size)) 
    #assign each sensor node a random 2D position and save it 
    for i in range(max_graph_size):

        r_radius = random.random()
        temp_radius = r_radius * radius

        r_angular = random.random()
        temp_angular = (r_angular * ang_seg + ang - ang_seg) * 2 * pi \
            / 360

        x = temp_radius * math.cos(temp_angular)
        y = temp_radius * math.sin(temp_angular)

        pos_matrix[0, i] = x
        pos_matrix[1, i] = y

    return pos_matrix


def SpeedSlideWindow(opt1,opt2,sectIncrAmount):
    global pi,d0,E, epson_fs, epson_mp, bitsnum
    aus_num = 6
    whole_buff = []
    round_buff = []
    CHIdx_buff = []
    
    addon=sectIncrAmount

    #d0 = 87.7085
    lambdap = 20

    t_start = time.time()

    #E = 50 * 1e-09
    #epson_fs = 1e-12 * 10 
    #epson_mp = 0.0013 * 1e-12
    #bitsnum = 4000

    
    for i in range(aus_num):
        ene_buff = numpy.ones(aus_num+addon)
        round_b = numpy.zeros(aus_num+addon)
        id_buff = numpy.zeros(aus_num+addon)
        round_buff.append(copy.copy(round_b))
        for j in range(aus_num+addon):
            id_buff[j] = 37
        CHIdx_buff.append(copy.copy(id_buff))
        if i>0:
        	for j in range(aus_num+addon):
            		if j % 2 == 0:
                		ene_buff[j] = 0.5
            		else:
                		ene_buff[j] = 1
        whole_buff.append(copy.copy(ene_buff))

    max_graph_size = 6

    theta_bound = [
        60,
        120,
        180,
        240,
        300,
        360,
        ]

    segAmount = 6

    pos_matrix_buff = []
    f = open("posdata.txt", "a+")
    for i in range(segAmount):
        ang = theta_bound[i]
        pos_matrix = assignSensorPos(i, segAmount, ang,addon)
        name='posdata_folder/posdata'+str(i)+'.txt'
        numpy.savetxt(name,pos_matrix)
        pos_matrix_buff.append(copy.copy(pos_matrix))

    run_flag = True
    
    

    count = 0

    P = 0.2
 
    round_limit = int(1 / P)
    alpha = 1
    m = 0.5
    E_in = 0.5

    P_nrm = P / (1 + alpha * m)
    P_adv = P * (1 + alpha) / (1 + alpha * m)

    #b = random.random()
    b=opt1
    #u = random.random()
    u=opt2
    v = 1 - u

    P_bnrm = b * P_nrm
    P_badv = b * P_adv

    FND_flag = False
    HND_flag = False
    LND_flag = False
    DeadRatio_1500R = False
    DR_1500 = 0

    FND_num = 0
    HND_num = 0
    LND_num = 0

    timeing_flag = True
    
    ch_formed_flag=True
    ch_formed_count=0
    ch_notformed_flag=True
    ch_notformed_count=0

    while run_flag == True and count==0:
        count = count + 1

           # ##select the cluster head
           # find the d_max

        d_max = 0
        d_sum = 0
        gcount = 0
        e_sum = 0

        # #print('round_buff is:'+str(round_buff))

        #print 'whole buff is:' + str(whole_buff)
        p_time_start = time.time()
        alive_node_count=0
        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num+addon):
                if round_bag[j] <= round_limit and whole_buff[i][j] > 0:

                    # #print('first loop in')

                    # #print('prev is:'+str(g_count))

                    gcount = gcount + 1

                    # #print('next is:'+str(gcount))

                    x_temp = pos_matrix_buff[i][0][j]
                    y_temp = pos_matrix_buff[i][1][j]
                    d_temp = math.sqrt(x_temp * x_temp + y_temp
                            * y_temp)
                    d_sum = d_sum + d_temp
                    e_sum = e_sum + whole_buff[i][j]
                    if d_temp > d_max:
                        d_max = d_temp
                if whole_buff[i][j]>0:
                	alive_node_count=alive_node_count+1

        if gcount == 0:
            #print 'loop num is:' + str(count)
            #print 'limit count is:' + str(gcount)
            break

        d_ave = d_sum / gcount

        d_bs = d_ave

        e_ave = e_sum / gcount
        CHID_buff = []
        
        T_n=P/(1-P*(count%(1/P)))
        
        b=opt1*(alive_node_count/198)
        
        P_bnrm = b * P_nrm
        P_badv = b * P_adv
        
        ch_formed_flag=True
        ch_notformed_flag=True
        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num+addon):
                if round_bag[j] <= round_limit and whole_buff[i][j] > 0:

                   # #print('loop1 in')

                    
                    if j % 2 == 0:
                        E_start = 0.5
                    else:
                        E_start = 1

                    E_cur = whole_buff[i][j]

                    x_temp = pos_matrix_buff[i][0][j]
                    y_temp = pos_matrix_buff[i][1][j]
                    d_cur = math.sqrt(x_temp * x_temp + y_temp * y_temp)

                    T_n_mnrm = P_bnrm * (u * (E_cur / E_start)
                            + v*(d_cur / d_max) ) / (1 - P_nrm
                            * (count % (1 / P_nrm)))

                    T_n_madv = P_badv * (u * (E_cur / E_start)
                            + v*(d_cur / d_max)) / (1 - P_adv
                            * (count % (1 / P_adv)))

                    T_n_nrm1 = T_n_mnrm * (E_cur / E_start)
                    T_n_adv1 = T_n_madv * (E_cur / E_start)

                    T_n_nrm2 = T_n_nrm1 * e_ave
                    T_n_adv2 = T_n_adv1 * e_ave

                    T_n_fnrm = T_n_nrm2 * (1 / d_ave + 1 / d_cur)
                    T_n_fadv = T_n_adv2 * (1 / d_ave + 1 / d_cur)

                    threshold_nrm = random.random()
                    threshold_adv = random.random()

                    #print('mnrm is:' + str(T_n_mnrm))
                    #print('madv is:' + str(T_n_madv))
                    #print('T_n is:'+str(T_n))
                    
                    #T_n_mnrm=0
                    #T_n_madv=0

                    if j % 2 == 0:
                        if threshold_nrm < T_n:
                            CHIdx_buff[i][j] = 1
                            round_buff[i][j] = 0
                            ID = i * (aus_num+addon) + j
                            CHID_buff.append(copy.copy(ID))
                            round_bag[j] = round_bag[j] + 1
                            if ch_formed_flag==True:
                            	ch_formed_count=ch_formed_count+1
                            	ch_formed_flag=False
                            	ch_notformed_flag=False
                        
                        
                    else:
                        if threshold_adv < T_n:
                            #CHIdx_buff[i][j] = 1
                            round_buff[i][j] = 0
                            ID = i * (aus_num+addon) + j
                            CHID_buff.append(copy.copy(ID))
                            round_bag[j] = round_bag[j] + 1
                            if ch_formed_flag==True:
                            	ch_formed_count=ch_formed_count+1
                            	ch_formed_flag=False
                            	ch_notformed_flag=False
                      
                        
    # #form the cluster
        if ch_notformed_flag==True:
        	ch_notformed_count=ch_notformed_count+1

        #print('ch list is:' + str(CHID_buff))
        #for each node within each sector excluding the cluster heads, assign it to its selected cluster head that is the least distance away from the node
        for j in range(aus_num):
            for k in range(max_graph_size+addon):
                if whole_buff[j][k] > 0:

                    # #print('loop2 in')

                    x = pos_matrix_buff[j][0][k]
                    y = pos_matrix_buff[j][1][k]

                    dist = 20000
                    CH_ID = 0
                    flag=True
                    #for this node indexed as j,k, for each cluster head within the cluster head array, find the one with the least distance from the node j,k
                    for p in range(len(CHID_buff)):
                        ID = CHID_buff[p]
                        node_id = ID % (max_graph_size+addon)
                        seg_id = (ID - node_id) / (max_graph_size+addon)
                        node_id = int(node_id)
                        seg_id = int(seg_id)

                        x_ch = pos_matrix_buff[seg_id][0][node_id]
                        y_ch = pos_matrix_buff[seg_id][1][node_id]
                        dist_temp = math.sqrt((x - x_ch) * (x - x_ch)
                                + (y - y_ch) * (y - y_ch))

                        if dist_temp < dist:
                            dist = dist_temp
                            CH_ID = ID
                            flag=False

                    #print 'dist is:' + str(dist)
                    if flag==False:
                    	CHIdx_buff[j][k] = int(CH_ID)
                    else:
                    	CHIdx_buff[j][k]=37
                else:
                    CHIdx_buff[j][k]=37

           # ##deduct the energy buff

        print('cluster head array is:' + str(CHID_buff))
        print('each node corresponding its cluster head ID:'+str(CHIdx_buff))
        #print('whole_buff is:'+str(whole_buff))
        
        flag_CHID=True
        flag_whole=True
   #flag_CHID indicates no cluster head assigned while True, flag_whole indicates network die-out when False
        for i in range(aus_num):
        	for j in range(aus_num+addon):
        		if CHIdx_buff[i][j] !=37:
        			flag_CHID=False
        			break
        
        for i in range(aus_num):
        	for j in range(aus_num+addon):
        		if whole_buff[i][j] >0:
        			flag_whole=False
        			break
        
#network dies out and no node assigned/selected cluster head, assign value to LND
        if flag_CHID == True \
            and flag_whole == True:
   
            if LND_flag==False:
            	LND_num = count

            	LND_flag = True
            	run_flag=False
            	break
 	#network dies out flase and no cluster head assigned/selected true, directiont transmission between the node and the sink..
        elif flag_CHID == True \
            and flag_whole == False:
   
            for i in range(aus_num):
            	for j in range(aus_num+addon):
            		if whole_buff[i][j]>0:
            			x=pos_matrix_buff[i][0][j]
            			y=pos_matrix_buff[i][1][j]
            			dist=math.sqrt(x*x+y*y)
            			if dist<87.6:
            				ene=bitsnum*E+epson_fs*bitsnum*dist*dist
            			else:
            				ene=bitsnum*E+epson_mp*bitsnum*dist*dist*dist*dist
            			whole_buff[i][j]=whole_buff[i][j]-ene
            			if whole_buff[i][j]<0:
            				whole_buff[i][j]=-1
            
           #cluster head selected/assigned 
        elif flag_CHID == False:
            for i in range(aus_num):
                for j in range(aus_num+addon):
                    id_ch = int(CHIdx_buff[i][j])
                    if id_ch != 37:
                        #print id_ch
                        CH_ID = id_ch
                        if CH_ID != i * (aus_num+addon) + j:
                            p_random = random.random()
                            if p_random <= 1:

                                CH_ID_y = CH_ID % (aus_num+addon)
                                CH_ID_x = (CH_ID - CH_ID_y) / (aus_num+addon)

                                CH_ID_y = int(CH_ID_y)
                                CH_ID_x = int(CH_ID_x)

                                x = pos_matrix_buff[i][0][j]
                                y = pos_matrix_buff[i][1][j]

                                x_ch = \
                                    pos_matrix_buff[CH_ID_x][0][CH_ID_y]
                                y_ch = \
                                    pos_matrix_buff[CH_ID_x][1][CH_ID_y]

                                dist = math.sqrt((x - x_ch) * (x
                                        - x_ch) + (y - y_ch) * (y
                                        - y_ch))

                                if dist < d0:
                                    ene = E * bitsnum + epson_fs \
    * bitsnum * dist * dist
                                else:
                                    ene = E * bitsnum + epson_mp \
    * bitsnum * dist * dist * dist * dist

                                whole_buff[i][j] = whole_buff[i][j] \
                                    - ene

                                if whole_buff[i][j] < 0:
                                    whole_buff[i][j] = -1

                                dist_ch = math.sqrt(x_ch * x_ch + y_ch
                                        * y_ch)
                                if dist_ch < d0:
                                    ene = E * bitsnum * 2 + epson_fs \
    * bitsnum * dist_ch * dist_ch
                                else:
                                    ene = E * bitsnum * 2 + epson_mp \
    * bitsnum * dist_ch * dist_ch * dist_ch * dist_ch

                                whole_buff[CH_ID_x][CH_ID_y] = \
                                    whole_buff[CH_ID_x][CH_ID_y] - ene

                                if whole_buff[CH_ID_x][CH_ID_y] < 0:
                                    whole_buff[CH_ID_x][CH_ID_y] = -1
                        else:
                            CH_ID_y = CH_ID % (aus_num+addon)
                            CH_ID_x = (CH_ID - CH_ID_y) / (aus_num+addon)

                            CH_ID_y = int(CH_ID_y)
                            CH_ID_x = int(CH_ID_x)

                            x_ch = pos_matrix_buff[CH_ID_x][0][CH_ID_y]
                            y_ch = pos_matrix_buff[CH_ID_x][1][CH_ID_y]

                            dist = math.sqrt(x_ch * x_ch + y_ch * y_ch)

                            if dist < d0:
                                ene = E * bitsnum + epson_fs * bitsnum \
                                    * dist * dist
                            else:
                                ene = E * bitsnum + epson_mp * bitsnum \
                                    * dist * dist * dist * dist

                            whole_buff[CH_ID_x][CH_ID_y] = \
                                whole_buff[CH_ID_x][CH_ID_y] - ene
                            if whole_buff[CH_ID_x][CH_ID_y] < 0:
                                whole_buff[CH_ID_x][CH_ID_y] = -1

        count1 = 0
        for i in range(aus_num):
            for j in range(aus_num+addon):
                if whole_buff[i][j] < 0:
                    count1 = count1 + 1
        #print(whole_buff)
        #count the first node die round amount
        if count1 >= 1 and FND_flag == False:
            FND_num = count
            FND_flag = True

            # break
#count the half node die round amount
        if count1 >= 99 and HND_flag == False:
            HND_num = count
            HND_flag = True

            # break
#count the whole network die out round amount
        if count1 >= 198 and LND_flag == False:
            LND_num = count
            LND_flag=True
            break

        if count == 1500:
            DeadRatio_1500R = True
            DR_1500 = count1 / 198
            

        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num+addon):
                if round_bag[j] > round_limit:
                    round_bag[j] = 0

        for i in range(aus_num):
            for j in range(aus_num+addon):
                CHIdx_buff[i][j] = 37

        #run_flag=False

    #print('FND_NUM is :' + str(FND_num))
    #print('HND_NUM is :' + str(HND_num))
    #print('LND_NUM is :' + str(LND_num))
    #print('cluster head found count is:'+str(ch_formed_count))
    #print('cluster head not found count is:'+str(ch_notformed_count))
    #print('dead ratio at #1500 rounds' + str(DR_1500))
    #t_stop = time.time()
    #t_lap = t_stop - t_start

    #print('average time per round is:'+str(t_lap/LND_num))

    #print('processing time is:' + str(t_lap))


#SpeedSlideWindow(1,0.08)
