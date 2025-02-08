#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy, math,copy, random,time
def assignSensorPos(i, segAmount, ang):
    radius=100
    max_graph_size = 6
    ang_seg=60
    pi=3.1415926
    pos_matrix=numpy.zeros((2,6))
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


def SpeedSlideWindow():
    aus_num = 6
    whole_buff = []
    round_buff = []
    CHIdx_buff = []

    d0 = 87.7085
    lambdap = 20
    
    t_start=time.time()
    
    E = 50 * 1e-09
    epson_fs = 1e-09 * 10 * 0.001
    epson_mp = 0.0013 * 0.000000000001 * 0.001
    bitsnum = 4000

    for i in range(aus_num):
        ene_buff = numpy.ones(aus_num)
        round_b = numpy.zeros(aus_num)
        id_buff=numpy.zeros(aus_num)
        round_buff.append(copy.copy(round_b))
        for j in range(aus_num):
            id_buff[j] = 37
        CHIdx_buff.append(copy.copy(id_buff))
        for j in range(aus_num):
            if j % 2 == 0:
                ene_buff[j] = 0.5
            else:
            	ene_buff[j]=1
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

    for i in range(segAmount):
        ang = theta_bound[i]
        pos_matrix = assignSensorPos(i, segAmount, ang)
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

    b = random.random()
    u = random.random()
    v = 1 - u

    P_bnrm = b * P_nrm
    P_badv = b * P_adv

    FND_flag = False
    HND_flag = False
    LND_flag=False
    DeadRatio_1500R = False
    DR_1500 = 0

    FND_num=0
    HND_num=0
    LND_num=0
    
    timeing_flag=True
    
    while run_flag == True:
        count = count + 1

           # ##select the cluster head
           # find the d_max

        d_max = 0
        d_sum = 0
        gcount = 0
        e_sum = 0
        #print('round_buff is:'+str(round_buff))
        print('whole buff is:'+str(whole_buff))
        p_time_start=time.time()
        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num):
                if round_bag[j] <= round_limit and whole_buff[i][j] > 0:
                    #print('first loop in')
                    
                    #print('prev is:'+str(g_count))
                    gcount=gcount+1
                    #print('next is:'+str(gcount))
                    x_temp = pos_matrix_buff[i][0][j]
                    y_temp = pos_matrix_buff[i][1][j]
                    d_temp = math.sqrt(x_temp*x_temp + y_temp*y_temp)
                    d_sum = d_sum + d_temp
                    e_sum = e_sum + whole_buff[i][j]
                    if d_temp > d_max:
                        d_max = d_temp
                   
        
        if gcount==0:
           print('loop num is:'+str(count))
           print('limit count is:'+str(gcount))
           break
        	               
        d_ave = d_sum / gcount

        d_bs = d_ave

        e_ave = e_sum / gcount
        CHID_buff = []
        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num):
                if round_bag[j] <= round_limit and whole_buff[i][j] > 0:
                   #print('loop1 in')
                   round_bag[j]=round_bag[j]+1
                   if j % 2 == 0:
                        E_start = 0.5
                   else:
                        E_start = 1

                   E_cur = whole_buff[i][j]

                   x_temp = pos_matrix_buff[i][0][j]
                   y_temp = pos_matrix_buff[i][1][j]
                   d_cur = math.sqrt(x_temp*x_temp + y_temp*y_temp)

                   T_n_mnrm = P_bnrm * u * (u * (E_cur / E_start)
                            + d_cur / d_max + 1 / d_bs) / (1 - P_nrm
                            * (count % (1 / P_nrm)))

                   T_n_madv = P_badv * u * (u * (E_cur / E_start)
                            + d_cur / d_max + 1 / d_bs) / (1 - P_adv
                            * (count % (1 / P_adv)))

                   T_n_nrm1 = T_n_mnrm * (E_cur / E_start)
                   T_n_adv1 = T_n_madv * (E_cur / E_start)

                   T_n_nrm2 = T_n_nrm1 * e_ave
                   T_n_adv2 = T_n_adv1 * e_ave

                   T_n_fnrm = T_n_nrm2 * (1 / d_ave + 1 / d_cur)
                   T_n_fadv = T_n_adv2 * (1 / d_ave + 1 / d_cur)

                   threshold_nrm = random.random()
                   threshold_adv = random.random()
                   
                   print('fnrm is:'+str(T_n_fnrm))
                   print('fadv is:'+str(T_n_fadv))

                   if j % 2 == 0:
                      if threshold_nrm < T_n_fnrm:
                            CHIdx_buff[i][j] = 1
                            round_buff[i][j] = 0
                            ID = i * aus_num + j
                            CHID_buff.append(copy.copy(ID))
                      else:
                            CHIdx_buff[i][j] = 0
                   else:
                       if threshold_adv < T_n_fadv:
                            CHIdx_buff[i][j] = 1
                            round_buff[i][j] = 0
                            ID = i * aus_num + j
                            CHID_buff.append(copy.copy(ID))
                       else:
                            CHIdx_buff[i][j] = 0

                   

    # #form the cluster
        print('ch list is:'+str(CHID_buff))
        for j in range(aus_num):
            for k in range(aus_num):
                if whole_buff[j][k] > 0:
                    #print('loop2 in')
                    x = pos_matrix_buff[j][0][k]
                    y = pos_matrix_buff[j][1][k]

                    dist = 20000
                    CH_ID = 0
                    for p in range(len(CHID_buff)):
                        ID = CHID_buff[p]
                        node_id = ID % 6
                        seg_id = (ID - node_id)/6
                        node_id=int(node_id)
                        seg_id=int(seg_id)

                        x_ch = pos_matrix_buff[seg_id][0][node_id]
                        y_ch = pos_matrix_buff[seg_id][1][node_id]
                        dist_temp = math.sqrt((x - x_ch)*(x-x_ch) + (y - y_ch)*(y-y_ch))
                        
                        if dist_temp < dist:
                           dist = dist_temp
                           CH_ID = ID
                       
                    print('dist is:'+str(dist))       
                    CHIdx_buff[j][k] = int(CH_ID)
	
           # ##deduct the energy buff
        print('cluster head array is:'+str(CHIdx_buff))
        if numpy.all(CHIdx_buff==37.)==True and LND_flag==False and numpy.all(whole_buff<0)==True:
            LND_num = count
            
            LND_flag=True
            break
        	 
        elif numpy.all(CHIdx_buff==37.)==True and numpy.all(whole_buff<0)==False:
        	count=count-1
        elif numpy.all(CHIdx_buff==37.)==False:
            for i in range(aus_num):
            	for j in range(aus_num):
            		id_ch=int(CHIdx_buff[i][j])
            		if id_ch !=37:
            	   		print(id_ch)
            	   		CH_ID = id_ch
            	   		if CH_ID != (i * aus_num + j):
                      			p_random = random.random()
                      			if p_random <= 1:
                            #print('loop3')
                            			CH_ID_y = CH_ID % aus_num
                            			CH_ID_x = (CH_ID - CH_ID_y)/aus_num
                            
                            			CH_ID_y=int(CH_ID_y)
                            			CH_ID_x=int(CH_ID_x)

                            			x = pos_matrix_buff[i][0][j]
                            			y = pos_matrix_buff[i][1][j]

                            			x_ch = pos_matrix_buff[CH_ID_x][0][CH_ID_y]
                            			y_ch = pos_matrix_buff[CH_ID_x][1][CH_ID_y]

                            			dist = math.sqrt((x - x_ch)*(x-x_ch) + (y - y_ch)*(y-y_ch))

                            			if dist < d0:
                                			ene = E * bitsnum + epson_fs * bitsnum* dist *dist
                            			else:
                            				ene = E * bitsnum + epson_mp * bitsnum * dist * dist *dist*dist

                            			whole_buff[i][j] = whole_buff[i][j] - ene

                            			if whole_buff[i][j] < 0:
                                			whole_buff[i][j] = -1

                            			dist_ch = math.sqrt(x_ch*x_ch + y_ch*y_ch)
                            			if dist_ch < d0:
                                			ene = E * bitsnum * 2 + epson_fs* bitsnum * dist_ch *dist_ch
                            			else:
                                			ene = E * bitsnum * 2 + epson_mp* bitsnum * dist_ch *dist_ch * dist_ch *dist_ch

                            			whole_buff[CH_ID_x][CH_ID_y] = whole_buff[CH_ID_x][CH_ID_y] - ene

                            			if whole_buff[CH_ID_x][CH_ID_y] < 0:
                                			whole_buff[CH_ID_x][CH_ID_y] = -1
                                else:
                                	CH_ID_y = CH_ID % aus_num
                            		CH_ID_x = (CH_ID - CH_ID_y)/aus_num
                            
                            		CH_ID_y=int(CH_ID_y)
                            		CH_ID_x=int(CH_ID_x)

                            			
                            		x_ch = pos_matrix_buff[CH_ID_x][0][CH_ID_y]
                            		y_ch = pos_matrix_buff[CH_ID_x][1][CH_ID_y]

                            		dist = math.sqrt(x_ch*x_ch + y_ch*y_ch)

                            		if dist < d0:
                                		ene = E * bitsnum + epson_fs * bitsnum* dist *dist
                            		else:
                            			ene = E * bitsnum + epson_mp * bitsnum * dist * dist *dist*dist

                         		whole_buff[CH_ID_x][CH_ID_y] = whole_buff[CH_ID_x][CH_ID_y] - ene
                         		if whole_buff[CH_ID_x][CH_ID_y] < 0:
                                			whole_buff[CH_ID_x][CH_ID_y] = -1
	
        count1=0                       
        for i in range(aus_num):
        	for j in range(aus_num):
        		if whole_buff[i][j]<0:
        			count1=count1+1

        if count1 >= 1 and FND_flag == False:
            FND_num = count
            FND_flag = True

            # break

        if count1 >= 18 and HND_flag == False:
            HND_num = count
            HND_flag = True

            # break

        if count1 == 36 and LND_flag==False:
            LND_num = count
            break

        if count == 1500:
            DeadRatio_1500R = True
            DR_1500 = count1 / 36
        
        for i in range(aus_num):
            round_bag = round_buff[i]
            for j in range(aus_num):
            	if round_bag[j]>round_limit:
            		round_bag[j]=0 
        
        for i in range(aus_num):
        	for j in range(aus_num):
        		CHIdx_buff[i][j]=37
        #run_flag=False

    print('FND_NUM is :' + str(FND_num))
    print('HND_NUM is :' + str(HND_num))
    print('LND_NUM is :' + str(LND_num))
    print('dead ratio at #1500 rounds' + str(DR_1500))
    t_stop = time.time()
    t_lap = t_stop - t_start
    #print('average time per round is:'+str(t_lap/LND_num))
    print('processing time is:'+str(t_lap))

SpeedSlideWindow()
