 
def assignSensorPos(i,segAmount,ang):
    max_graph_size=6
    for i in range(max_graph_size):
        r_radius = random.random()
        temp_radius = r_radius * radius

        r_angular = random.random()
        temp_angular = (r_angular * ang_seg + ang - ang_seg) * 2 * pi/ 360

        x = temp_radius * math.cos(temp_angular)
        y = temp_radius * math.sin(temp_angular)

        pos_matrix[0, i] = x
        pos_matrix[1, i] = y
     
    return pos_matrix

def SpeedSlideWindow():
    aus_num = 6
    whole_buff = []
    round_buff=[]
    CHIdx_buff=[]
    
    d0 = 87.7085
    lambdap = 20

    E = 50 * 1e-09
    epson_fs = 1e-09 * 10 * 0.001
    epson_mp = 0.0013 * 0.000000000001 * 0.001
    bitsnum = 4000

    
    for i in range(aus_num):
    	ene_buff = numpy.ones(aus_num)
    	round_b=numpy.zeros(aus_num)
    	for j in range(aus_num):
            if j%2 ==0:
               ene_buff[j] = 0.5
        whole_buff.append(copy.copy(ene_buff))
        round_buff.append(copy.copy(round_b))
        CHIdx_buff.append(copy.copy(round_b))
        
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
        pos_matrix = assignSensorPos(i,segAmount,ang)
        pos_matrix_buff.append(copy.copy(pos_matrix))
        
    run_flag=True
   
    count=0
    
    P=0.2
    round_limit=int(1/P)
    alpha=1
    m=0.5
    E_in=0.5
    
    P_nrm=P/(1+alpha*m)
    P_adv=P*(1+alpha)/(1+alpha*m)
    
    b=random.random()
    u=random.random()
    v=1-u
    
    P_bnrm=b*P_nrm
    P_badv=b*P_adv
    
    
    FND_flag = False
    HND_flag = False
    DeadRatio_1500R=False
    DR_1500=0
    
    while run_flag==True:
    	   count=count+1
    	   ###select the cluster head
    	   #find the d_max
    	   d_max=0
    	   d_sum=0
    	   limit_count=0
    	   e_sum=0
    	   for i in range(aus_num):
    	   	round_bag= round_buff[i]
    	   	for j in range(aus_num):
    	   		if round_bag[j]>=round_limit and whole_buff[i][j]>0:
    	   			limit_count=limit_count+1
    	   		        x_temp=pos_matrix_buff[i][j][0]
    	   		        y_temp=pos_matrix_buff[i][j][1]
    	   		        d_temp=math.sqrt(x_temp^2+y_temp^2)
    	   		        d_sum=d_sum+d_temp
    	   		        e_sum=e_sum+whole_buff[i][j]
    	   		        if d_temp>d_max:
    	   		           d_max=d_temp
  	   d_ave=d_sum/limit_count  	   		      
  	   
  	   d_bs=d_ave     
  	   
  	   e_ave=e_sum/limit_count
    	   
    	   CHID_buf=[]           
    	   for i in range(aus_num):
    	   	round_bag= round_buff[i]
    	   	for j in range(aus_num):
    	   		if round_bag[j]>=round_limit and whole_buff[i][j]>0:	
    	   			if j%2==0:
    	   				E_start=0.5
    	   			else:
    	   				E_start=1
    	   			
    	   			E_cur=whole_buff[i][j]
    	   			
    	   			x_temp=pos_matrix_buff[i][j][0]
    	   		        y_temp=pos_matrix_buff[i][j][1]
    	   		        d_cur=math.sqrt(x_temp^2+y_temp^2)
    	   		        		        
    	   		        T_n_mnrm=(P_bnrm*u*(u*(E_cur/E_start)+(d_cur/d_max)+(1/d_bs)))/(1-P_nrm*(count%(1/P_nrm)))
    	   		        
    	   		        T_n_madv=(P_badv*u*(u*(E_cur/E_start)+(d_cur/d_max)+(1/d_bs)))/(1-P_adv*(count%(1/P_nrm)))
    	   		        
    	   		        T_n_nrm1=T_n_mnrm*(E_cur/E_start)
    	   		        T_n_adv1=T_n_madv*(E_cur/E_start)
    	   
    	   			T_n_nrm2=T_n_nrm1*e_ave
    	   			T_n_adv2=T_n_adv1*e_ave
    	   			
    	   			T_n_fnrm=T_n_nrm2*(1/d_ave+1/d_cur)
    	   			T_n_fadv=T_n_adv2*(1/d_ave+1/d_cur)
    	   			
    	   			threshold_nrm=random.random()
    	   			threshold_adv=random.random()
    	   			
    	   			if j%2==0:
    	   		           if threshold_nrm<T_n_fnrm:
    	   		           	CHIdx_buff[i][j]=1
    	   		           	round_buff[i][j]=round_buff[i][j]+1
    	   		           else:
    	   		           	CHIdx_buff[i][j]=0
    	   		        else:
    	   		            if threshold_adv<T_n_fadv:
    	   		            	CHIdx_buff[i][j]=1
    	   		            	round_buff[i][j]=round_buff[i][j]+1
    	   		            else:
    	   		            	CHIdx_buff[i][j]=0
				
				ID=i*aus_num+j
				CHID_buff.append(copy.copy(ID))
	##form the cluster
	
	   for j in range(aus_num):
	       for k in range(aus_num):
	           if whole_buff[j][k]>0:
		   	x=pos_matrix_buff[j][k][0]
		   	y=pos_matrix_buff[j][k][1]
		   
		   	dist=20000
		   	CH_ID=0
		   	for p in range(len(CHID_buff)):
		   		ID=CHID_buff[i]
	       		seg_id=ID%6
	       		node_id=ID-seg_id*6
	       
	       		x_ch=pos_matrix_buff[seg_id][node_id][0]
	       		y_ch=pos_matrix_buff[seg_id][node_id][1]

				dist_temp=math.sqrt((x-x_ch)^2+(y-y_ch)^2)
			
				if dist_temp<dist:
					dist=dist_temp
					CH_ID=ID
		      else:
		      		CH_ID=37
		    
		 CHIdx_buff[j][k]=CH_ID		   		   		            	
    	   			
    	   ###deduct the energy buff
    	   
    	   for i in range(aus_num):
    	   	for j in range(aus_num):
    	   		if CHIdx_buff[j][k] !=37:
    	   		   CH_ID = CHIdx_buff[j][k]
    	   		   if CH_ID != i*aus_num+j:
    	   		        p_random=random.random()
    	   		        if p_random<0.5:
    	   		   		CH_ID_x=CH_ID%aus_num
    	   		   		CH_ID_y=CH_ID-aus_num*CH_ID_x
    	   		   	
    	   		   		x=pos_matrix_buff[i][j][0]
    	   		   		y=pos_matrix_buff[i][j][1]
    	   		   		
    	   		   		x_ch=pos_matrix_buff[CH_ID_x][CH_ID_y][0]
    	   		   		y_ch=pos_matrix_buff[CH_ID_x][CH_ID_y][1]
    	   		   	
    	   		   		dist=math.sqrt((x-x_ch)^2+(y-y_ch)^2)
    	   		   	
    	   		   		if dist<d0:
    	   		   	   		ene=E*bitsnum+epson_fs*bitsnum*dist^2
    	   		   		else:
    	   		   	   		ene=E*bitsnum+epson_mf*bitsnum*dist^2*dist^2
    	   		   	   
    	   		   		whole_buff[i][j]=whole_buff[i][j]-ene
    	   		   		
    	   		   		if whole_buff[i][j]<0:
    	   		   			whole_buff[i][j]=-1
    	   		   	
    	   		   	
    	   		        	dist_ch=math.sqrt(x_ch^2+y_ch^2)
    	   		        	if dist_ch<d0:
    	   		        		ene=E*bitsnum*2+epson_fs*bitsnum*dist_ch^2
    	   		        	else:
    	   		        		ene=E*bitsnum*2+epson_mf*bitsnum*dist_ch^2*dist_ch^2		   
    	   		   	
    	   		   		whole_buff[CH_ID_x][CH_ID_y]=whole_buff[CH_ID_x][CH_ID_y]-ene
    	   		   		
    	   		   		if whole_buff[CH_ID_x][CH_ID_y]<0:
    	   		   			whole_buff[CH_ID_x][CH_ID_y]=-1
    	   		   			
    	   		   	
    	   
    	   ###decide the amount of alive nodes and dead nodes
    	   ###decide whether run_flag=Flase or break the loop
    	  if count1 >= 1 and FND_flag == False:
            FND_num = count
            FND_flag = True

            # break

         if count1 >= 18 and HND_flag == False:
            HND_num = count
            HND_flag = True

            # break

         if count1 == 36:
            LND_num = count
            break
        
         if count==1500:
            DeadRatio_1500R=True
            DR_1500=(36-count1)/36
            
             print('FND_NUM is :' + str(FND_num))
    print('HND_NUM is :' + str(HND_num))
    print('LND_NUM is :' + str(LND_num))
    print('dead ratio at #1500 rounds'+str(DR_1500))
    
SpeedSlideWindow()
