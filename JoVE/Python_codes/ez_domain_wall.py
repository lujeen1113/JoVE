import numpy as np

#written by Nicholas Chancellor
#Based on the work in Nicholas Chancellor 2019 Quantum Sci. Technol. 4 045004 https://doi.org/10.1088/2058-9565/ab33c2
#feel free to resuse/modify but please attribute the source and cite the paper in any published work

# the function of this module is to provide an easy utility to convert integer problems with arbitrary two variable interactions into the domain wall encoding
# the main function is 'make_domain_wall_encoding', it takes three inputs, 'var_sizes' is a list of the sizes of all of the integer variables, 'penalties' is the list of of the penalties assigned to each integer variable based on it's value, it is a vector where the length has to equal the sum of 'var_sizes', and 'interactions' is an upper triangular array of interactions, if it is not upper triangular, all elements in the lower triangle are ignored

# make_domain_wall_encode returns four variables, Jcore and hcore, which are the Ising couplings (upper triangular) and fieds to enforce the domain wall constraints, and Jprob and hprob, which actually enoforce the problem

def make_domain_wall_encoding(var_sizes,penalties,interactions):
    # creates a domain wall encoding as described above	
    ##print(var_sizes)
    Ntot=sum(var_sizes)-len(var_sizes)
    Ntot=int(Ntot)
    ##print('dim is'+str(Ntot))
    J_core=np.zeros((Ntot,Ntot))
    h_core=np.zeros(Ntot) # make matrix of zeros
    startNum=0 # keeps track of where each variable core starts
    varStarts=np.zeros(len(var_sizes),dtype=int) # make matrix of zeros to keep track of where each variable starts
    for iCore in range(len(var_sizes)): # create core for each variable 
        qubits=list(range(startNum,startNum+var_sizes[iCore]-1)) # qubits involved in core
        ##print(qubits)
        (J,h)=make_var_encode_core_dw(qubits,Ntot) # make domain wall core
        J_core=J_core+J # add to total core coupling
        h_core=h_core+h # add to total core fields
        varStarts[iCore]=int(startNum) # track start of variable
        startNum=startNum+var_sizes[iCore]-1 # increment starting qubit number
    J_prob=np.zeros([Ntot,Ntot]) # make matrix of zeros
    h_prob=np.zeros(Ntot) # make matrix of zeros
    for iProb in range(len(var_sizes)): # outer loop for problem creation
        startNum_i=varStarts[iProb] # starting qubit for variable
        qubits_i=list(range(startNum_i,startNum_i+var_sizes[iProb]-1)) # qubits involved in core
        for iValue in range(var_sizes[iProb]): # loop over values for each variable
            penNumber=startNum_i+iProb+iValue # position in penalty vector
            if not penalties[int(penNumber)]==0: # if there is a penalty
                h=penalties[penNumber]*value_penalty_dw(qubits_i,iValue,Ntot) # encoding of penalty
                h_prob=h_prob+h # add to total problem
        for jProb in range(iProb+1,len(var_sizes)): # second index for qubit interactions
            startNum_j=varStarts[jProb] # starting qubit for variable
            qubits_j=list(range(startNum_j,startNum_j+var_sizes[jProb]-1)) # qubits involved in core
            for iValue in range(var_sizes[iProb]): # loop over values for first variable
                intNumber_i=startNum_i+iProb+iValue # i position in interaction matrix
                for jValue in range(var_sizes[jProb]): # loop over values for second variable
                    intNumber_j=startNum_j+jProb+jValue # j position in interaction matrix
                    if not interactions[int(intNumber_i),int(intNumber_j)]==0: # if there is an interaction
                        (J,h)=two_var_interaction_dw(qubits_i,qubits_j,iValue,jValue,Ntot) # interactions
                        J_prob=J_prob+interactions[intNumber_i,intNumber_j]*J # multiply by strength and add
                        h_prob=h_prob+interactions[intNumber_i,intNumber_j]*h # multiply by strength and add
    return J_core,h_core,J_prob,h_prob


def make_var_encode_core_dw(qubits,Ntot):
	# makes a domain wall encoding for a single variable over the qubits 'qubits'
	h=np.zeros(Ntot) # empty list for putting fields in
	J=np.zeros([Ntot,Ntot]) # make matrix of zeros
	if len(qubits)>1: # only if there is more than one qubit in the variable
		h[qubits[0]]=1 # encourage first variable to be one
		h[qubits[len(qubits)-1]]=-1 # encourages last variable to be zero
		J[qubits[0:(len(qubits)-1)],qubits[1:len(qubits)]]=-1 # ferromagnetic interactions
		J[qubits[1:len(qubits)],qubits[0:(len(qubits)-1)]]=-1 # ferromagnetic interactions
		J=np.triu(J) # upper triangular part
	return (J,h)


def value_penalty_dw(qubits,val,Ntot):
    # encodes a penalty on a single value of a domain wall encoded integer variable
    h=np.zeros(Ntot) # empty list for putting fields in
    if val==0:	# special case of minimum possible value
        #print(qubits)
        h[qubits[0]]=1 # renforces penalty on first qubit
    elif val==len(qubits):	# special case of maximum possible value
        h[qubits[len(qubits)-1]]=-1 # reinforces penalty on last qubit
    else: # case where no interactions with val2 are replaced by fields
        h[qubits[val-1]]=-1 # increased cost of having domain wall above this value
        h[qubits[val]]=1 # increased cost of having domain wall below this value
    h=h/2 # factor of two between Ising and Qubo
    return h



def two_var_interaction_dw(qubits1,qubits2,val1,val2,Ntot):
	# encodes an interaction penalizing the simultaneous case where variable 1 takes val1 and variable 2 takes val 2
	h=np.zeros(Ntot) # empty list for putting fields in
	J=np.zeros([Ntot,Ntot]) # make matrix of zeros
	if val1==0: # special case where some couplers are replaced by fields
		if val2==0:	# special case where some coupers are replaced by fields
			J[qubits1[0],qubits2[0]]=1 # anti-ferro coupling
			h[qubits1[0]]=1 # encourages qubit to be 1
			h[qubits2[0]]=1 # encourages qubit to be 1
		elif val2==len(qubits2):	# special case where some coupers are replaced by fields
			J[qubits1[0],qubits2[len(qubits2)-1]]=-1 # ferro coupling
			h[qubits1[0]]=1 # encourages qubit to be 1
			h[qubits2[len(qubits2)-1]]=-1 # encourages qubit to be 0
		else: # case where no interactions with val2 are replaced by fields
			J[qubits1[0],qubits2[val2-1]]=-1 # ferro coupling
			J[qubits1[0],qubits2[val2]]=1 # anti-ferro coupling
			h[qubits2[val2-1]]=-1 # encourages qubit to be 0
			h[qubits2[val2]]=1 # encourages qubit to be 1
	elif  val1==len(qubits1): # special case where some couplers are replaced by fields
		if val2==0:	# special case where some coupers are replaced by fields
			J[qubits1[len(qubits1)-1],qubits2[0]]=-1 # ferro coupling
			h[qubits1[len(qubits1)-1]]=-1 # encourages qubit to be 0
			h[qubits2[0]]=1 # encourages qubit to be 1
		elif val2==len(qubits2):	# special case where some coupers are replaced by fields
			J[qubits1[len(qubits1)-1],qubits2[len(qubits2)-1]]=1 # anti-ferro coupling
			h[qubits1[len(qubits1)-1]]=-1 # encourages qubit to be 0
			h[qubits2[len(qubits2)-1]]=-1 # encourages qubit to be 0
		else: # case where no interactions with val2 are replaced by fields
			J[qubits1[len(qubits1)-1],qubits2[val2-1]]=1 # anti-ferro coupling
			J[qubits1[len(qubits1)-1],qubits2[val2]]=-1 # ferromagnetic coupling
			h[qubits2[val2-1]]=-1 # encourages qubit to be 0
			h[qubits2[val2]]=1 # encourages qubit to be 1
	else: # val1 is not at either end point of the chain
		if val2==0:	# special case where some coupers are replaced by fields
			J[qubits1[val1],qubits2[0]]=1 # anti-ferro coupling
			J[qubits1[val1-1],qubits2[0]]=-1 # ferro coupling
			h[qubits1[val1]]=1 # encourages qubit to be 1
			h[qubits1[val1-1]]=-1 # encourages qubit to be 0
		elif val2==len(qubits2):	# special case where some coupers are replaced by fields
			J[qubits1[val1],qubits2[len(qubits2)-1]]=-1 # ferro coupling
			J[qubits1[val1-1],qubits2[len(qubits2)-1]]=1 # anti-ferro coupling
			h[qubits1[val1]]=1 # encourages qubit to be 1
			h[qubits1[val1-1]]=-1 # encourages qubit to be 0
		else: # case where no values are at the end of chains
			J[qubits1[val1],qubits2[val2]]=1 # anti-ferro coupling
			J[qubits1[val1-1],qubits2[val2-1]]=1 # anti-ferro coupling
			J[qubits1[val1],qubits2[val2-1]]=-1 # ferro coupling
			J[qubits1[val1-1],qubits2[val2]]=-1 # ferro coupling	
	J=J+J.T # make symetric
	J=np.triu(J) # make upper triangular
	J=J/4 # factor of two between Ising and Qubo
	h=h/4 # factor of two between Ising and Qubo
	return (J,h)


