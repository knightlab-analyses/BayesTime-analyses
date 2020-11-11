#-------------------------------------------
# simulation setting for SFPCA
# N: number of subjects
# N_T_max: number of maximum time points
# N_T_mean: number of average time points
#-------------------------------------------
sim_data = function(N, N_T_max, N_T_mean){
	T_min=0; T_max=1; T_range=c(T_min,T_max)
	P=2
	nknots=c(1,2)
	params=list()
	orth=TRUE

	Q_true=c(2,2) # number of modes of variation (i.e. number of PCs)
	R_true=diag(rep(1,sum(Q_true)));R_true[1,4]=.75;R_true[4,1]=R_true[1,4] # correlation of FPCA scores
	D_true=diag(rep(1,sum(Q_true)));D_true[1,1]=2;D_true[3,3]=2; #D_true=0.1*D_true # standard deviations of FPCA scores
	SIGMA_ALPHA_true=D_true%*%R_true%*%D_true # Variance-covariance matrix of FPCA scores
	SIGMA_OMEGA_true=0.01*rep(1,P) # Error variances
	MU_true=c(c(1.8, .8, .4, .2, -1.2), 
			  c(0.125689927, 0.085802413, -0.107942359, -0.164469410, -0.009916578, -0.081660406))


	## set up FPC_true
	time=seq(T_range[1],T_range[2],(T_range[2]-T_range[1])/(N_T_max-1))
	time_sparse=list()
	for(i in 1:N){
		time_sparse[[i]]=list()
		for(p in 1:P){
			time_sparse[[i]][[p]]=time
		}
	}
	basis_stuff=basis_setup_sparse(time_sparse,nknots, orth=orth)
	knots=basis_stuff[[1]]
	time_sparse=basis_stuff[[2]]
	time_sparse_combined=basis_stuff[[3]]
	phi_t=basis_stuff[[4]]
	time_cont=basis_stuff[[5]]
	phi_t_cont=	basis_stuff[[6]]
	Phi_D=basis_stuff[[7]]
	Phi_D_T=basis_stuff[[8]]

	# FPCA loadings
	if(nknots[2]==1){
		THETA_true=cbind(c(0.40529575,0.82201395,0.39930063,-0.02187301,-0.01044030,rep(0,5)),
			c(0.06137363,0.10733627,-0.33973487,-0.81361591,-0.45532579,rep(0,5)),
			c(rep(0,5),0.3490362,0.7295541,0.5036934,0.2608049,0.1555581),
			c(rep(0,5),0.2280605,0.3930486,-0.2056171,-0.7682501,-0.4012662))
	}
	if(nknots[2]==2){
		THETA_true=cbind(c(0.40529575,0.82201395,0.39930063,-0.02187301,-0.01044030,rep(0,6)),
			c(0.06137363,0.10733627,-0.33973487,-0.81361591,-0.45532579,rep(0,6)),
			c(rep(0,5), 0.1734486, 0.2471184, 0.5260293, 0.5316424, 0.4836753, 0.2758460),
			c(rep(0,5), 0.79569747, 0.31964112, -0.16640001, -0.17366537, -0.04200097, -0.05197922))
	}

	FPC_true=array(0,dim=c(N_T_max,sum(Q_true)))
	for(t in 1:N_T_max){
		ind=0
		for(q in 1:length(Q_true)){
			FPC_true[t,(ind+1):(ind+Q_true[q])]=(Phi_D[[1]][[t]]%*%THETA_true)[q,(ind+1):(ind+Q_true[q])]
			ind=ind+Q_true[q]
		}	
	}

	## Set up param list

	params[[1]]=nknots
	params[[2]]=Q_true
	params[[3]]=R_true
	params[[4]]=D_true
	params[[5]]=SIGMA_OMEGA_true
	params[[6]]=THETA_true
	params[[7]]=MU_true
	params[[8]]=FPC_true

	return(list('params'=params, 'N'=N, 'Missing'=(1-N_T_mean/N_T_max)*100,
		        'T_range'=T_range, 'N_T_max'=N_T_max, 'N_T_mean'=N_T_mean,
		        'P'=P, 'orth'=orth))
}


#-------------------------------------------
# reformat simulation data
#-------------------------------------------
sim_data_reformat = function(Y_SPARSE, TIME_SPARSE, MU_SPARSE, F_SPARSE, OMEGA_SPARSE, 
							 ALPHA, params, nsim){
	J = 2
	N = length(Y_SPARSE[[nsim]])
	nknots = params[[1]][[J]] 
	Q = nknots + 4 # number of basis
	K = params[[2]][[J]] # number of PC
	sigma_eps = sqrt(params[[5]][[J]]) # SIGMA_OMEGA_true is error variance

	Q1 = params[[1]][[1]] + 4
	K1 = params[[2]][[1]]

	Theta = params[[6]][(1+Q1):(Q+Q1), (1+K1):(K+K1)] # THETA_true
	theta_mu = params[[7]][(1+Q1):(Q+Q1)] # MU_true
	ALPHA = t(ALPHA[[nsim]][(1+K1):(K+K1), ]) # ALPHA

	# basis
	phi_t = phi_t[[nsim]][[J]]
	phi_t_stacked=NULL
	for(i in 1:N){
		phi_t_stacked = rbind(phi_t_stacked, t(phi_t[[i]]))
	}


	## simulated data
	time = Y = MU = F = Omega = Alpha = list() 
	ids=rep(1:N,each=1)

	for (i in 1:N){
		time[[i]] = TIME_SPARSE[[nsim]][[i]][[J]]
		Y[[i]] = Y_SPARSE[[nsim]][[i]][[J]]
		MU[[i]] = MU_SPARSE[[nsim]][[i]][[J]]
		F[[i]] = F_SPARSE[[nsim]][[i]][[J]]
		Omega[[i]] = OMEGA_SPARSE[[nsim]][[i]][[J]]
	}

	## convert data to data frame
	nrows = length(unlist(time))
	nvisits = sapply(time, length)
	df = data.frame(matrix(rep(0, nrows * 3), nrow=nrows))
	colnames(df) = c('id', 'time', 'response')
	for (j in 1:nrows){
		df$id = rep(1:N, nvisits)
		df$time = unlist(time)
		df$response = unlist(Y)
	}

	return(list('df'= df, 'nknots'=nknots, 'Q'=Q, 'K'=K, 'Q1'=Q1, 'K1'=K1, 'J'=J, 'N'=N,
		        'phi_t'=phi_t, 'phi_t_stacked'=phi_t_stacked))

}
