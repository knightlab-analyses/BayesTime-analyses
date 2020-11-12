simulate_data_sparse=function(T_range, N_T_max, N_T_mean=N_T_max, P, N, params, orth=FALSE, nseeds){
	set.seed(nseeds)
		
	library(MASS)

	nknots=params[[1]]
	ids=rep(1:N,each=1)
	possible_times=seq(T_range[1],T_range[2],(T_range[2]-T_range[1])/(N_T_max-1))
	K_true=nknots+4

	# select number and placing of possible_times points
	ind=list()
	for(p in 1:P){
		ind[[p]]=array(0,dim=c(N_T_max,N))
	}	
	if(N_T_mean<N_T_max){
		for(i in 1:N){
			for(p in 1:P){
				lambda=max(0,N_T_mean-3)
				tmp1=sample(N_T_max)[1:(rpois(1,lambda)+3)]
				ind[[p]][tmp1,i]=1
			}	
		}	
	}
	if(N_T_mean==N_T_max){
		for(p in 1:P){
			ind[[p]]=array(1,dim=c(N_T_max,N))
		}
	}
	time_sparse=list()
	for(i in 1:N){
		time_sparse[[i]]=list()
		for(p in 1:P){
			time_sparse[[i]][[p]]=possible_times[ind[[p]][,i]==1]
		}
	}
	
	# set up basis functions
	basis_stuff=basis_setup_sparse(time_sparse,nknots,plotit=TRUE,orth=orth)
	knots=basis_stuff[[1]]
	time_sparse=basis_stuff[[2]]
	phi_t=basis_stuff[[3]]
	time_cont=basis_stuff[[4]]
	phi_t_cont=	basis_stuff[[5]]
	#time_sparse_combined=basis_stuff[[6]]
	#Phi_D_combined=basis_stuff[[7]]
	#Phi_D_T_combined=basis_stuff[[8]]

	if(length(params)==8){
		Q_true=params[[2]]
		R_true=params[[3]]
		D_true=params[[4]]
		SIGMA_ALPHA_true=D_true%*%R_true%*%D_true
		SIGMA_OMEGA_true=params[[5]]
		THETA_true=params[[6]]
		MU_true=params[[7]]
		FPC_true=params[[8]]		
		ALPHA_true=array(0,dim=c(sum(Q_true),N))
		for(i in 1:N){
			ALPHA_true[,i]=mvrnorm(1,rep(0,sum(Q_true)),SIGMA_ALPHA_true)
		}
	}										

	Y_sparse=list()
	MU_sparse_true=list()
	F_sparse_true=list()
	Omega_sparse_true=list()
	for(i in 1:N){
		Y_sparse[[i]]=list()
		MU_sparse_true[[i]]=list()
		F_sparse_true[[i]]=list()
		Omega_sparse_true[[i]]=list()
		indk=1
		indq=1
		for(p in 1:P){
			MU_p=MU_true[indk:sum(K_true[1:p])]
			ALPHA_ip=ALPHA_true[,i][indq:sum(Q_true[1:p])]
			MU_sparse_true[[i]][[p]]=t(phi_t[[p]][[i]])%*%MU_p	
			THETA_p=THETA_true[indk:sum(K_true[1:p]),indq:sum(Q_true[1:p])]
			F_sparse_true[[i]][[p]]=t(phi_t[[p]][[i]])%*%THETA_p%*%ALPHA_ip
			Omega_sparse_true[[i]][[p]]=rnorm(length(time_sparse[[i]][[p]]),0,SIGMA_OMEGA_true[p]^.5)
			Y_sparse[[i]][[p]]=MU_sparse_true[[i]][[p]]+F_sparse_true[[i]][[p]]+Omega_sparse_true[[i]][[p]]
			indk=indk+K_true[p]
			indq=indq+Q_true[p]
		}
	}	
	
	simdata=list()
	simdata[[1]]=time_sparse
	simdata[[2]]=Y_sparse
	simdata[[3]]=MU_sparse_true
	simdata[[4]]=F_sparse_true
	simdata[[5]]=Omega_sparse_true
	simdata[[6]]=ALPHA_true
	simdata[[7]]=possible_times
	simdata[[8]]=params

	## add output of basis
	simdata[[9]]=basis_stuff[[3]] # phi_t

	return(simdata)

}
