basis_setup_sparse=function(time_sparse,nknots=rep(1,length(time_sparse[[1]])),
		orth=FALSE,delta=1/10000,plotit=TRUE){


	#orth=FALSE
	#delta=1/10000
	#plotit=TRUE

	library(splines)
	library(pracma)

	N=length(time_sparse)
	P=length(time_sparse[[1]])
	time_unique=list()
	for(p in 1:P){
		for(i in 1:N){
			time_sparse[[i]][[p]]=round(time_sparse[[i]][[p]]/delta)*delta
		}	
		time_unique[[p]]=time_sparse[[1]][[p]]
		for(i in 2:N){
			time_unique[[p]]=c(time_unique[[p]],time_sparse[[i]][[p]])
		}
		time_unique[[p]]=round(sort(unique(time_unique[[p]]))/delta)*delta
	}	
	time_unique_combined=time_unique[[1]]
	for(p in 2:P){
		time_unique_combined=c(time_unique_combined,time_unique[[p]])
	}
	time_unique_combined=sort(unique(time_unique_combined))
	
	T_min=min(time_unique_combined)
	T_max=max(time_unique_combined)
	time_cont=seq(T_min,T_max/delta)*delta
	time_cont=round(time_cont/delta)*delta
	
	knots=list()
	for(p in 1:P){
		K=nknots[p]+4
		qs=1/(nknots[p]+1)
		knots[[p]]=quantile(time_unique[[p]],qs)
		if(nknots[p]>1){
			for(q in 2:nknots[p]){
				knots[[p]]=c(knots[[p]],q*quantile(time_unique[[p]],qs))
			}
		}	
		knots[[p]]=as.vector(knots[[p]])
	}

	phi_t_cont=list()
	for(p in 1:P){
		phi_t_cont[[p]]=bs(time_cont,knots=knots[[p]],degree=3,intercept=TRUE)
		temp=phi_t_cont[[p]]
		for(k in 1:(nknots[p]+4)){
			if(orth==TRUE){
				if(k>1){
					for(q in 1:(k-1)){
						temp[,k]=temp[,k]-(sum(temp[,k]*temp[,k-q])/
							sum(temp[,k-q]^2))*temp[,k-q];
					}
				}
			}		
		    temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
		}
		phi_t_cont[[p]]=t(sqrt(1/delta)*temp)
	}

	if(plotit==TRUE){
		par(mfrow=c(1,P))
		for(p in 1:P){
			plot(time_cont,phi_t_cont[[p]][1,],type="l",
						ylim=c(min(phi_t_cont[[p]]),max(phi_t_cont[[p]])))
			for(k in 2:dim(phi_t_cont[[p]])[1]){
				lines(time_cont,phi_t_cont[[p]][k,],type="l")
			}
		}	
	}	

	phi_t=list()
	for(p in 1:P){
		phi_t[[p]]=list()
		for(i in 1:N){
			phi_t[[p]][[i]]=array(0,dim=c((nknots[p]+4),length(time_sparse[[i]][[p]])))
			for(k in 1:(nknots[p]+4)){
				for(t in 1:length(time_sparse[[i]][[p]])){
					phi_t[[p]][[i]][k,t]=phi_t_cont[[p]][k,abs(time_cont-time_sparse[[i]][[p]][t])==
						min(abs(time_cont-time_sparse[[i]][[p]][t]))][1]
				}
			}
		}		
	}

	time_sparse_combined=list()
	for(i in 1:N){
		time_sparse_combined[[i]]=time_sparse[[i]][[1]]
		for(p in 2:P){
			time_sparse_combined[[i]]=c(time_sparse_combined[[i]],time_sparse[[i]][[p]])
		}
		time_sparse_combined[[i]]=sort(unique(time_sparse_combined[[i]]))
	}
	phi_t_combined=list()
	for(p in 1:P){
		phi_t_combined[[p]]=list()
		for(i in 1:N){
			phi_t_combined[[p]][[i]]=array(0,dim=c((nknots[p]+4),length(time_sparse_combined[[i]])))
			for(k in 1:(nknots[p]+4)){
				for(t in 1:length(time_sparse_combined[[i]])){
					phi_t_combined[[p]][[i]][k,t]=phi_t_cont[[p]][k,abs(time_cont-
						time_sparse_combined[[i]][t])==min(abs(time_cont-
							time_sparse_combined[[i]][t]))][1]
				}
			}
		}		
	}
	Phi_D_combined=list()
	Phi_D_T_combined=list()
	for(i in 1:N){
		Phi_D_combined[[i]]=list()
		for(t in 1:length(time_sparse_combined[[i]])){
		    Phi_D_combined[[i]][[t]]=rbind(phi_t_combined[[1]][[i]][,t])
		    for(p in 2:P){
		    	Phi_D_combined[[i]][[t]]=blkdiag(Phi_D_combined[[i]][[t]],rbind(phi_t_combined[[p]][[i]][,t]))
		    }
		}
		Phi_D_T_combined[[i]]=Phi_D_combined[[i]][[1]]
		for(t in 2:length(time_sparse_combined[[i]])){
		    Phi_D_T_combined[[i]]=rbind(Phi_D_T_combined[[i]],Phi_D_combined[[i]][[t]])
		}
	}
			

	results=list()
	results[[1]]=knots
	results[[2]]=time_sparse
	results[[3]]=phi_t
	results[[4]]=time_cont
	results[[5]]=phi_t_cont
	results[[6]]=time_sparse_combined	
	results[[7]]=Phi_D_combined
	results[[8]]=Phi_D_T_combined

	return(results)	
}	
	
	
	
