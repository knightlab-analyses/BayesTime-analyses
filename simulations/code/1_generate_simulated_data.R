################################
### generate simulated data ###
###############################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

source('sfpca_sim_setting.R')
source('simulate_data_sparse.R')
source('basis_setup_sparse.R')

#-------------------------
# simulation parameters
#-------------------------
param_list <- list(list(N=100, N_T_mean=10), # 100 subjects with 0% missing
                   list(N=100, N_T_mean=8), # 100 subjects with 20% missing
                   list(N=100, N_T_mean=5), # 100 subjects with 50% missing
                   list(N=100, N_T_mean=2), # 100 subjects with 80% missing
                   list(N=50, N_T_mean=10), # 50 subjects with 0% missing
                   list(N=50, N_T_mean=8), # 50 subjects with 20% missing
                   list(N=50, N_T_mean=5), # 50 subjects with 50% missing
                   list(N=50, N_T_mean=2), # 50 subjects with 80% missing
                   list(N=25, N_T_mean=10), # 25 subjects with 0% missing
                   list(N=25, N_T_mean=8), # 25 subjects with 20% missing
                   list(N=25, N_T_mean=5), # 25 subjects with 50% missing
                   list(N=25, N_T_mean=2), # 25 subjects with 80% missing
                   list(N=10, N_T_mean=10), # 10 subjects with 0% missing
                   list(N=10, N_T_mean=8), # 10 subjects with 20% missing
                   list(N=10, N_T_mean=5), # 10 subjects with 50% missing
                   list(N=10, N_T_mean=2)) # 10 subjects with 80% missing
                   
i <- as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(i)
sim_param = sim_data(N=param_list[[i]]$N, N_T_max=10, N_T_mean=param_list[[i]]$N_T_mean)

#-------------------------
# generate data
#-------------------------
TIME_SPARSE=list()
Y_SPARSE=list()
MU_SPARSE=list()
F_SPARSE=list()
OMEGA_SPARSE=list()
ALPHA=list()
phi_t = list()

nsims=1000
for(ii in 1:nsims){
	print(paste0('nsims: ', ii))
    simdata=simulate_data_sparse(T_range=sim_param$T_range, N_T_max=sim_param$N_T_max, N_T_mean=sim_param$N_T_mean,
    	                         P=sim_param$P, N=sim_param$N, params=sim_param$params, orth=sim_param$orth, nseeds=ii)
	TIME_SPARSE[[ii]]=simdata[[1]]
	Y_SPARSE[[ii]]=simdata[[2]]
	MU_SPARSE[[ii]]=simdata[[3]]
	F_SPARSE[[ii]]=simdata[[4]]
	OMEGA_SPARSE[[ii]]=simdata[[5]]
	ALPHA[[ii]]=simdata[[6]]
	phi_t[[ii]] = simdata[[9]]
}

save(nsims, TIME_SPARSE, Y_SPARSE, MU_SPARSE, F_SPARSE, OMEGA_SPARSE, ALPHA, phi_t, sim_param, 
	 file=paste0(dir, '/sim_data/sim_dat_', 'N', sim_param$N, '_M', sim_param$Missing, '.Rdata'))



