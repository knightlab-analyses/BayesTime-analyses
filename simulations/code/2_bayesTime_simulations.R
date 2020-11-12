args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

library(BayesTime)
library(Matrix) # need to load this for t() to work on sparse matrix
source('sfpca_sim_setting.R')

sim_list <- c('N100_M0', 'N100_M20', 'N100_M50', 'N100_M80',
              'N50_M0', 'N50_M20', 'N50_M50', 'N50_M80',
              'N25_M0', 'N25_M20', 'N25_M50', 'N25_M80',
              'N10_M0', 'N10_M20', 'N10_M50', 'N10_M80')
sim_setting <- expand.grid(repetitions = 1:1000, scenarios = sim_list)

M <- as.numeric(Sys.getenv("PBS_ARRAYID"))
set.seed(M)

idx_sim = sim_setting[M, 1]
dat_name = sim_setting[M, 2]
load(paste0(dir, '/sim_data/sim_dat_', dat_name, '.Rdata'))


## data preparation
sim_data <- sim_data_reformat(Y_SPARSE=Y_SPARSE, TIME_SPARSE=TIME_SPARSE, MU_SPARSE=MU_SPARSE,
           F_SPARSE=F_SPARSE, OMEGA_SPARSE=OMEGA_SPARSE, ALPHA=ALPHA, 
           params=sim_param$params, nsim=idx_sim)
df = sim_data$df
Q1 = sim_data$Q1
Q = sim_data$Q
K1 = sim_data$K1
K = sim_data$K

prepared_data <- prepare_data(data = df, unique_subject_id = 'id', time_name = 'time',
                              response_name = 'response', scale_time = FALSE,
                              transform_y = 'no')

sfpca_stan_results_bayes <- stan_fit(sfpca_data = prepared_data, Nsamples = 1000, Nchain = 3,
                                     PC_range = sim_data$K, nknot_range = sim_data$nknots, seed=M)

#--------------------------------------------
# check model results
#--------------------------------------------
optimal_model <- sfpca_stan_results_bayes[[1]] # 1 model only

# plot_k_diagnostic(prepared_data, optimal_model)$figure 
# ggsave(paste0('bayesTime_N100_M80_trueModel_nsim', nsim, '_k_diagnostic.png'))

# plot_posterior_diagnostic(prepared_data, optimal_model)$figure # nice
# ggsave(paste0('bayesTime_N100_M80_trueModel_nsim', nsim, '_posterior.png'))

model_output <- output_results(sfpca_data = prepared_data, model = optimal_model)

#----------------------------------------------
# extract error in estimation of MU & THETA
#----------------------------------------------
MU_array = model_output$rotation$theta_mu_new
THETA_array = model_output$rotation$Theta_new

nloop=dim(MU_array)[2]
first=1
last=nloop

MU_mean = MU_array[, first]
THETA_mean = THETA_array[,,first] 
for(iter in 2:nloop){
  MU_mean = MU_mean + MU_array[, iter]
  THETA_mean = THETA_mean + THETA_array[,,iter]
}
MU_mean=cbind(MU_mean/(last-first+1))
THETA_mean = cbind(THETA_mean/(last-first+1)) 

MU_true = sim_param$params[[7]][(1+Q1):(Q+Q1)]
THETA_true = sim_param$params[[6]][(1+Q1):(Q1+Q), (1+K1): (K+K1)]

MU_error = mean((MU_mean - MU_true)^2)
THETA_error = mean((THETA_mean - THETA_true)^2)

#----------------------------------------------
# extract predicted vals
#----------------------------------------------
sfpca_data = prepared_data
unique_subject_id = 'id'
data = df
time_name = 'time'
response_name = 'response'
error_predict = NULL
for (idx_each in unique(df$id)){
  subject_starts = sfpca_data$visits.start
  subject_stops = sfpca_data$visits.stop
  subject_idx = which(unique(sfpca_data$data$ID) %in% idx_each)
  fit_curve = rstan::summary(optimal_model$sa, "Ynew")$summary[
                     subject_starts[subject_idx]: subject_stops[subject_idx], c('mean')]

  sub_data = data[data[, unique_subject_id] %in% idx_each,]
  sub_data = sub_data[order(sub_data[, time_name]), ]
  response_obs = sub_data[, response_name]

  error_each_subject = mean((fit_curve - response_obs)^2) # mean squared error
  error_predict = cbind(error_predict, error_each_subject)
}

# save minimum results
save(error_predict, MU_error, THETA_error,
     file= paste0(dir, "/summary_results/concise_sfpca_results_", dat_name, '_', idx_sim, ".RData"))  

# save more results
# save(sim_data, sfpca_stan_results_bayes, sfpca_data, data, error_predict,
#      MU_true, THETA_true, MU_mean, THETA_mean, MU_error, THETA_error,
#      file= paste0(dir, "/summary_results/sfpca_results_", dat_name, '_', idx_sim, ".RData"))    
