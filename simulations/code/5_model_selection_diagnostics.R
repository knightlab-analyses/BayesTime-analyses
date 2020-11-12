args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

library(BayesTime)
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
                                     PC_range = c(1:3), nknot_range = c(1:3), seed=idx_sim)

#-------------------------------
# save results
#-------------------------------
optimal_model_idx <- optimal(sfpca_stan_results_bayes)

optimal_model <- sfpca_stan_results_bayes[[optimal_model_idx]]
save(prepared_data, optimal_model, file=paste0(dir, '/bayesTime_modelSelection_', dat_name, '_nsim', idx_sim, '.RData'))

pdf(paste0('figures/k_diagnostics_', dat_name, '_nsim', idx_sim, '.pdf'), width=4, height=4)
plot_k_diagnostic(prepared_data, optimal_model)$figure
dev.off()
#ggsave(paste0('figures/k_diagnostics_', dat_name, '_nsim', idx_sim, '.png'), width=4, height=4, dpi=300)

pdf(paste0('figures/posterior_', dat_name, '_nsim', idx_sim, '.pdf'), width=4, height=4)
plot_posterior_diagnostic(prepared_data, optimal_model)$figure
dev.off()
#ggsave(paste0('figures/posterior_', dat_name, '_nsim', idx_sim, '.png'), width=4, height=4, dpi=300)

