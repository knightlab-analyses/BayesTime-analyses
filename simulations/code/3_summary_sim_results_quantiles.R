#----------------------------------------------------------------------
# summarize simulation results (with quantiles)
# estimation of mean curve (MU)
# estimation of FPC curve (THETA)
#----------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

nsims = as.numeric(Sys.getenv("PBS_ARRAYID"))
sim_list <- c('N100_M0', 'N100_M20', 'N100_M50', 'N100_M80',
              'N50_M0', 'N50_M20', 'N50_M50', 'N50_M80',
              'N25_M0', 'N25_M20', 'N25_M50', 'N25_M80',
              'N10_M0', 'N10_M20', 'N10_M50', 'N10_M80')

for (dat_name in sim_list){
    print(paste0('dat_name:', dat_name))
    error_MU_list = error_THETA_list = NULL
    
    for (idx_sim in 1:nsims){
        print(paste0('idx_sim:', idx_sim))
        if (dat_name %in% c('N100_M0', 'N100_M20', 'N100_M50', 'N100_M80', 'N50_M0')){
            tmp <- try(load(paste0(dir, "/summary_results/concise_sfpca_results_", dat_name, '_', idx_sim, ".RData")), TRUE)
        }else{    
            tmp <- try(load(paste0(dir, "/summary_results/sfpca_results_", dat_name, '_', idx_sim, ".RData")), TRUE)
        }    
        if (class(tmp) == 'try-error'){
            print(paste0('error in: ', dat_name, ' idx_sim: ', idx_sim))
            MU_error = THETA_error = NA 
        }

        error_MU_list = cbind(error_MU_list, MU_error)
        error_THETA_list = cbind(error_THETA_list, THETA_error)
    }

    assign(paste0('mean_MU_', dat_name), mean(error_MU_list, na.rm=T))
    assign(paste0('q025_MU_', dat_name), quantile(error_MU_list, probs = 0.025, na.rm=T)) 
    assign(paste0('q975_MU_', dat_name), quantile(error_MU_list, probs = 0.975, na.rm=T)) 
    assign(paste0('mean_THETA_', dat_name), mean(error_THETA_list, na.rm=T))
    assign(paste0('q025_THETA_', dat_name), quantile(error_THETA_list, probs = 0.025, na.rm=T)) 
    assign(paste0('q975_THETA_', dat_name), quantile(error_THETA_list, probs = 0.975, na.rm=T)) 
}


#----------------------------------------------------------------------
# summarize estimation of mean curve (MU)
#----------------------------------------------------------------------
table_results_MU = matrix(NA, nrow=4, ncol=4)
colnames(table_results_MU) = c("N100", "N50", "N25", "N10")
rownames(table_results_MU) = c('MU_M0', 'MU_M20', 'MU_M50', 'MU_M80')

table_results_MU[, 1] = c(paste0(round(mean_MU_N100_M0,3), ' (', 
                                 round(q025_MU_N100_M0,3), ', ',
                                 round(q975_MU_N100_M0,3), ')'),
                          paste0(round(mean_MU_N100_M20,3), ' (', 
                                 round(q025_MU_N100_M20,3), ', ',
                                 round(q975_MU_N100_M20,3), ')'),
                          paste0(round(mean_MU_N100_M50,3), ' (', 
                                 round(q025_MU_N100_M50,3), ', ',
                                 round(q975_MU_N100_M50,3), ')'),
                          paste0(round(mean_MU_N100_M80,3), ' (', 
                                 round(q025_MU_N100_M80,3), ', ',
                                 round(q975_MU_N100_M80,3), ')'))

table_results_MU[, 2] = c(paste0(round(mean_MU_N50_M0,3), ' (', 
                                 round(q025_MU_N50_M0,3), ', ',
                                 round(q975_MU_N50_M0,3), ')'),
                          paste0(round(mean_MU_N50_M20,3), ' (', 
                                 round(q025_MU_N50_M20,3), ', ',
                                 round(q975_MU_N50_M20,3), ')'),
                          paste0(round(mean_MU_N50_M50,3), ' (', 
                                 round(q025_MU_N50_M50,3), ', ',
                                 round(q975_MU_N50_M50,3), ')'),
                          paste0(round(mean_MU_N50_M80,3), ' (', 
                                 round(q025_MU_N50_M80,3), ', ',
                                 round(q975_MU_N50_M80,3), ')'))

table_results_MU[, 3] = c(paste0(round(mean_MU_N25_M0,3), ' (', 
                                 round(q025_MU_N25_M0,3), ', ',
                                 round(q975_MU_N25_M0,3), ')'),
                          paste0(round(mean_MU_N25_M20,3), ' (', 
                                 round(q025_MU_N25_M20,3), ', ',
                                 round(q975_MU_N25_M20,3), ')'),
                          paste0(round(mean_MU_N25_M50,3), ' (', 
                                 round(q025_MU_N25_M50,3), ', ',
                                 round(q975_MU_N25_M50,3), ')'),
                          paste0(round(mean_MU_N25_M80,3), ' (', 
                                 round(q025_MU_N25_M80,3), ', ',
                                 round(q975_MU_N25_M80,3), ')'))

table_results_MU[, 4] = c(paste0(round(mean_MU_N10_M0,3), ' (', 
                                 round(q025_MU_N10_M0,3), ', ',
                                 round(q975_MU_N10_M0,3), ')'),
                          paste0(round(mean_MU_N10_M20,3), ' (', 
                                 round(q025_MU_N10_M20,3), ', ',
                                 round(q975_MU_N10_M20,3), ')'),
                          paste0(round(mean_MU_N10_M50,3), ' (', 
                                 round(q025_MU_N10_M50,3), ', ',
                                 round(q975_MU_N10_M50,3), ')'),
                          paste0(round(mean_MU_N10_M80,3), ' (', 
                                 round(q025_MU_N10_M80,3), ', ',
                                 round(q975_MU_N10_M80,3), ')'))

print('summarize estimation of mean curve finished')

#----------------------------------------------------------------------
# estimation of FPC curve (THETA)
#----------------------------------------------------------------------
table_results_THETA = matrix(NA, nrow=4, ncol=4)
colnames(table_results_THETA) = c("N100", "N50", "N25", "N10")
rownames(table_results_THETA) = c('THETA_M0', 'THETA_M20', 'THETA_M50', 'THETA_M80')

table_results_THETA[, 1] = c(paste0(round(mean_THETA_N100_M0,3), ' (', 
                                 round(q025_THETA_N100_M0,3), ', ',
                                 round(q975_THETA_N100_M0,3), ')'),
                          paste0(round(mean_THETA_N100_M20,3), ' (', 
                                 round(q025_THETA_N100_M20,3), ', ',
                                 round(q975_THETA_N100_M20,3), ')'),
                          paste0(round(mean_THETA_N100_M50,3), ' (', 
                                 round(q025_THETA_N100_M50,3), ', ',
                                 round(q975_THETA_N100_M50,3), ')'),
                          paste0(round(mean_THETA_N100_M80,3), ' (', 
                                 round(q025_THETA_N100_M80,3), ', ',
                                 round(q975_THETA_N100_M80,3), ')'))

table_results_THETA[, 2] = c(paste0(round(mean_THETA_N50_M0,3), ' (', 
                                 round(q025_THETA_N50_M0,3), ', ',
                                 round(q975_THETA_N50_M0,3), ')'),
                          paste0(round(mean_THETA_N50_M20,3), ' (', 
                                 round(q025_THETA_N50_M20,3), ', ',
                                 round(q975_THETA_N50_M20,3), ')'),
                          paste0(round(mean_THETA_N50_M50,3), ' (', 
                                 round(q025_THETA_N50_M50,3), ', ',
                                 round(q975_THETA_N50_M50,3), ')'),
                          paste0(round(mean_THETA_N50_M80,3), ' (', 
                                 round(q025_THETA_N50_M80,3), ', ',
                                 round(q975_THETA_N50_M80,3), ')'))

table_results_THETA[, 3] = c(paste0(round(mean_THETA_N25_M0,3), ' (', 
                                 round(q025_THETA_N25_M0,3), ', ',
                                 round(q975_THETA_N25_M0,3), ')'),
                          paste0(round(mean_THETA_N25_M20,3), ' (', 
                                 round(q025_THETA_N25_M20,3), ', ',
                                 round(q975_THETA_N25_M20,3), ')'),
                          paste0(round(mean_THETA_N25_M50,3), ' (', 
                                 round(q025_THETA_N25_M50,3), ', ',
                                 round(q975_THETA_N25_M50,3), ')'),
                          paste0(round(mean_THETA_N25_M80,3), ' (', 
                                 round(q025_THETA_N25_M80,3), ', ',
                                 round(q975_THETA_N25_M80,3), ')'))

table_results_THETA[, 4] = c(paste0(round(mean_THETA_N10_M0,3), ' (', 
                                 round(q025_THETA_N10_M0,3), ', ',
                                 round(q975_THETA_N10_M0,3), ')'),
                          paste0(round(mean_THETA_N10_M20,3), ' (', 
                                 round(q025_THETA_N10_M20,3), ', ',
                                 round(q975_THETA_N10_M20,3), ')'),
                          paste0(round(mean_THETA_N10_M50,3), ' (', 
                                 round(q025_THETA_N10_M50,3), ', ',
                                 round(q975_THETA_N10_M50,3), ')'),
                          paste0(round(mean_THETA_N10_M80,3), ' (', 
                                 round(q025_THETA_N10_M80,3), ', ',
                                 round(q975_THETA_N10_M80,3), ')'))

print('summarize estimation of FPC curve finished')


save(table_results_MU, table_results_THETA, 
     file= paste0(dir, "/simulation_results_quantiles_nsims_", nsims, ".RData"))
