library(gridExtra)
library(cowplot) # for plot_grid()
library(gridGraphics)
library(ggplot2)

ggplot_results <- function(dir, dat_name, idx_sim, title_mean, fpc_range){
    load(paste0(dir, "/summary_results/sfpca_results_", dat_name, '_', idx_sim, ".RData"))
    load(paste0(dir, '/sim_data/sim_dat_', dat_name, '.Rdata'))
    
    results_basis = BayesTime::basis_setup_sparse(sfpca_data=sfpca_data, nknots=sim_data$nknots, orth=TRUE)
    time_cont = results_basis$time_cont
    phi_t_cont = results_basis$orth_spline_basis_cont
    phi_transform = Matrix::bdiag(cbind(phi_t_cont))

    optimal_model <- sfpca_stan_results_bayes[[1]]
    model_output <- BayesTime::output_results(sfpca_data = sfpca_data, model = optimal_model)
    Mu_functions=model_output$Mu_functions

    Mu_true_functions=Matrix::t(phi_transform)%*%MU_true
    Mu1_true=Mu_true_functions[1:length(time_cont)]
    Mu1=Mu_functions[1:length(time_cont)]

    FPC1_mean=model_output$FPC_mean
    time_unique = unique(sim_data$df$time)
    time_unique = time_unique[order(time_unique)]

    K1 = sim_data$K1
    K = sim_data$K
    FPC_true = sim_param$params[[8]][, (1+K1):(K+K1)]

    N = sim_data$N
    J = sim_data$J
    time_sparse = Y_sparse = list()
    for (i in 1:N){
      time_sparse[[i]] = TIME_SPARSE[[idx_sim]][[i]][[J]]
      Y_sparse[[i]] = Y_SPARSE[[idx_sim]][[i]][[J]]
    }
    
    #--------------------------
    # population mean curve
    #--------------------------
    figure_mean <- ggplot() +
        ylim(-10, 10) +
        labs(title= title_mean, x = 'Time', y = 'Response') +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
              axis.text.x = element_text(size = 10, face = "bold"),
              axis.text.y = element_text(size = 10, face = "bold"),
              axis.title.x = element_text(size = 12, face = "bold"),
              axis.title.y = element_text(size = 12, face = "bold"),
              legend.title = element_blank(), legend.position ='top')

      for (i in 1:N){
        figure_mean <- figure_mean + geom_line(aes_string(x = time_sparse[[i]],y = Y_sparse[[i]]), lwd = 0.1, color=1)
      }

    figure_mean <- figure_mean + geom_line(aes(x = time_cont, y = Mu1_true, color='True'), lwd = 1)
    figure_mean <- figure_mean + geom_line(aes(x = time_cont, y = Mu1, color='Estimated'), lwd = 1)
    figure_mean <- figure_mean + 
                     scale_color_manual(name='Color', values = c("True" = "blue", "Estimated" = "red"))

    #--------------------------
    # FPC curve
    #--------------------------
    figure_fpc <- ggplot() +
    ylim(fpc_range) +
    labs(x = 'Time', y = 'PC curve') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), legend.position ='top') 

    figure_fpc <- figure_fpc + geom_line(aes(x = time_unique, y = FPC_true[,1], color='True', linetype = "pc1"), lwd = 1)
    figure_fpc <- figure_fpc + geom_line(aes(x = time_unique, y = FPC_true[,2], color='True', linetype = "pc2"), lwd = 1)
    figure_fpc <- figure_fpc + geom_line(aes(x = time_cont, y = FPC1_mean[,1], 
                                             color='Estimated', linetype = "pc1"), lwd = 1)
    figure_fpc <- figure_fpc + geom_line(aes(x = time_cont, y = FPC1_mean[,2], 
                                             color='Estimated', linetype = "pc2"), lwd = 1) 
    figure_fpc <- figure_fpc + scale_color_manual(name='Color', values = c("True" = "black", "Estimated" = "red")) +
                             scale_linetype_manual(name='Linetype', 
                                                   values = c("pc1"="solid", "pc2"="dashed"),
                                                   guide = guide_legend(override.aes = list(
                                                       color = c("black", "black"),size = c(.5, .5)))) 

    return(list('figure_mean'=figure_mean, 'figure_fpc'=figure_fpc))
}
