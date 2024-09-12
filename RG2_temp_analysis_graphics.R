
packages <- c("magrittr", "dplyr", "here", "metafor", "patchwork")

# check, whether library already installed or not - install and load as needed:
apply(as.matrix(packages), MARGIN = 1, FUN = function(x) {
  
  pkg_avail <- nzchar(system.file(package = x))   # check if library is installed on system
  
  if(pkg_avail){
    require(x, character.only = TRUE)             # load the library, if already installed
    
  }else{
    install.packages(x)                           # install the library, if missing
    require(x, character.only = TRUE)             # load after installation
  }
})


full_paths <- list.files(here("Data/Extracted (Project) Data"), full.names = T)

MASC_names <- substr(list.files(here("Data/Extracted (Project) Data")), 1, nchar(list.files(here("Data/Extracted (Project) Data")))-4) 

data.list <- lapply(full_paths, read.csv)

agg_L <- lapply(data.list, FUN = function(x){
  
  if("group" %in% names(x)){
    
    MASC_MD_L <- lapply(unique(x$source), FUN = function(lab){
      
      d <- x[x$source == lab,]
      
      d1 <- na.omit(d[d$group == 1, -c(grep("group", names(d)), grep("source", names(d)))])
      d0 <- na.omit(d[d$group == 0, -c(grep("group", names(d)), grep("source", names(d)))])
      
      n1 <- nrow(d1)
      n0 <- nrow(d0)
      
      C1 <- cov(d1)
      C0 <- cov(d0)
      
      j <- dim(C1)[1]
      
      alpha1 <- (1 - sum(diag(C1))/sum(C1)) * (j/(j - 1))
      alpha0 <- (1 - sum(diag(C0))/sum(C0)) * (j/(j - 1))
      
      SE_Balpha1 <- sqrt((2 * j)/((j - 1) * (n1 - 2))) 
      SE_Balpha0 <- sqrt((2 * j)/((j - 1) * (n1 - 2))) 
      
      mv1 <- rowMeans(d1)
      mv0 <- rowMeans(d0)
      
      m1 <- mean(mv1)
      m0 <- mean(mv0)
      
      sd1 <- sqrt(mean((mv1 - m1)^2))
      sd0 <- sqrt(mean((mv0 - m0)^2))
      
      sigma1 <- sd(mv1)
      sigma0 <- sd(mv0)
      
      MD <- m1-m0
      
      pooled_sd <- sqrt(((n1-1)*sd1^2 + (n0-1)*sd0^2)/(n1+n0-2))
      
      pooled_sd_corr <- sqrt(((n1-1)*((sd1^2)*alpha1) + (n0-1)*((sd0^2)*alpha0))/(n1+n0-2))
      
      d_raw <- MD/pooled_sd
      d_corr <- MD/pooled_sd_corr
      
      SE_d <- sqrt(((n0 + n1)/(n0 * n1)) + ((d_raw^2) / (2 * (n0 + n1))))
      SE_d_corr <- sqrt(SE_d^2 * (d_corr/d_raw)^2)
      
      return(data.frame(n1 = n1, 
                        n0 = n0,
                        j = j,
                        m1 = m1,
                        m0 = m0,
                        sd1 = sd1,
                        sd0 = sd0,
                        sigma1 = sigma1,
                        sigma0 = sigma0,
                        alpha1 = alpha1,
                        alpha0 = alpha0,
                        SE_Balpha1 = SE_Balpha1,
                        SE_Balpha0 = SE_Balpha0,
                        MD = MD,
                        pooled_sd = pooled_sd,
                        pooled_sd_corr = pooled_sd_corr,
                        d_raw = d_raw,
                        d_corr = d_corr,
                        SE_d = SE_d,
                        SE_d_corr = SE_d_corr))
      
    })
    
    MASC_MD_df <- do.call(rbind, MASC_MD_L)
    
    return(MASC_MD_df)
    
  }else{
    
    MASC_al_L <- lapply(unique(x$source), FUN = function(lab){
      
      d <- x[x$source == lab,]
      
      d1 <- na.omit(d[, -grep("source", names(d))])
      
      n1 <- nrow(d1)
      
      C1 <- cov(d1)
      
      j <- dim(C1)[1]
      
      alpha1 <- (1 - sum(diag(C1))/sum(C1)) * (j/(j - 1))
      
      SE_Balpha1 <- sqrt((2 * j)/((j - 1) * (n1 - 2))) 
      
      mv1 <- rowMeans(d1)
      
      m1 <- mean(mv1)
      
      sd1 <- sqrt(mean((mv1 - m1)^2))
      
      sigma1 <- sd(mv1)
      
      
      return(data.frame(n1 = n1, 
                        n0 = NA,
                        j = j,
                        m1 = m1,
                        m0 = NA,
                        sd1 = sd1,
                        sd0 = NA,
                        sigma1 = sigma1,
                        sigma0 = NA,
                        alpha1 = alpha1,
                        alpha0 = NA,
                        SE_Balpha1 = SE_Balpha1,
                        SE_Balpha0 = NA,
                        MD = NA,
                        pooled_sd = NA,
                        pooled_sd_corr = NA,
                        d_raw = NA,
                        d_corr = NA,
                        SE_d = NA,
                        SE_d_corr = NA))
      
    })
    
    MASC_al_df <- do.call(rbind, MASC_al_L)
    
    return(MASC_al_df)
    
  }
  
})


var_Bonnett_backtransformed <- function(rma_obj){
  (((-exp(rma_obj$b[1]))^2) * rma_obj$tau2) + (.5*((-exp(rma_obj$b[1]))^2)*(rma_obj$tau2^2)) + ((-exp(rma_obj$b[1])) * (-exp(rma_obj$b[1])) * (rma_obj$tau2^2))
}

mean_Bonnett_backtransformed <- function(rma_obj){
  1 - exp(rma_obj$b[1]) + ((-exp(rma_obj$b[1])) / 2) * rma_obj$tau2
}


B_alpha_rma <- lapply(agg_L, FUN = function(x){
  
  if(!is.na(x$n0[1])){
    
    alpha_pooled <- (x$n1 / (x$n0 + x$n1)) * x$alpha1 + (x$n0 / (x$n0 + x$n1)) * x$alpha0
    
    SE_B_alpha_pooled <- sqrt((2 * x$j)/((x$j - 1) * (x$n1 + x$n0 - 2))) 
    
    BA_rma <- metafor::rma(yi = log(1 - alpha_pooled),
                           sei = SE_B_alpha_pooled,
                           measure = "GEN",
                           method = "REML")
    
  }else{
    
    BA_rma <- metafor::rma(yi = log(1 - x$alpha1),
                           sei = x$SE_Balpha1,
                           measure = "GEN",
                           method = "REML")
    
    
  }
  
  
  
  return(data.frame(tau2_BA = BA_rma$tau2,
                    mu_BA = BA_rma$b[1],
                    tau2_alpha = var_Bonnett_backtransformed(BA_rma),
                    mu_alpha = mean_Bonnett_backtransformed(BA_rma),
                    k = nrow(x),
                    I2 = BA_rma$I2,
                    H2 = BA_rma$H2))
  
})

B_alpha_rma_df <- do.call(rbind, B_alpha_rma) %>% 
  mutate(MASC = MASC_names)


library(ggplot2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(3)

ggplot(B_alpha_rma_df) +
  geom_violin(aes(x = 1, y = sqrt(tau2_alpha)), fill = "#00BA38", colour = "#00BA38", alpha = .3) +
  geom_boxplot(aes(x = 1, y = sqrt(tau2_alpha)), fill = "#00BA38", alpha = .3, width = .2, outliers = FALSE) +
  geom_point(aes(x = 1, y = sqrt(tau2_alpha)), position = position_jitter(width = .15), shape = 21,
             fill = "#00BA38", alpha = .5, size = 3) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 15)) +
  labs(x = "Reliability", y = expression(tau[r["xx'"]]))

ggsave(here("Graphics/rel_violin.png"),
       plot = last_plot(),
       height = 5.5, width = 5)

ggplot(B_alpha_rma_df) +
  geom_violin(aes(x = 1, y = mu_alpha), fill = "#00BA38", colour = "#00BA38", alpha = .3) +
  geom_boxplot(aes(x = 1, y = mu_alpha), fill = "#00BA38", alpha = .3, width = .2, outliers = FALSE) +
  geom_point(aes(x = 1, y = mu_alpha), position = position_jitter(width = .15), shape = 21,
             fill = "#00BA38", alpha = .5, size = 3) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(colour = "grey"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 15)) +
  labs(x = "Reliability", y = expression(mu[r["xx'"]]))


ggsave(here("Graphics/rel_violin2.png"),
       plot = last_plot(),
       height = 5.5, width = 5)


MASC_names

agg_L[[37]]$d_raw <- -agg_L[[37]]$d_raw 
agg_L[[37]]$d_corr <- -agg_L[[37]]$d_corr 

d_rma_L <- lapply(agg_L, FUN = function(x){
  
  if(!is.na(x$n0[1])){
    
    d_rma_raw <- metafor::rma(yi = x$d_raw,
                              sei = x$SE_d,
                              measure = "GEN",
                              method = "REML")
    
    d_rma_corr <- metafor::rma(yi = x$d_cor,
                               sei = x$SE_d_corr,
                               measure = "GEN",
                               method = "REML")
    
    return(data.frame(mu_d_raw = d_rma_raw$b[1],
                      mu_d_corr = d_rma_corr$b[1],
                      tau_d_raw = sqrt(d_rma_raw$tau2),
                      tau_d_corr = sqrt(d_rma_corr$tau2),
                      k = nrow(x),
                      I2_raw = d_rma_raw$I2,
                      I2_corr = d_rma_corr$I2,
                      H2_raw = d_rma_raw$H2,
                      H2_corr = d_rma_corr$H2,
                      p_raw = d_rma_raw$QEp,
                      p_corr = d_rma_corr$QEp))
    
  }else{
    
    return(data.frame(mu_d_raw = NA,
                      mu_d_corr = NA,
                      tau_d_raw = NA,
                      tau_d_corr = NA,
                      k = NA,
                      I2_raw = NA,
                      I2_corr = NA,
                      H2_raw = NA,
                      H2_corr = NA,
                      p_raw = NA,
                      p_corr = NA))
    
  }
  
})

d_rma_df <- do.call(rbind, d_rma_L) %>% 
  mutate(MASC = MASC_names)

d_rma_df2 <- data.frame(corr = as.factor(c(rep(0, sum(!is.na(d_rma_df$mu_d_raw))), 
                                 rep(1, sum(!is.na(d_rma_df$mu_d_raw))))),
                        tau = c(d_rma_df$tau_d_raw[!is.na(d_rma_df$mu_d_raw)], d_rma_df$tau_d_corr[!is.na(d_rma_df$mu_d_raw)]),
                        mu = abs(c(d_rma_df$mu_d_raw[!is.na(d_rma_df$mu_d_raw)], d_rma_df$mu_d_corr[!is.na(d_rma_df$mu_d_raw)])),
                        sig = as.factor(c(d_rma_df$p_raw[!is.na(d_rma_df$mu_d_raw)], d_rma_df$p_corr[!is.na(d_rma_df$mu_d_raw)]) < .05))

ggplot(d_rma_df2) +
  geom_violin(aes(x = corr, y = mu, colour = corr)) +
  geom_boxplot(aes(x = corr, y = mu, colour = corr)) +
  geom_point(aes(x = corr, y = mu), position = position_jitter(width = .1))


ggplot(d_rma_df2) +
  geom_violin(aes(x = corr, y = tau, colour = corr)) +
  geom_boxplot(aes(x = corr, y = tau, colour = corr)) +
  geom_point(aes(x = corr, y = tau, shape = sig), position = position_jitter(width = .1))

d_rma_df3 <- d_rma_df %>% 
  filter(d_rma_df$p_raw <= .05 | d_rma_df$p_corr <= .05) %>% 
  mutate(mu_d_raw = abs(mu_d_raw),
         mu_d_corr = abs(mu_d_corr))

d_rma_df3_l <- data.frame(sig = c(d_rma_df3$p_raw <= .05, d_rma_df3$p_corr <= .05),
                          mu_d = c(d_rma_df3$mu_d_raw, d_rma_df3$mu_d_corr),
                          tau_d = c(d_rma_df3$tau_d_raw, d_rma_df3$tau_d_corr),
                          corr = c(rep(0, nrow(d_rma_df3)), rep(1, nrow(d_rma_df3))),
                          MASC = c(d_rma_df3$MASC, d_rma_df3$MASC)) 


ggplot(data = d_rma_df3_l, aes(x = as.factor(corr), y = mu_d)) +
  geom_violin() +
  geom_boxplot() +
  geom_point(aes(shape = sig), position = position_jitter(width = .1))
  

ggplot(data = d_rma_df3_l, aes(x = corr, y = mu_d)) +
  geom_point(aes(colour = as.factor(MASC), shape = sig)) +
  geom_line(aes(colour = as.factor(MASC)))

ggplot(data = d_rma_df3_l, aes(x = corr, y = tau_d)) +
  geom_point(aes(colour = as.factor(MASC), shape = sig)) +
  geom_line(aes(colour = as.factor(MASC)))


d_rma_df_fin_l <- d_rma_df3_l %>% 
  filter(MASC %in% c("Nosek_Explicit_Art", "PSACR002_neg_photo", "PSACR001_anxiety_int", "PSACR001_behav_int"))


ggplot(data = d_rma_df_fin_l, aes(x = corr, y = tau_d)) +
  geom_point(aes(colour = as.factor(MASC))) +
  geom_line(aes(colour = as.factor(MASC)))




dich_index <- sapply(agg_L, FUN = function(x){!is.na(x$n0[1])})
MASC_names[dich_index]

effect_index <- MASC_names %in% c("Albarracin_Priming_SAT", "Alter_Analytic_Processing",
                                  "Carter_Flag_Priming", "Caruso_Currency_Priming",
                                  "Dijksterhuis_trivia", "Finkel_Exit_Forgiveness", 
                                  "Finkel_Neglect_Forgiveness",
                                  "Giessner_Vertical_Position", "Hart_Criminal_Intentionality",    
                                  "Hart_Detailed_Processing", 
                                  # "Hart_Intention_Attribution",       
                                  "Husnu_Imagined_Contact", "Nosek_Explicit_Art",
                                  "Nosek_Explicit_Math", "PSACR001_anxiety_int", 
                                  "PSACR001_behav_int", 
                                  "PSACR002_neg_photo", #"Shnabel_Willingness_Reconcile_Rev",
                                  "Shnabel_Willingness_Reconcile_RPP", "Srull_Behaviour_Hostility",
                                  "Srull_Ronald_Hostility", "Tversky_Directionality_Similarity1"
                                  #, "Zhong_Desirability_Cleaning"
                                  )


d_rma_full_L <- lapply(agg_L[effect_index], FUN = function(x){
  
  d_rma_raw <- metafor::rma(yi = x$d_raw,
                            sei = x$SE_d,
                            measure = "GEN",
                            method = "REML")
  
  d_rma_corr <- metafor::rma(yi = x$d_cor,
                             sei = x$SE_d_corr,
                             measure = "GEN",
                             method = "REML")
  
  return(list(rma_raw = d_rma_raw,
              rma_corr = d_rma_corr))
})


names(d_rma_full_L) <- MASC_names[effect_index]


ggplot() +
  geom_histogram(aes(x = y$rma_raw$yi), binwidth = .1, fill = "darkgrey") +
  geom_histogram(aes(x = y$rma_corr$yi), binwidth = .1, fill = "black") +
  geom_line(aes(x = seq(from = min(y$rma_raw$yi), to = max(y$rma_raw$yi), length.out = 100),
                y = dnorm(seq(from = min(y$rma_raw$yi), to = max(y$rma_raw$yi), length.out = 100),
                          mean = y$rma_raw$b[1], sd = sqrt(y$rma_raw$tau2))),
            linewidth = 2, colour = "red") +
  geom_line(aes(x = seq(from = min(y$rma_corr$yi), to = max(y$rma_corr$yi), length.out = 100),
                y = dnorm(seq(from = min(y$rma_corr$yi), to = max(y$rma_corr$yi), length.out = 100),
                          mean = y$rma_corr$b[1], sd = sqrt(y$rma_corr$tau2))),
            linewidth = 2, colour = "green")


library(patchwork)

plots <- lapply(1:length(d_rma_full_L), FUN = function(idx){
  
  name <- MASC_names[effect_index][idx]
  
  y <- d_rma_full_L[[idx]]
  
  minmaxsequence_raw <- seq(from = min(y$rma_raw$yi), to = max(y$rma_raw$yi), length.out = 100)
  minmaxsequence_corr <- seq(from = min(y$rma_corr$yi), to = max(y$rma_corr$yi), length.out = 100)
  
  weight_raw <- y$rma_raw$tau2 / (y$rma_raw$tau2 + y$rma_raw$vb[1])
  shrunk_raw <- weight_raw * y$rma_raw$yi + (1-weight_raw) * y$rma_raw$b[1]
  
  p <- ggplot() +
    scale_shape_identity() +
    geom_point(aes(x = shrunk_raw, y = 0), colour = "darkgrey", shape = 108, size = 3) +
    # geom_histogram(aes(x = y$rma_corr$yi), binwidth = .1, fill = "black") +
    geom_line(aes(x = minmaxsequence_raw,
                  y = dnorm(minmaxsequence_raw, mean = y$rma_raw$b[1], sd = sqrt(y$rma_raw$tau2))),
              linewidth = 1, colour = "darkgrey") +
    # geom_line(aes(x = minmaxsequence_corr,
    #               y = dnorm(minmaxsequence_corr, mean = y$rma_corr$b[1], sd = sqrt(y$rma_corr$tau2))),
    #           linewidth = 2, colour = "green") +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = "white")) +
    labs(subtitle = name)
  
  if(y$rma_raw$tau2 == 0){
    p + ylim(c(0,1))
  }else{
    p
  }
  
  
})

combined_plot <- Reduce('+', plots) + plot_layout(ncol = 4)
combined_plot

ggsave(filename = here("Graphics/densities_full.png"),
       plot = last_plot(),
       width = 11, height = 5)


plots2 <- lapply((1:length(d_rma_full_L))[c(12,14,15,16)], FUN = function(idx){
  
  name <- MASC_names[effect_index][idx]
  
  y <- d_rma_full_L[[idx]]
  
  minmaxsequence_raw <- seq(from = min(y$rma_raw$yi), to = max(y$rma_raw$yi), length.out = 100)
  minmaxsequence_corr <- seq(from = min(y$rma_corr$yi), to = max(y$rma_corr$yi), length.out = 100)
  
  weight_raw <- y$rma_raw$tau2 / (y$rma_raw$tau2 + y$rma_raw$vb[1])
  shrunk_raw <- weight_raw * y$rma_raw$yi + (1-weight_raw) * y$rma_raw$b[1]
  weight_corr <- y$rma_corr$tau2 / (y$rma_corr$tau2 + y$rma_corr$vb[1])
  shrunk_corr <- weight_corr * y$rma_corr$yi + (1-weight_corr) * y$rma_corr$b[1]
  
  
  p <- ggplot() +
    scale_shape_identity() +
    geom_point(aes(x = y$rma_raw$yi, y = 0), colour = "darkgrey", shape = 108, size = 3) +
    geom_point(aes(x = y$rma_corr$yi, y = 0), colour = "black", shape = 108, size = 3) +
    geom_point(aes(x = y$rma_raw$b[1], y = 0), colour = "darkgrey", shape = 108, size = 8) +
    geom_point(aes(x = y$rma_corr$b[1], y = 0), colour = "black", shape = 108, size = 8) +
    geom_line(aes(x = minmaxsequence_raw,
                  y = dnorm(minmaxsequence_raw, mean = y$rma_raw$b[1], sd = sqrt(y$rma_raw$tau2))),
              linewidth = 1, colour = "darkgrey") +
    geom_line(aes(x = minmaxsequence_corr,
                  y = dnorm(minmaxsequence_corr, mean = y$rma_corr$b[1], sd = sqrt(y$rma_corr$tau2))),
              linewidth = 1, colour = "black") +
    theme(axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = "white")
    ) +
    labs(subtitle = name)
  
  if(y$rma_raw$tau2 == 0){
    p + ylim(c(0,1))
  }else{
    p
  }
  
})

combined_plot2 <- Reduce('+', plots2) + plot_layout(ncol = 2)
combined_plot2

ggsave(filename = here("Graphics/densities_four.png"),
       plot = last_plot(),
       width = 11, height = 5)




agg_L[[37]]



forest_plot_rel <- function(rma.fit_raw, rma.fit_cor, rma.data, ci.lvl = .975){
  rma.data %>% 
    mutate(cil_raw = d_raw - (1-(1-ci.lvl)/2) * (SE_d_raw),
           ciu_raw = d_raw + (1-(1-ci.lvl)/2) * (SE_d_raw),
           cil_cor = d_cor - (1-(1-ci.lvl)/2) * (SE_d_cor),
           ciu_cor = d_cor + (1-(1-ci.lvl)/2) * (SE_d_cor)) %>% 
    arrange(desc(d_raw)) %>% 
    ggplot() +
    # point estimate of raw d
    geom_point(aes(x = d_raw, y = 1:nrow(rma.data)), colour = "darkgrey", size = 2) +
    # CI of point estimate of raw d
    geom_segment(aes(x = cil_raw, y = 1:nrow(rma.data), xend = ciu_raw, yend = 1:nrow(rma.data)), colour = "darkgrey") +
    geom_segment(aes(x = cil_raw, xend = cil_raw, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "darkgrey") +
    geom_segment(aes(x = ciu_raw, xend = ciu_raw, y =( 1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "darkgrey") +
    # point estiamte of corrected d
    geom_point(aes(x = d_cor, y = 1:nrow(rma.data)), colour = "black", size = 2) +
    # CI of point estimate of corrected d
    geom_segment(aes(x = cil_cor, y = 1:nrow(rma.data), xend = ciu_cor, yend = 1:nrow(rma.data)), colour = "black") +
    geom_segment(aes(x = cil_cor, xend = cil_cor, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "black") +
    geom_segment(aes(x = ciu_cor, xend = ciu_cor, y = (1:nrow(rma.data))+.3, yend = (1:nrow(rma.data))-.3), colour = "black") +
    # solid black line separating RMA-estimates from point estimates
    geom_abline(slope = 0, intercept = -1, colour = "black") +
    # adding diamond showing rma-estimate of raw d
    geom_polygon(data = data.frame(x = c(rma.fit_raw$ci.ub, rma.fit_raw$b[1], rma.fit_raw$ci.lb, rma.fit_raw$b[1]),
                                   y = c(-3, -3-.7, -3, -3+.7)),
                 aes(x = x, y = y),
                 fill = "darkgrey", colour = "darkgrey") +
    # adding diamond showing rma-estimate of corrected d
    geom_polygon(data = data.frame(x = c(rma.fit_cor$ci.ub, rma.fit_cor$b[1], rma.fit_cor$ci.lb, rma.fit_cor$b[1]),
                                   y = c(-3, -3-.7, -3, -3+.7)),
                 aes(x = x, y = y),
                 colour = "black", fill = "black") +
    # defining theme of plot (transparent background etc.)
    theme(legend.position = "bottom", 
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", colour = "transparent"), 
          panel.grid.major.y = element_line(colour = "transparent"),
          panel.grid.major.x = element_line(colour = "grey"),
          panel.grid.minor = element_line(colour = "transparent"),
          axis.ticks = element_line(colour = "grey"),
          strip.background = element_rect(fill = "transparent"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 18)) +
    scale_y_continuous(breaks = c(-3, 1:nrow(rma.data)), labels = c("RMA-estimate", rma.data$country)) +
    labs(x = "Cohen's d",
         y = "Country") 
}


example_names <- unique(read.csv(here("Data/Extracted (Project) Data/PSACR002_neg_photo.csv"))$source)

example_codes <- substr(example_names, nchar(example_names) - 2, nchar(example_names)) 

fp_rma_dat <- agg_L[[37]] %>% 
  mutate(SE_d_raw = SE_d,
         SE_d_cor = SE_d_corr,
         d_cor = d_corr,
         country = example_codes)
fp_rma.fit_raw <- d_rma_full_L$PSACR002_neg_photo$rma_raw
fp_rma.fit_corr <- d_rma_full_L$PSACR002_neg_photo$rma_corr


forest_plot_rel(fp_rma.fit_raw, fp_rma.fit_corr, fp_rma_dat)

ggsave(filename = here("Graphics/forest_plot.png"),
       plot = last_plot(),
       width = 10,
       height = 6)





minmaxsequence_raw <- seq(from = min(fp_rma.fit_raw$yi), to = max(fp_rma.fit_raw$yi), length.out = 100)
minmaxsequence_corr <- seq(from = min(fp_rma.fit_corr$yi), to = max(fp_rma.fit_corr$yi), length.out = 100)

weight_raw <- fp_rma.fit_raw$tau2 / (fp_rma.fit_raw$tau2 + fp_rma.fit_raw$vb[1])
shrunk_raw <- weight_raw * fp_rma.fit_raw$yi + (1-weight_raw) * fp_rma.fit_raw$b[1]
weight_corr <- fp_rma.fit_corr$tau2 / (fp_rma.fit_corr$tau2 + fp_rma.fit_corr$vb[1])
shrunk_corr <- weight_corr * fp_rma.fit_corr$yi + (1-weight_corr) * fp_rma.fit_corr$b[1]


p <- ggplot() +
  scale_shape_identity() +
  geom_point(aes(x = fp_rma.fit_raw$yi, y = 0), colour = "darkgrey", shape = 108, size = 3) +
  geom_point(aes(x = fp_rma.fit_corr$yi, y = 0), colour = "black", shape = 108, size = 3) +
  geom_point(aes(x = fp_rma.fit_raw$b[1], y = 0), colour = "darkgrey", shape = 108, size = 8) +
  geom_point(aes(x = fp_rma.fit_corr$b[1], y = 0), colour = "black", shape = 108, size = 8) +
  geom_line(aes(x = minmaxsequence_raw,
                y = dnorm(minmaxsequence_raw, mean = fp_rma.fit_raw$b[1], sd = sqrt(fp_rma.fit_raw$tau2))),
            linewidth = 1, colour = "darkgrey") +
  geom_line(aes(x = minmaxsequence_corr,
                y = dnorm(minmaxsequence_corr, mean = fp_rma.fit_corr$b[1], sd = sqrt(fp_rma.fit_corr$tau2))),
            linewidth = 1, colour = "black") +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white")
  ) 


p


ggsave(filename = here("Graphics/density_plot.png"),
       plot = last_plot(),
       width = 10,
       height = 5.5)
