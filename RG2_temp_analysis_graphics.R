
packages <- c("magrittr", "dplyr", "here", "metafor")

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


ggplot(B_alpha_rma_df) +
  geom_violin(aes(x = 1, y = sqrt(tau2_alpha))) +
  geom_boxplot(aes(x = 1, y = sqrt(tau2_alpha))) +
  geom_point(aes(x = 1, y = sqrt(tau2_alpha)), position = position_jitter(width = .1))


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


