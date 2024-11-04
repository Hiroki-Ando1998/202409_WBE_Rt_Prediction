
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(rstan)
library(stats)


gamma_function <- function(x, shape, rate) {
  return((rate^shape) * (x^(shape - 1)) * exp(-rate * x) / gamma(shape))
}
gamma_probability <- function(lower, upper, shape, rate) {
  integrate(gamma_function, lower = lower, upper = upper, shape = shape, rate = rate)$value
}
shape <- 4  
rate <- 1   
prob <- NULL
for(i in 0:9){
  probability <- gamma_probability(i, i + 1, shape, rate)
  prob <- c(prob, probability) 
}
gamma_prob <- data.frame(prob)



#Data import and smoothing
data_CA <- read.csv("20241025_SARS-CoV-2_COVID-19case.csv")

data_stan <- data_CA

sample_size_1 <- nrow(data_stan)
data_stan_removed <- data_stan %>% drop_na()



sample_size_1 <- nrow(data_stan)
data_stan_removed <- data_stan %>% drop_na()
incidence <- data_stan_removed$case

data_row <- data.frame(true = which(!is.na(data_stan$case)))
sample_size_2 <- nrow(data_row)
nrow(data_stan_removed)
data_list_ww <- list(n1 = sample_size_1, incidence = incidence)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mcmc_4 <- stan(
  file = "20241025_COVID-19case.stan",
  data = data_list_ww,
  seed = 1,
  chain = 4,
  iter = 40000,
  warmup = 5000,
  thin = 1
)

mcmc_sample <- rstan::extract(mcmc_4)
state_name <- "Re" 
data_ef <- data.frame(mcmc_sample["Re"])
result_ca_Re <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_ca_Re) <- c("lower", "median", "upper")
#write.csv(x = result_ca_Re, file = "/home/u31/hirokiando/result_Re_CA_SARS.csv")






WBE <- data_stan_removed$Log 
data_row <- data.frame(true = which(!is.na(data_stan$Log)))
sample_size_2 <- nrow(data_row)
nrow(data_stan_removed)
data_list_ww <- list(n1 = sample_size_1, n2 = sample_size_2, wbe = WBE, row = data_row$true)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mcmc_3 <- stan(
  file = "20240505_stat_space_IAV.stan",
  data = data_list_ww,
  seed = 1,
  chain = 4,
  iter = 40000,
  warmup = 5000,
  thin = 1
)

mcmc_sample <- rstan::extract(mcmc_3)
state_name <- "mu" 
data_ef <- data.frame(mcmc_sample["mu"])
result_ca_1 <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_ca_1) <- c("lower", "median", "upper")
#write.csv(x = result_ca_1, file = "/home/u31/hirokiando/result_wastewaterconcentration_CA_SARS.csv")


data_ww_plot <- cbind(data_stan, result_ca_1)
data_stan <- data_ww_plot

#Rt calculation
data_ef_2 <- 10^(data_ef)
size <- nrow(data_stan) 
size_2 <- nrow(data_ef_2)
size_3 <- size - 9

Rt_estimate_3 <- data.frame(A = rep(0, (40000-5000)*4)) #the number of mcmc sample
for (i in 11:size){
  Rt_estimate_2 <- NULL   
  for(m in 1:size_2){
    IAV_0 <- data_ef_2[m,i] #m行目i列目
    IAV_1 <- data_ef_2[m,i-1]
    IAV_2 <- data_ef_2[m,i-2]
    IAV_3 <- data_ef_2[m,i-3]
    IAV_4 <- data_ef_2[m,i-4]
    IAV_5 <- data_ef_2[m,i-5]
    IAV_6 <- data_ef_2[m,i-6]
    IAV_7 <- data_ef_2[m,i-7]
    IAV_8 <- data_ef_2[m,i-8]
    IAV_9 <- data_ef_2[m,i-9]
    IAV_10 <- data_ef_2[m,i-10]
    
    Rt_estimate <- IAV_0/(IAV_1*gamma_prob$prob[1] + IAV_2*gamma_prob$prob[2] + IAV_3*gamma_prob$prob[3] +
                            IAV_4*gamma_prob$prob[4] + IAV_5*gamma_prob$prob[5] + IAV_6*gamma_prob$prob[6] + IAV_7*gamma_prob$prob[7] +
                            IAV_8*gamma_prob$prob[8]+IAV_9*gamma_prob$prob[9]+IAV_10*gamma_prob$prob[10])
    Rt_estimate_2 <- c(Rt_estimate_2, Rt_estimate)
  }
  Rt_estimate_2 <- data.frame(Rt_estimate_2)
  Rt_estimate_3 <- cbind(Rt_estimate_3, Rt_estimate_2)
}
Rt_estimate_4 <- Rt_estimate_3[,2:size_3]
Rt_result_1 <- data.frame(t(apply(
  X = Rt_estimate_4,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(Rt_result_1) <- c("lw", "median", "upr")
data_Date <- data_stan[11:nrow(data_stan),] %>% select(Date)
data_eff_fig <- cbind(data_Date, Rt_result_1)

#write.csv(x = Rt_result_1, file = "/home/u31/hirokiando/result_Reww_CA_SARS.csv")


