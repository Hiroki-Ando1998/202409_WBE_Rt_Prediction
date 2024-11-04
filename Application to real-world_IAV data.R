setwd("C:/wastewater_reproduction_number")


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
shape <- 3  
rate <- 1   
prob <- NULL
for(i in 0:7){
  probability <- gamma_probability(i, i + 1, shape, rate)
  prob <- c(prob, probability) 
}
gamma_prob <- data.frame(prob)

#Sac, 11/27(275), 12/4(282), 12/11 (289)
#Nor, 11/30(278), 12/7(285), 12/14(292), 12/21(299)
#Sou, 11/24(272), 12/1(279), 12/8 (286), 12/15 (293)
#Los, 11/23(271), 11/30(278), 12/7(285)

#Los, 11/16(261), 11/23(268), 11/30(275), 12/7(282), 12/12(287) 
#Nor, 11/23(268), 12/21(296), 12/29(304)
#Sou, 11/17(262), 11/24(269), 12/6(281), 12/22(297)
#Sac, 11/16(261), 11/23(268), 12/21(296)

#Data import and smoothing
data_CA <- read.csv("20240606_IAV_CA_analysis.csv")
data_CA <- data_CA[(1:296),]　#change values of row number
data_stan <- data_CA %>% select(Date, PMMoV_Sac,IAV_Sac) #areas
colnames(data_stan) <- c("Date", "PMMoV", "IAV") #PMMoV: normamlzied, IAV: not normalized

sample_size_1 <- nrow(data_stan)
data_stan_removed <- data_stan %>% drop_na()
WBE <- data_stan_removed$IAV #PMMoV or IAV

data_row <- data.frame(true = which(!is.na(data_stan$IAV)))##PMMoV or IAV
sample_size_2 <- nrow(data_row)
nrow(data_stan_removed)
data_list_ww <- list(n1 = sample_size_1, n2 = sample_size_2, wbe = WBE, row = data_row$true)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mcmc_4 <- stan(
  file = "20240505_stat_space_IAV.stan",
  data = data_list_ww,
  seed = 1,
  chain = 4,
  iter = 40000,
  warmup = 5000,
  thin = 1
)

mcmc_sample <- rstan::extract(mcmc_4)
state_name <- "mu" 
data_ef <- data.frame(mcmc_sample["mu"])
result_ca_1 <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_ca_1) <- c("lower", "median", "upper")
data_ww_plot <- cbind(data_stan, result_ca_1)
data_stan <- data_ww_plot

#Rt calculation
data_ef_2 <- 10^(data_ef)
size <- nrow(data_stan) 
size_2 <- nrow(data_ef_2)
size_3 <- size - 7

Rt_estimate_3 <- data.frame(A = rep(0, (40000-5000)*4)) #the number of mcmc sample
for (i in 9:size){
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
    
    Rt_estimate <- IAV_0/(IAV_1*gamma_prob$prob[1] + IAV_2*gamma_prob$prob[2] + IAV_3*gamma_prob$prob[3] +
                            IAV_4*gamma_prob$prob[4] + IAV_5*gamma_prob$prob[5] + IAV_6*gamma_prob$prob[6] + IAV_7*gamma_prob$prob[7] +
                            IAV_8*gamma_prob$prob[8])
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
data_Date <- data_stan[9:nrow(data_stan),] %>% select(Date)
data_eff_fig <- cbind(data_Date, Rt_result_1)



#data_stan <- read.csv("20240626_1222_South_IAV_95%.csv")#
data_Rt <- Rt_result_1
data_Rt_mcmc <- Rt_estimate_4

data_ef_2 <- 10^(data_ef)
data_eff_fig <- data_Rt
Rt_estimate_4 <- data_Rt_mcmc
print(data_eff_fig)



#Define initial point and initial value of Ws
int_p <- nrow(data_stan) -2
WW_Rt <- data_eff_fig
tail(data_eff_fig, 60) # confirm maximum value of Rt

sp <- which.max(data_eff_fig$median)+8 # row number of maximum value of Rt
ww <- data_stan %>% select(Date, median) %>% mutate(median = 10^(median))
ww <- ww[sp:int_p,]
Ws_0 <- sum(ww$median) 
Ws_1 <- ww[1,2] #highest value day ago
Rt_0 <- WW_Rt[int_p-8, 2] 
Rt_1 <- WW_Rt[sp-8, 2]#1 day ago
W0 <- (-Ws_0*Rt_1)/(Rt_0-Rt_1)
Ws_ini <- W0 - Ws_0

#Ws_0 <- sum(ww$median) 
#Ws_1 <- sum(ww$median) - ww[int_p-sp,2] #1 day ago
#Rt_0 <- WW_Rt[int_p-8, 2] 
#Rt_1 <- WW_Rt[int_p-9, 2]#1 day ago
#W0 <- (Ws_1*Rt_0-Ws_0*Rt_1)/(Rt_1-Rt_0)
#Ws_ini <- W0 - Ws_1

#Prediction
sample_size_3 <- nrow(data_ef_2)/10
data_ww_2 <- data.frame(A = rep(0, 29)) 
data_WS_2 <- data.frame(A = rep(0, 29)) 
data_Rt_2 <- data.frame(A = rep(0, 29)) 

m1 <- int_p
m1_1 <- m1 - 8
m2 <- m1 + 28 
for(n in 1:sample_size_3){
  Rt_estimate <- Rt_estimate_4[10*n, m1_1+1] #231
  WS_estimate <- Ws_ini
  data_ef_3 <- data.frame(data_ef_2[10*n, (1:m1+1)])
  p_ww_estimate_2 <- NULL
  p_WS_estimate_2 <- NULL
  p_Rt_estimate_2 <- NULL
  for(k in m1:m2){
    A <- WS_estimate
    ww_estimate <- Rt_estimate*(data_ef_3[(k-1)]*gamma_prob$prob[1] + data_ef_3[k-2]*gamma_prob$prob[2] + data_ef_3[k-3]*gamma_prob$prob[3] +
                                  data_ef_3[k-4]*gamma_prob$prob[4] + data_ef_3[k-5]*gamma_prob$prob[5] + data_ef_3[k-5]*gamma_prob$prob[6] + 
                                  data_ef_3[k-7]*gamma_prob$prob[7] + data_ef_3[k-8]*gamma_prob$prob[8])
    WS_estimate <- WS_estimate - ww_estimate
    Rt_estimate <- WS_estimate/A*Rt_estimate
    
    data_ef_3 <- cbind(data_ef_3, data.frame(ww_estimate))
    p_ww_estimate_2 <- c(p_ww_estimate_2, ww_estimate)
    p_WS_estimate_2 <- c(p_WS_estimate_2, WS_estimate)
    p_Rt_estimate_2 <- c(p_Rt_estimate_2, Rt_estimate)
  }
  data_ww <- data.frame(p_ww_estimate_2)
  data_ww <- t(data_ww)
  data_WS <- data.frame(p_WS_estimate_2)
  data_WS <- t(data_WS)
  data_Rt <- data.frame(p_Rt_estimate_2)
  data_Rt <- t(data_Rt)
  
  data_ww_2 <- cbind(data_ww_2, data_ww)
  data_WS_2 <- cbind(data_WS_2, data_WS)
  data_Rt_2 <- cbind(data_Rt_2, data_Rt)
}

data_ww_2 <- data_ww_2[,2:sample_size_3]
data_WS_2 <- data_WS_2[,2:sample_size_3]
data_Rt_2 <- data_Rt_2[,2:sample_size_3]


#A_1 <- int_p + 1
#A_2 <- A_1 + 28
#data_time <- data_sim_final %>% select(time, wastewater) 
#data_time <- data_time[(A_1:A_2),]

data_ww_3 <- apply(
  X = data_ww_2, 
  MARGIN = 1,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975))

data_ww_3 <- as.data.frame(t(data_ww_3))
colnames(data_ww_3) <- c("ww_low", "ww_med", "ww_upr")
#data_ww_3 <- cbind(data_time, data_ww_3)

data_ws_3 <- apply(
  X = data_WS_2, 
  MARGIN = 1,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975))  

data_ws_3 <- as.data.frame(t(data_ws_3))
colnames(data_ws_3) <- c("ws_low", "ws_med", "ws_upr")


data_Rt_3 <- apply(
  X = data_Rt_2, 
  MARGIN = 1,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975))

data_Rt_3 <- as.data.frame(t(data_Rt_3))
colnames(data_Rt_3) <- c("Rt_low", "Rt_med", "Rt_upr")

data_result_7 <- cbind(data_ww_3, data_ws_3, data_Rt_3)

write.csv(x = data_result_7, file = "/home/u31/hirokiando/20240629_Prediction_1221_Sac_IAV.csv")
write.csv(x = data_ww_plot, file = "/home/u31/hirokiando/20240629_WWcon95%_1221_Sac_IAV.csv")






#Prediction_2 if Rtww is increasing

#Prediction (data_ef #204)
sample_size_3 <- nrow(data_ef_2)/10
data_ww_2 <- data.frame(A = rep(0, 29))
data_WS_2 <- data.frame(A = rep(0, 29))
data_Rt_2 <- data.frame(A = rep(0, 29)) 

m1 <- int_p
m1_1 <- m1 - 8
m2 <- m1 + 28 
for(n in 1:sample_size_3){
  Rt_estimate <- Rt_estimate_4[n, m1_1] #231
  data_ef_3 <- data.frame(data_ef_2[10*n, (1:m1)])
  p_ww_estimate_2 <- NULL
  for(k in m1:m2){
    ww_estimate <- Rt_estimate*(data_ef_3[(k-1)]*gamma_prob$prob[1] + data_ef_3[k-2]*gamma_prob$prob[2] + data_ef_3[k-3]*gamma_prob$prob[3] +
                                  data_ef_3[k-4]*gamma_prob$prob[4] + data_ef_3[k-5]*gamma_prob$prob[5] + data_ef_3[k-5]*gamma_prob$prob[6] + 
                                  data_ef_3[k-7]*gamma_prob$prob[7] + data_ef_3[k-8]*gamma_prob$prob[8])
    data_ef_3 <- cbind(data_ef_3, data.frame(ww_estimate))
    p_ww_estimate_2 <- c(p_ww_estimate_2, ww_estimate)
  }
  data_ww <- data.frame(p_ww_estimate_2)
  data_ww <- t(data_ww)
  
  data_ww_2 <- cbind(data_ww_2, data_ww)
}

data_ww_3 <- apply(
  X = data_ww_2, 
  MARGIN = 1,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975))

data_ww_3 <- as.data.frame(t(data_ww_3))
colnames(data_ww_3) <- c("ww_low", "ww_med", "ww_upr")
