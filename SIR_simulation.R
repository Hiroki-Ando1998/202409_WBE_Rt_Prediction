setwd("C:/wastewater_reproduction_number")


library(deSolve)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
library(rstan)
library(stats)

# SIRモデル関数を定義する(βの値が途中で変わる)
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # tが100未満の場合、100から150の間、150以上の場合でbetaの値を変更する
    if (time < 45) {
      beta_value <- beta
    } else if (time >= 35 & time < 50) {
      beta_value <- beta/3
    } else {
      beta_value <- beta/1
    }
    dS <- -beta_value * S * I
    dI <- beta_value * S * I - gamma * I
    dR <- gamma * I
    return(list(c(dS, dI, dR)))
  })
}

# 初期値
A <- 0.00001
initial_state <- c(S = 1 - A, I = A, R = 0)

# パラメータ
parameters <- c(beta = 0.50, gamma = 0.25) 
times <- seq(0, 250, by = 1)

# モデルを解く
output <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

data_sim <- data.frame(output)


#新規感染者数に変換
population <- 100000
incidence <- NULL
incidence <- A
run <- nrow(data_sim)

  for(i in 1:run-1){
    n_incidence <- data_sim$S[i] - data_sim$S[i+1]
  incidence <- c(incidence, n_incidence) 
}

data_sim_incidence <- data.frame(incidence)
data_sim_2 <- cbind(data_sim, incidence)
data_sim_2 <- data_sim_2 %>% mutate(incidence_2 = incidence*population)
data_sim_2 <- data_sim_2 %>% mutate(incidence_2 = round(incidence_2))

zeros <- rep(0, 25)
data_b <- data.frame(time = zeros, S = zeros, I = zeros, R = zeros, incidence = zeros, incidence_2 = zeros)
data_c <- rbind(data_b, data_sim_2)
time <- data.frame(Time_2 = seq(0, nrow(data_c)-1, by = 1))
data_sim_3 <- cbind(data_c, time) 


#下水濃度に変換
data_shedding <- read.csv("Assumed_shedding.csv")
wastewater <- numeric()
run <- nrow(data_sim_3)-1
for(i in 26:run){
  n_wastewater_2 <-0
  for(m in 1:25){
  n_wastewater_1 <- data_sim_3$incidence_2[i-m]*data_shedding$super_rapid[m] #rapid, middle, prolong, super_rapid
  n_wastewater_2 <- n_wastewater_2 + n_wastewater_1
  }
  wastewater <- c(wastewater, n_wastewater_2) 
}
data_waste <- data.frame(wastewater)

data_sim_incidence <- data.frame(incidence)
data_sim_2 <- cbind(data_sim, incidence)
data_sim_2 <- data_sim_2 %>% mutate(incidence_2 = incidence*population)


#測定誤差を組み込む
measured_ww <- c(0)
run_2 <- nrow(data_waste)
for(i in 2:run_2){
    measured_2 <- rnorm(n = 1, mean = log10(data_waste$wastewater[i]), sd = 0.5)
  measured_ww <- c(measured_ww, measured_2) 
}

data_measured_ww <- data.frame(measured_ww)
data_waste_2 <- cbind(data_waste, data_measured_ww)

#データの統合
zeros_1 <- rep(0,23)
zeros_1_b <- rep(0,25)
zeros_2 <- rep(0,1)
data_b_2 <- data.frame(wastewater = zeros_1_b)
data_b_2_ww <- data.frame(measured_ww = zeros_1_b)
data_b_2_2 <- cbind(data_b_2, data_b_2_ww)

data_b_3 <- data.frame(wastewater = zeros_2)
data_b_3_ww <- data.frame(measured_ww = zeros_2)
data_b_3_2 <- cbind(data_b_3, data_b_3_ww)

data_b_3 <- data.frame(wastewater = zeros_2)
data_waste_3 <- rbind(data_b_2_2, data_waste_2, data_b_3_2)

data_sim_final <- cbind(data_sim_3, data_waste_3)
data_sim_final <- data_sim_final %>% mutate(measured_ww_2 = if_else(measured_ww <= 0, 0, 10^(measured_ww)))
data_sim_final <- data_sim_final[26:250, ]

#報告誤差をいれる
report_case <- numeric()
run_3 <- nrow(data_sim_final)
for(i in 1:run_3){
  reported_2 <- rpois(n = 1, lambda = data_sim_final$incidence_2[i])
  report_case <- c(report_case, reported_2) 
}

data_report_error <- data.frame(report_case)
data_sim_final <- cbind(data_sim_final, data_report_error)



#incidence
plot <- ggplot(data_sim_final, aes(x = time, y = incidence_2))
plot <- plot + geom_col(fill = "lightskyblue3", colour = "lightskyblue2")
plot <- plot + scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50))
plot <- plot + labs(x = "Time (day)", y = "incidence (case)")
plot <- plot + theme_bw()
plot <- plot + theme(
  axis.line = element_line(size = 1.0, lineend = "square"),
  text = element_text(colour ="black", size = 14),
  legend.position = "none",
  axis.ticks = element_line(linewidth = 1.5),
  axis.ticks.length = unit(-2, "mm"))
plot


#wastewater
plot_ww <- ggplot(data_sim_final, aes(x = time, y = measured_ww_2))
plot_ww <- plot_ww + geom_col(fill = "lightskyblue3", colour = "lightskyblue2")
plot_ww <- plot_ww + scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, 50))
plot_ww <- plot_ww + labs(x = "Time (day)", y = "concentration (copies/g)")
plot_ww <- plot_ww + theme_bw()
plot_ww <- plot_ww + theme(
  axis.line = element_line(size = 1.0, lineend = "square"),
  text = element_text(colour ="black", size = 14),
  legend.position = "none",
  axis.ticks = element_line(linewidth = 1.5),
  axis.ticks.length = unit(-2, "mm"))
plot_ww


#data_sim_final_a <- read.csv("20240520_prolong_data.csv")

#データの一部をNAに変換(7の倍数の行をNAに変換する)
data_sim_final_2 <- data_sim_final %>% filter(wastewater > 0)
nrow(data_sim_final_2)

#three day sampling
data_stan_1 <- data_sim_final_2
data_stan_1[data_stan_1$time %% 7 == 1|data_stan_1$time %% 7 == 3|data_stan_1$time %% 7 == 5 |data_stan_1$time %% 7 == 6, ] <- NA

#one day sampling
data_stan_2 <- data_sim_final_2
data_stan_2[data_stan_2$time %% 7 == 1|data_stan_2$time %% 7 == 2
            |data_stan_2$time %% 7 == 3 |data_stan_2$time %% 7 == 4 
            |data_stan_2$time %% 7 == 5 |data_stan_2$time %% 7 == 6,] <- NA


#状態空間モデル
#


data_stan <- data_sim_final_2 
sample_size_1 <- nrow(data_stan)
data_stan_removed <- data_stan %>% drop_na()
WBE <- data_stan_removed$measured_ww #data_stan_removed$measured_WW or log10(data_stan_removed$wastewaterwastewater)

data_row <- data.frame(true = which(!is.na(data_stan$measured_ww)))#ここのdataファイル名と数字を変える
sample_size_2 <- nrow(data_row)
nrow(data_stan_removed)
data_list_ww <- list(n1 = sample_size_1, n2 = sample_size_2, wbe = WBE, row = data_row$true)

mcmc_1 <- stan(
  file = "20240505_stat_space_1.stan",
  data = data_list_ww,
  seed = 1,
  chain = 1,
  iter = 10000,
  warmup = 4000,
  thin = 1
)

#print(mcmc_1, probe = c(0.025, 0.50, 0.975))


mcmc_sample <- rstan::extract(mcmc_1)
state_name <- "mu" #状態の名前：stan fileを参照すること
result_df_1 <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
  )))

colnames(result_df_1) <- c("lw_7", "median_7", "upr_7")


#three day sampling
data_stan <- data_stan_1 

sample_size_1 <- nrow(data_stan)

data_stan_removed <- data_stan %>% drop_na()
WBE <- data_stan_removed$measured_ww #measured_WW or log10(wastewater)

data_row <- data.frame(true = which(!is.na(data_stan$measured_ww)))#ここのdataファイル名と数字を変える
sample_size_2 <- nrow(data_row)

data_list_ww_2 <- list(n1 = sample_size_1, n2 = sample_size_2, wbe = WBE, row = data_row$true)

mcmc_2 <- stan(
  file = "20240505_stat_space_1.stan",
  data = data_list_ww_2,
  seed = 1,
  chain = 1,
  iter = 10000,
  warmup = 4000,
  thin = 1
)
#print(mcmc_2, probe = c(0.025, 0.50, 0.975))

mcmc_sample <- rstan::extract(mcmc_2)
state_name <- "mu" #状態の名前：stan fileを参照すること
result_df_2 <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_df_2) <- c("lw_3", "median_3", "upr_3")


#One day sampling
data_stan <- data_stan_2 

sample_size_1 <- nrow(data_stan)

data_stan_removed <- data_stan %>% drop_na()
WBE <- data_stan_removed$measured_ww #measured_WWW or log10(wastewater)

data_row <- data.frame(true = which(!is.na(data_stan$measured_ww)))#ここのdataファイル名と数字を変える
sample_size_2 <- nrow(data_row)

data_list_ww_3 <- list(n1 = sample_size_1, n2 = sample_size_2, wbe = WBE, row = data_row$true)

mcmc_3 <- stan(
  file = "20240505_stat_space_1.stan",
  data = data_list_ww_3,
  seed = 1,
  chain = 1,
  iter = 10000,
  warmup = 4000,
  thin = 1
)
#print(mcmc_3, probe = c(0.025, 0.50, 0.975))

#library(bayesplot)
#mcmc_combo(mcmc_3, pars = c("mu[1]", "s1", "s2"))



mcmc_sample <- rstan::extract(mcmc_3)
state_name <- "mu" #状態の名前：stan fileを参照すること
result_df_3 <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_df_3) <- c("lw_1", "median_1", "upr_1")



data_m <- cbind(result_df_1, result_df_2, result_df_3)

#write.csv(x = data_m, file = "C:/wastewater_reproduction_number/2022_data.csv")


#時間軸の追加
data_x <- data_sim_final_2 %>% select(time, measured_ww, wastewater)
data_x_2 <- cbind(data_x, data_m)

data_subset <- data_x_2 %>% mutate(wastewater_2 = log10(wastewater))







#Validation figure (summary of all sampling condition)
data_sampling_7 <- data_subset %>% select("time", "lw_7", "median_7", "upr_7")
data_name_7 <- data.frame(c = "7day_sampling")
data_sampling_7 <- cbind(data_sampling_7, data_name_7)
colnames(data_sampling_7) <- c("time", "min", "median", "max", "frequency")

data_sampling_3 <- data_subset %>% select("time", "lw_3", "median_3", "upr_3")
data_name_3 <- data.frame(c = "3day_sampling")
data_sampling_3 <- cbind(data_sampling_3, data_name_3)
colnames(data_sampling_3) <- c("time", "min", "median", "max", "frequency")

data_sampling_1 <- data_subset %>% select("time", "lw_1", "median_1", "upr_1")
data_name_1 <- data.frame(c = "1day_sampling")
data_sampling_1 <- cbind(data_sampling_1, data_name_1)
colnames(data_sampling_1) <- c("time", "min", "median", "max", "frequency")

data_true <- data_subset %>% select("time", "wastewater_2")
data_name <- data.frame(c = "true_value")
data_true_2 <- data.frame(time = data_true$time, min = data_true$wastewater_2,
                          median = data_true$wastewater_2, max = data_true$wastewater_2)
data_true_2 <- cbind(data_true_2, data_name)
colnames(data_true_2) <- c("time", "min", "median", "max", "frequency")

data_sub_file <- rbind(data_true_2, data_sampling_1, data_sampling_3, data_sampling_7)





#Fiugre
plot <- ggplot(data = data_sub_file, aes(x = time, y = median, colour = frequency, fill = frequency)) +
  geom_ribbon(aes(ymin = min, ymax = max),colour = "NA", alpha = 0.5) +
  geom_line(size = 1.2)
plot <- plot + scale_fill_manual(values = c("#D0D1E6", "#A6BDDB", "#74A9CF", "#D7301F"))
plot <- plot + scale_colour_manual(values = c("#8C96C6", "#3690C0", "#0570B0", "#B30000"))
plot <- plot + scale_y_continuous(limits = c(-1,7), breaks = seq(-1, 7, 2))
plot <- plot + scale_x_continuous(limits = c(0,135), breaks = seq(0, 135, 50))
plot <- plot +theme_bw()
plot <- plot +  xlab("Time (day)")
plot <- plot + theme(
  axis.line = element_line(size = 1.0, lineend = "square"),
  text = element_text(colour ="black", size = 14),
  legend.position = "NA",
  axis.ticks = element_line(linewidth = 1.5),
  axis.ticks.length = unit(2, "mm"))
plot

#Calculation of average
data_ave <- data_subset %>% select(time, median_7, median_3, median_1, wastewater_2) #138行目(rapid)まで
#data_ave <- data_subset[1:146, ] #Bust epidemic

#one day sampling
df_1 <- data_ave %>% select(time, median_1, wastewater_2)
data_ave_1 <- df_1 %>% mutate(error = (median_1 - wastewater_2)^2)
data_ave_1 <- data_ave_1 %>% mutate(standard_error = abs(median_1/wastewater_2 - 1))
average_error_1 <- sqrt(sum(data_ave_1$error)/nrow(data_ave_1))
print(average_error_1)
average_error_1_b <- sum(data_ave_1$standard_error)/nrow(data_ave_1)
print(average_error_1_b)


#three day sampling
df_3 <- data_ave %>% select(time, median_3, wastewater_2)
data_ave_3 <- df_3 %>% mutate(error = (median_3 - wastewater_2)^2)
data_ave_3 <- data_ave_3 %>% mutate(standard_error = abs(median_3/wastewater_2 - 1))
average_error_3 <- sqrt(sum(data_ave_3$error)/nrow(data_ave_3))
print(average_error_3)
average_error_3_b <- sum(data_ave_3$standard_error)/nrow(data_ave_3)
print(average_error_3_b)


#everyday day sampling
df_7 <- data_ave %>% select(time, median_7, wastewater_2)
data_ave_7 <- df_7 %>% mutate(error = (median_7 - wastewater_2)^2)
data_ave_7 <- data_ave_7 %>% mutate(standard_error = abs(median_7/wastewater_2 - 1))
average_error_7 <- sqrt(sum(data_ave_7$error)/nrow(data_ave_7))
print(average_error_7)
average_error_7_b <- sum(data_ave_7$standard_error)/nrow(data_ave_7)
print(average_error_7_b)









#Effective reproduction number 
#average: shape / scale
#variance: shape / scale^2
#samples <- rgamma(1000, shape = 4, scale = 1)
#hist(samples, breaks = 30, col = "skyblue", main = "Gamma Distribution")
data_nrow_1 <- data_sim_final %>% filter(wastewater > 0) #(144行目のファイルを毎回確認する)
data_nrow_2 <- nrow(data_nrow_1) + 1
data_sim_final_2 <- data_sim_final[1:data_nrow_2, ]
data_eff_1 <- data_sim_final_2 %>% select(time, wastewater, incidence_2, report_case, measured_ww)

zeros_2 <- rep(0,1)
data_c_0 <- data.frame(lw_7 = zeros_2)
data_c_1 <- data.frame(median_7 = zeros_2)
data_c_1_a <- data.frame(upr_7 = zeros_2)
data_c_2 <- data.frame(lw_3 = zeros_2)
data_c_3 <- data.frame(median_3 = zeros_2)
data_c_3_b <- data.frame(upr_3 = zeros_2)
data_c_4 <- data.frame(lw_1 = zeros_2)
data_c_5 <- data.frame(median_1 = zeros_2)
data_c_6 <- data.frame(upr_1 = zeros_2)

data_c_2_2 <- cbind(data_c_0, data_c_1, data_c_1_a, data_c_2, data_c_3, data_c_3_b, data_c_4, data_c_5, data_c_6)

data_m_1 <- rbind(data_c_2_2, data_m)


data_eff_1 <- cbind(data_eff_1, data_m_1)
data_eff_1 <- data_eff_1 %>% select(time, wastewater, incidence_2, report_case, median_7, median_3, median_1, measured_ww)

data_eff <- data_eff_1 %>% select(time, incidence_2, report_case, median_7, median_3, median_1, measured_ww, wastewater)
colnames(data_eff) <- c("time", "incidence_2", "report_case", "ww_7", "ww_3", "ww_1", "measured_7", "ww_true" )

data_eff <- data_eff %>% mutate(ww_7 = 10^(ww_7), ww_3 = 10^(ww_3),ww_1 = 10^(ww_1), measured_7 = 10^(measured_7))
data_eff <- data_eff %>% select(time, incidence_2, report_case, ww_true, measured_7, ww_7, ww_3, ww_1)

#諸々データファイルを調整
data_a_1 <- data_eff[1:1, ] %>% select(time, incidence_2, report_case)
zeros_2 <- rep(0,1)
data_d_0 <- data.frame(ww_true = zeros_2)
data_d_1 <- data.frame(measured_7 = zeros_2)
data_d_1_a <- data.frame(ww_7 = zeros_2)
data_d_2 <- data.frame(ww_3 = zeros_2)
data_d_3 <- data.frame(ww_1 = zeros_2)  

data_d_2_2 <- cbind(data_a_1, data_d_0, data_d_1, data_d_1_a, data_d_2, data_d_3)

data_a_2 <- data_eff[2:nrow(data_eff), ]

#解析に使うファイルがほぼ完成
data_eff <- rbind(data_d_2_2, data_a_2)


zeros_1 <- rep(0,25)
data_b_0 <- data.frame(time = zeros_1)
data_b_1 <- data.frame(incidence_2 = zeros_1)
data_b_1_a <- data.frame(report_case = zeros_1)
data_b_2 <- data.frame(ww_true = zeros_1)
data_b_3 <- data.frame(measured_7 = zeros_1)
data_b_3_b <- data.frame(ww_7 = zeros_1)
data_b_4 <- data.frame(ww_3 = zeros_1)
data_b_5 <- data.frame(ww_1 = zeros_1)

data_b_2_2 <- cbind(data_b_0, data_b_1, data_b_1_a, data_b_2, data_b_3, data_b_3_b, data_b_4, data_b_5)

data_eff_2 <- rbind(data_b_2_2, data_eff)




#clinical_Re
data_stan_eff <- data_eff_2 
sample_size_1 <- nrow(data_stan_eff)

#generation distribution

# ガンマ関数の定義
gamma_function <- function(x, shape, rate) {
  return((rate^shape) * (x^(shape - 1)) * exp(-rate * x) / gamma(shape))
}

# 確率密度関数の積分を計算する関数
gamma_probability <- function(lower, upper, shape, rate) {
  integrate(gamma_function, lower = lower, upper = upper, shape = shape, rate = rate)$value
}

shape <- 4  # ガンマ関数の形状パラメータ
rate <- 1   # ガンマ関数のレートパラメータ

prob <- NULL
for(i in 0:15){
  probability <- gamma_probability(i, i + 1, shape, rate)
  prob <- c(prob, probability) 
}
generation <- data.frame(generation = prob)

incidence <- data_eff$incidence_2
data_list_eff_1 <- list(n1 = sample_size_1, n2 = nrow(generation)-1, incidence = incidence, generation = generation$generation)



#data_expect_true
case_all <- numeric()
run <- nrow(data_eff_2)
for(i in 25:run){
  case_2 <-0
  for(m in 1:nrow(generation)){
    case_1 <- data_eff_2$incidence_2[i-m]*generation$generation[m] #incidence_2, report_case
    case_2 <- case_2 + case_1
  }
  case_all <- c(case_all, case_2) 
}
data_expected_case <- data.frame(case_all)

zeros_1 <- rep(0,24)
data_b_0 <- data.frame(case_all = zeros_1)
data_eff_2_b <- rbind(data_b_0, data_expected_case)
data_eff_3 <- cbind(data_eff_2, data_eff_2_b)
data_eff_3 <- data_eff_3 %>% filter(case_all > 0)


#incidence_2, data_report_error
data_list_eff_2 <- list(n1 = nrow(data_eff_3), incidence = data_eff_3$incidence_2, expect = data_eff_3$case_all)

mcmc_eff_cli_incidence <- stan(
  file = "20240515_clinical_Re_2.stan",
  data = data_list_eff_2,
  seed = 1,
  chain = 1,
  iter = 10000,
  warmup = 4000,
  thin = 1
)
print(mcmc_eff_cli_incidence, probe = c(0.025, 0.50, 0.975))

mcmc_sample <- rstan::extract(mcmc_eff_cli_incidence)
state_name <- "Re" #状態の名前：stan fileを参照すること
result_eff_cli_1_incidence <- data.frame(t(apply(
  X = mcmc_sample[[state_name]],
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.5, 0.975)
)))

colnames(result_eff_cli_1_incidence) <- c("low_incidenc_1", "med_incidenc_1", "upr_incidenc_1")


#write.csv(x = result_eff_cli_1_incidence, file = "C:/wastewater_reproduction_number/20240520_1_7day.csv")


#expected wastewater concentration
#ww_1
waste_all <- numeric()
run <- nrow(data_eff_2)
for(i in 25:run){
  waste_2 <-0
  for(m in 1:nrow(generation)){
    waste_1 <- data_eff_2$ww_1[i-m]*generation$generation[m] #median_7, ww_3, ww_true
    waste_2 <- waste_2 + waste_1
  }
  waste_all <- c(waste_all, waste_2) 
}
data_expected_waste <- data.frame(waste_all)

zeros_1 <- rep(0,24)
data_b_0 <- data.frame(waste_all = zeros_1)
data_eff_2_waste <- rbind(data_b_0, data_expected_waste)
data_eff_waste <- cbind(data_eff_2, data_eff_2_waste)
data_eff_waste <- data_eff_waste %>% filter(waste_all > 0 & incidence_2 > 0)


#log10(data_eff_waste$ww_3)
data_list_eff_1 <- list(n1 = nrow(data_eff_waste), waste = log10(data_eff_waste$ww_1), expect = data_eff_waste$waste_all)

mcmc_eff_ww_1 <- stan(
  file = "20240515_wastewater_Re_2.stan",
  data = data_list_eff_1,
  seed = 1,
  chain = 1,
  iter = 10000,
  warmup = 4000,
  thin = 1
)
#print(mcmc_eff_ww_1, probe = c(0.025, 0.50, 0.975))
  
  mcmc_sample_eff_ww_1 <- rstan::extract(mcmc_eff_ww_1)
  state_name <- "Re" #状態の名前：stan fileを参照すること
  result_eff_ww_1 <- data.frame(t(apply(
    X = mcmc_sample_eff_ww_1[[state_name]],
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.025, 0.5, 0.975)
  )))
  
  colnames(result_eff_ww_1) <- c("low_ww_1", "med_1", "upr_1")
  
  
  
  
  #ww_3
  waste_all <- numeric()
  run <- nrow(data_eff_2)
  for(i in 25:run){
    waste_2 <-0
    for(m in 1:nrow(generation)){
      waste_1 <- data_eff_2$ww_3[i-m]*generation$generation[m] #median_7, ww_3, estimated
      waste_2 <- waste_2 + waste_1
    }
    waste_all <- c(waste_all, waste_2) 
  }
  data_expected_waste <- data.frame(waste_all)
  
  zeros_1 <- rep(0,24)
  data_b_0 <- data.frame(waste_all = zeros_1)
  data_eff_2_waste <- rbind(data_b_0, data_expected_waste)
  data_eff_waste <- cbind(data_eff_2, data_eff_2_waste)
  data_eff_waste <- data_eff_waste %>% filter(incidence_2 > 0 & waste_all > 0)
  
  data_1 <- data_eff_waste %>% mutate(Re = ww_3/waste_all) #median_7, wastewater, ww_3
  
  #waste: log10(estimated_7), log10(ww_7) or log10(ww_3) or log10(ww_1) or log10(wastewater)
  data_list_eff_3 <- list(n1 = nrow(data_eff_waste), waste = log10(data_eff_waste$ww_3), expect = data_eff_waste$waste_all)
  
  mcmc_eff_ww_3 <- stan(
    file = "20240515_wastewater_Re_2.stan",
    data = data_list_eff_3,
    seed = 1,
    chain = 1,
    iter = 10000,
    warmup = 4000,
    thin = 1
  )
  #print(mcmc_eff_ww, probe = c(0.025, 0.50, 0.975))
  
  mcmc_sample_eff_ww_3 <- rstan::extract(mcmc_eff_ww_3)
  state_name <- "Re" #状態の名前：stan fileを参照すること
  result_eff_ww_3 <- data.frame(t(apply(
    X = mcmc_sample_eff_ww_3[[state_name]],
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.025, 0.5, 0.975)
  )))
  
  colnames(result_eff_ww_3) <- c("low_ww_3", "med_ww_3", "upr_ww_3")
  
  
  
  #ww_7
  waste_all <- numeric()
  run <- nrow(data_eff_2)
  for(i in 25:run){
    waste_2 <-0
    for(m in 1:nrow(generation)){
      waste_1 <- data_eff_2$ww_7[i-m]*generation$generation[m] #median_7, ww_3, estimated
      waste_2 <- waste_2 + waste_1
    }
    waste_all <- c(waste_all, waste_2) 
  }
  data_expected_waste <- data.frame(waste_all)
  
  zeros_1 <- rep(0,24)
  data_b_0 <- data.frame(waste_all = zeros_1)
  data_eff_2_waste <- rbind(data_b_0, data_expected_waste)
  data_eff_waste <- cbind(data_eff_2, data_eff_2_waste)
  data_eff_waste <- data_eff_waste %>% filter(incidence_2 > 0 & waste_all > 0)
  
  data_1 <- data_eff_waste %>% mutate(Re = ww_7/waste_all) #median_7, wastewater, ww_3
  
  #waste: log10(estimated_7), log10(ww_7) or log10(ww_3) or log10(ww_1) or log10(wastewater)
  data_list_eff_7 <- list(n1 = nrow(data_eff_waste), waste = log10(data_eff_waste$ww_7), expect = data_eff_waste$waste_all)
  
  mcmc_eff_ww_7 <- stan(
    file = "20240515_wastewater_Re_2.stan",
    data = data_list_eff_7,
    seed = 1,
    chain = 1,
    iter = 10000,
    warmup = 4000,
    thin = 1
  )
  #print(mcmc_eff_ww, probe = c(0.025, 0.50, 0.975))
  
  mcmc_sample_eff_ww_7 <- rstan::extract(mcmc_eff_ww_7)
  state_name <- "Re" #状態の名前：stan fileを参照すること
  result_eff_ww_7 <- data.frame(t(apply(
    X = mcmc_sample_eff_ww_7[[state_name]],
    MARGIN = 2,
    FUN = quantile,
    probs = c(0.025, 0.5, 0.975)
  )))
  
  colnames(result_eff_ww_7) <- c("low_ww_7", "med_ww_7", "upr_ww_7")
  
  
  

  
  
  
  #行の長さを調整
  result_eff_cli_1_incidence_a <- result_eff_cli_1_incidence[1:nrow(result_eff_ww_1) + 1, ]
  #時間の追加
  data_time <- data.frame(time = seq(from = 2, to = nrow(result_eff_ww_1) +1))
  result_eff_ww_1 <-result_eff_ww_1[1:nrow(result_eff_ww_1), ]
  result_eff_ww_3 <-result_eff_ww_3[1:nrow(result_eff_ww_3), ]
  result_eff_ww_7 <-result_eff_ww_7[1:nrow(result_eff_ww_7), ]
  
  data_re_subset <- cbind(data_time, result_eff_cli_1_incidence_a, result_eff_ww_1, result_eff_ww_3, result_eff_ww_7)

  #zeros_1 <- rep(0,1)
  #data_b_0 <- data.frame( = zeros_1)
  #data_eff_2_b <- rbind(data_b_0, data_expected_case)

  
  #calculation(最初の五行は計算から除外する)
  data_cal_re_1 <- data_re_subset[14:nrow(result_eff_ww_1),] #14 or 5
  #one day sampling
  df_re_1 <- data_cal_re_1 %>% select(time, med_incidenc_1, med_1)
  data_ave_re_1 <- df_re_1 %>% mutate(error = (med_incidenc_1 - med_1)^2)
  data_ave_re_1 <- data_ave_re_1 %>% mutate(standard_error = abs(med_1/med_incidenc_1 - 1))
  average_error_re_1 <- sqrt(sum(data_ave_re_1$error)/nrow(data_cal_re_1))
  print(average_error_1)
  average_error_1_b <- (sum(data_ave_re_1$standard_error)/nrow(data_cal_re_1))
  print(average_error_1_b)
  
  
  #three day sampling
  df_re_3 <- data_cal_re_1 %>% select(time, med_incidenc_1, med_ww_3)
  data_ave_re_3 <- df_re_3 %>% mutate(error = (med_incidenc_1 - med_ww_3)^2)
  data_ave_re_3 <- data_ave_re_3 %>% mutate(standard_error = abs(med_ww_3/med_incidenc_1 - 1))
  average_error_re_3 <- sqrt(sum(data_ave_re_3$error)/nrow(data_cal_re_1))
  print(average_error_3)
  average_error_3_b <- (sum(data_ave_re_3$standard_error)/nrow(data_cal_re_1))
  print(average_error_3_b)
  
  
  #everyday day sampling
  df_re_7 <- data_cal_re_1 %>% select(time, med_incidenc_1, med_ww_7)
  data_ave_re_7 <- df_re_7 %>% mutate(error = (med_incidenc_1 - med_ww_7)^2)
  data_ave_re_7 <- data_ave_re_7 %>% mutate(standard_error = abs(med_ww_7/med_incidenc_1 - 1))
  average_error_re_7 <- sqrt(sum(data_ave_re_7$error)/nrow(data_cal_re_1))
  print(average_error_7)
  average_error_7_b <- (sum(data_ave_re_7$standard_error)/nrow(data_cal_re_1))
  print(average_error_7_b)
  
  
  
  
  
  
  
    
  #Validation figure (summary of all sampling condition)
  data_re_1_incidence <- data_re_subset %>% select("time", "low_incidenc_1", "med_incidenc_1", "upr_incidenc_1")
  data_name_re_1_inc <- data.frame(c = "1day_cli_incidence")
  data_re_1_incidence <- cbind(data_re_1_incidence, data_name_re_1_inc)
  colnames(data_re_1_incidence) <- c("time", "min", "median", "max", "Re_type")  
  
  data_re_ww_1 <- data_re_subset %>% select("time", "low_ww_1", "med_1", "upr_1")
  data_name_ww_1 <- data.frame(c = "1day_ww")
  data_re_ww_1 <- cbind(data_re_ww_1, data_name_ww_1)
  colnames(data_re_ww_1) <- c("time", "min", "median", "max", "Re_type")

  data_re_ww_3 <- data_re_subset %>% select("time", "low_ww_3", "med_ww_3", "upr_ww_3")
  data_name_ww_3 <- data.frame(c = "3day_ww")
  data_re_ww_3 <- cbind(data_re_ww_3, data_name_ww_3)
  colnames(data_re_ww_3) <- c("time", "min", "median", "max", "Re_type")
  
  data_re_ww_7 <- data_re_subset %>% select("time", "low_ww_7", "med_ww_7", "upr_ww_7")
  data_name_ww_7 <- data.frame(c = "7day_ww")
  data_re_ww_7 <- cbind(data_re_ww_7, data_name_ww_7)
  colnames(data_re_ww_7) <- c("time", "min", "median", "max", "Re_type")
    
  
  data_sub_re_file <- rbind(data_re_1_incidence, data_re_ww_1, data_re_ww_3, data_re_ww_7)
  
  
  #Fiugre
  plot <- ggplot(data = data_sub_re_file, aes(x = time, y = median, colour = Re_type, fill = Re_type)) +
    geom_ribbon(aes(ymin = min, ymax = max),colour = "NA", alpha = 0.5) +
    geom_line(size = 1.2)
  plot <- plot + scale_fill_manual(values = c( "#FC8D59","#D0D1E6", "#A6BDDB", "#74A9CF"))
  plot <- plot + scale_colour_manual(values = c( "#D7301F", "#8C96C6", "#3690C0", "#0570B0"))
  plot <- plot + scale_y_continuous(limits = c(0,5), breaks = seq(0, 5, 1))
  plot <- plot + scale_x_continuous(limits = c(0,130), breaks = seq(0, 130, 50))
  plot <- plot +theme_bw()
  plot <- plot +  xlab("Time (day)") + ylab("Wastewater concetration (Log10 copies/L)")
  plot <- plot + theme(
    axis.line = element_line(size = 1.0, lineend = "square"),
    text = element_text(colour ="black", size = 14),
    legend.position = "NA",
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(2, "mm"))
  plot
  


