library(tidyverse)
library(here)
options(scipen = 99999)


sim_exp <- data.frame()
for(i in 1:5){
  name <- paste0("simulations_exp_",i,".rds")
  sim_exp <- rbind(sim_exp, readRDS(here("SimulationData",name)))
}

sim_lnormal <- data.frame()
for(i in 1:5){
  name <- paste0("simulations_lnormal_",i,".rds")
  sim_lnormal <- rbind(sim_lnormal, readRDS(here("SimulationData",name)))
}

sim_normal <- data.frame()
for(i in 1:5){
  name <- paste0("simulations_normal_",i,".rds")
  sim_normal <- rbind(sim_normal, readRDS(here("SimulationData",name)))
}

all_sim <- rbind(sim_exp,sim_lnormal, sim_normal) %>%
  mutate(Z_p.value= 2*pnorm(abs(S_statistic), lower.tail = FALSE),
         t_p.value = 2*pt(abs(S_statistic), lower.tail = FALSE,df = (number_individuals * number_measurements) - 3),
         KR_sig = ifelse(KR_p.value < .05,TRUE,FALSE),
         S_sig = ifelse(S_p.value < .05,TRUE,FALSE),
         Z_sig = ifelse(Z_p.value < .05,TRUE,FALSE),
         t_sig = ifelse(t_p.value < .05,TRUE,FALSE))

#saveRDS(all_sim,"all_sim.rds")


all_sim_sum <- all_sim %>%
  group_by(distribution,
                    number_individuals,
                    params,
                    number_measurements,
                    rslope,
                    KR_effect, KR_term, skew, kurtosis) %>%
  summarize(
            KR_t1err = mean(KR_sig),
            S_t1err = mean(S_sig),
            Z_t1err = mean(Z_sig),
            t_t1err = mean(t_sig))

all_sim_sum %>%
  filter(KR_term == "treatment") %>%
  pivot_wider(names_from = c("number_individuals","rslope"),
              values_from = "KR_t1err")%>% group_by(distribution, params,number_measurements) %>%
  summarise_all(funs(.[!is.na(.)])) %>%
  slice(1) %>%
  ungroup() %>%
  select(-c(KR_effect,KR_term,S_t1err,Z_t1err,t_t1err, distribution)) %>%
  kbl() %>%
  pack_rows("Exponential", 1, 4) %>%
  pack_rows("Normal", 5, 6) %>%
  pack_rows("Lognormal", 7, 10) %>%
  add_header_above(c(" " = 2, "Random Intercept" = 1, "Random Slope" = 1,"Random Intercept" = 1,
                     "Random Slope" = 1,"Random Intercept" = 1,"Random Slope" = 1)) %>%
  add_header_above(c(" " = 2, "10" = 2, "18" = 2, "26" = 2)) %>%
  add_header_above(c(" " = 2, "Sample Size" = 6))  


all_sim_sum %>%
  filter(KR_term == "treatment") %>%
  pivot_longer(cols = c("KR_t1err","S_t1err","Z_t1err","t_t1err"),names_to = "DF_method", values_to = "error_rate") %>%
  filter(DF_method %in% c("KR_t1err","S_t1err")) %>%
  pivot_wider(names_from = c("rslope"),
              values_from = "error_rate", names_prefix = "slope") %>%
  group_by(number_individuals,skew,kurtosis,DF_method) %>%
  summarize(random_intercept = mean(slopeFALSE), random_slope = mean(slopeTRUE)) %>%
  pivot_wider(names_from = c("number_individuals"),
              values_from = c("random_intercept","random_slope")) %>%
  select(skew,kurtosis,DF_method,random_intercept_10,random_slope_10,random_intercept_18,
         random_slope_18,random_intercept_26,random_slope_26) %>%
  ungroup() %>%
  kbl() %>%
  pack_rows("Exponential", 7, 8) %>%
  pack_rows("Normal", 1, 2) %>%
  pack_rows("Lognormal", 3, 6) %>%
  add_header_above(c(" " = 3, "Random Intercept" = 1, "Random Slope" = 1,"Random Intercept" = 1,
                     "Random Slope" = 1,"Random Intercept" = 1,"Random Slope" = 1)) %>%
  add_header_above(c(" " = 3, "10" = 2, "18" = 2, "26" = 2)) %>%
  add_header_above(c(" " = 3, "Sample Size" = 6))  





# OTHER VISUALIZATIONS
# P- value histograms? # qq plot?
# do i have to do difference testing to test signficance between exp/normal/lnormal
hist(sim_exp$KR_p.value )