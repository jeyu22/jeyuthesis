library(tidyverse)
library(here)

sim_exp <- data.frame()
for(i in 1:10){
  name <- paste0("simulations_exp_",i,".rds")
  sim_exp <- rbind(sim_exp, readRDS(here("SimulationData",name)))
}

sim_lnormal <- data.frame()
for(i in 1:10){
  name <- paste0("simulations_lnormal_",i,".rds")
  sim_lnormal <- rbind(sim_lnormal, readRDS(here("SimulationData",name)))
}

sim_normal <- data.frame()
for(i in 1:10){
  name <- paste0("simulations_normal_",i,".rds")
  sim_normal <- rbind(sim_normal, readRDS(here("SimulationData",name)))
}


# KR VS S: same
sim_exp %>%
  filter(KR_p.value < .05) %>% nrow()
83/3600

sim_exp %>%
  filter(S_p.value < .05) %>% nrow()

sim_lnormal %>%
  filter(KR_p.value < .05) %>% nrow()
117/4800

# Number of individuals: less individuals more error rate
sim_exp %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_individuals) %>%
  count() %>%
  mutate(type1 = n/1200)

sim_lnormal %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_individuals) %>%
  count() %>%
  mutate(type1 = n/1600)

sim_normal %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_individuals) %>%
  count() %>%
  mutate(type1 = n/1200)

# Params 
sim_exp %>%
  filter(KR_p.value < .05) %>% 
  group_by(params) %>%
  count() %>%
  mutate(type1 = n/1200)

sim_lnormal %>%
  filter(KR_p.value < .05) %>% 
  group_by(params) %>%
  count() %>%
  mutate(type1 = n/1200)

sim_normal %>%
  filter(KR_p.value < .05) %>% 
  group_by(params) %>%
  count() %>%
  mutate(type1 = n/1200)

# Number of measurements 
sim_exp %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_measurements) %>%
  count() %>%
  mutate(type1 = n/1800)

sim_lnormal %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_measurements) %>%
  count() %>%
  mutate(type1 = n/2400)

sim_normal %>%
  filter(KR_p.value < .05) %>% 
  group_by(number_measurements) %>%
  count() %>%
  mutate(type1 = n/1800)

# RSlope vs intercept only 
sim_exp %>%
  filter(KR_p.value < .05) %>% 
  group_by(rslope) %>%
  count() %>%
  mutate(type1 = n/1800)

sim_lnormal %>%
  filter(KR_p.value < .05) %>% 
  group_by(rslope) %>%
  count() %>%
  mutate(type1 = n/2400)

sim_normal %>%
  filter(KR_p.value < .05) %>% 
  group_by(rslope) %>%
  count() %>%
  mutate(type1 = n/1800)




all_sim <- rbind(sim_exp,sim_lnormal, sim_normal)



# OTHER VISUALIZATIONS
# P- value histograms? # qq plot?
# do i have to do difference testing to test signficance between exp/normal/lnormal
hist(sim_exp$KR_p.value )