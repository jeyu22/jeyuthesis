
library(tidyverse)
library(lmerTest)
library(SimMultiCorrData)
library(purrr)
library(broom)
library(broom.mixed)

## Function that returns output from LME model using KR & Satterthwaite DF
sim_mixedeffects <- function(dist, parameters, random_slopes = TRUE,
                             n.indiv, n.measurements, corr) {
  corr.data <- c(1, corr, corr, 1)
  corr1 <- matrix(corr.data, nrow = 2)

  # save kurtosis and skew values
  m <- round(calc_theory(Dist = dist, params = parameters), 8)

  if (random_slopes == TRUE) {
    r_effects <- rcorrvar(
      n = n.indiv, k_cont = 2, method = "Fleishman", 
      means = c(m[1], m[1]), vars = c(m[2]^2, (m[2]^2) / 18), 
      skews = c(m[3], m[3]),
      skurts = c(m[4], m[4]),
      errorloop = TRUE, rho = corr1, seed = sample(1000:9999, 1)
    )
  
    # centering the mean at 0 
    r_effects1 <- r_effects$continuous_variables$V1 - m[1]
    r_effects2 <- r_effects$continuous_variables$V2 - m[1]

    errors <- rnorm(n.indiv * n.measurements, 0, m[2] * sqrt(.168))
  } else { 
    # random intercept model 
    r_effects1 <- nonnormvar1("Fleishman",
      means = m[1], vars = m[2]^2, skews = m[3],
      skurts = m[4], n = n.indiv, seed = sample(1000:9999, 1)
    )$continuous_variable$V1

    r_effects1 <- r_effects1 - m[1]
    r_effects2 <- rep(0, n.indiv)
    
    errors <- rnorm(n.indiv * n.measurements, 0, m[2] * sqrt(.4))
  }

  data <- expand.grid(time = 1:n.measurements, indiv = as.factor(1:n.indiv))
  data$treatment <- rep(0:1, each = (n.indiv * n.measurements) / 2)
  
  # model matrix for fixed effects
  matrix_fixedeffects <- model.matrix(~ time + treatment, data) 

  # fixed effects coefficient vector
  betas <- c(3.1, 0, 0) 

  # model matrix for random intercepts
  matrix_intercept <- model.matrix(~ 0 + indiv, data)

  # model matrix for random slopes
  matrix_slope <- model.matrix(~ 0 + indiv, data) * data$time


  data$Y_val <- (matrix_fixedeffects %*% betas +
    matrix_intercept %*% r_effects1 +
    matrix_slope %*% r_effects2 + errors)[, 1]



  # construct model
  if (random_slopes == FALSE) {
    model <- lmerTest::lmer(Y_val ~ time + (1 | indiv) + treatment, data)
  } else {
    model <- lmerTest::lmer(Y_val ~ time + (time | indiv) + treatment, data)
  }

  # extract p-values/summary output from KR and Satterthwaite
  KR <- broom.mixed::tidy(model, ddf = "Kenward-Roger")
  colnames(KR) <- paste("KR", colnames(KR), sep = "_")

  S <- broom.mixed::tidy(model, ddf = "Satterthwaite") %>%
    rename_with(.cols = everything(), .fn = ~ str_c("S_", .x)) %>%
    add_column(
      .before = "S_effect",
      distribution = dist, number_individuals = n.indiv,
      params = paste0(parameters, collapse = ","),
      number_measurements = n.measurements,
      rslope = random_slopes, mean = m[1], sd = m[2], skew = m[3], 
      kurtosis = m[4],
      fifth = m[5], sixth = m[6]
    )


  return(cbind(S, KR))
}



## Setting up factors for lognormal condition
dist <- c("Lognormal")
corr <- c(-.38)
parameters <- list(c(0, .25), c(1, .5))
random_slopes <- c(TRUE, FALSE)
n.indiv <- c(10, 18, 26)
n.measurements <- c(4, 8)

# all possible combinations
combos_lnormal <- crossing(
  dist, parameters, random_slopes,
  n.indiv, n.measurements, corr
)


## Setting up factors for exponential condition
dist <- c("Exponential")
corr <- c(-.38)
parameters <- c(.2, .9)
random_slopes <- c(TRUE, FALSE)
n.indiv <- c(10, 18, 26)
n.measurements <- c(4, 8)

# all possible combinations
combos_exp <- crossing(
  dist, parameters, random_slopes,
  n.indiv, n.measurements, corr
)


numsim <- 80 

###### Lognormal Distribution: 400 total simulations
set.seed(2021)

simulations_lnormal_1 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_lnormal$dist, combos_lnormal$parameters,
    combos_lnormal$random_slopes, combos_lnormal$n.indiv,
    combos_lnormal$n.measurements, combos_lnormal$corr
  ),
  .f = sim_mixedeffects
  ))

# saveRDS(simulations_lnormal_1, "simulations_lnormal_1.rds")
set.seed(2005)

simulations_lnormal_2 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_lnormal$dist, combos_lnormal$parameters,
    combos_lnormal$random_slopes, combos_lnormal$n.indiv,
    combos_lnormal$n.measurements, combos_lnormal$corr
  ),
  .f = sim_mixedeffects
  ))
# saveRDS(simulations_lnormal_2, "simulations_lnormal_2.rds")

set.seed(2011)

simulations_lnormal_3 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_lnormal$dist, combos_lnormal$parameters,
    combos_lnormal$random_slopes, combos_lnormal$n.indiv,
    combos_lnormal$n.measurements, combos_lnormal$corr
  ),
  .f = sim_mixedeffects
  ))
# saveRDS(simulations_lnormal_3, "simulations_lnormal_3.rds")

set.seed(1012)
simulations_lnormal_4 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_lnormal$dist, combos_lnormal$parameters,
    combos_lnormal$random_slopes, combos_lnormal$n.indiv,
    combos_lnormal$n.measurements, combos_lnormal$corr
  ),
  .f = sim_mixedeffects
  ))
# saveRDS(simulations_lnormal_4, "simulations_lnormal_4.rds")

set.seed(3434)
simulations_lnormal_5 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_lnormal$dist, combos_lnormal$parameters,
    combos_lnormal$random_slopes, combos_lnormal$n.indiv,
    combos_lnormal$n.measurements, combos_lnormal$corr
  ),
  .f = sim_mixedeffects
  ))
# saveRDS(simulations_lnormal_5, "simulations_lnormal_5.rds")


########### Normal distribution: 400 simulations

set.seed(7006)
simulations_normal_1 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(
    combos_normal$dist, combos_normal$parameters,
    combos_normal$random_slopes, combos_normal$n.indiv,
    combos_normal$n.measurements, combos_normal$corr
  ),
  .f = sim_mixedeffects
  ))
# saveRDS(simulations_normal_1, "simulations_normal_1.rds")



set.seed(3692)
simulations_normal_2 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_normal$dist, combos_normal$parameters, 
                         combos_normal$random_slopes, combos_normal$n.indiv, 
                         combos_normal$n.measurements, combos_normal$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_normal_2, "simulations_normal_2.rds")

set.seed(0814)
simulations_normal_3 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_normal$dist, combos_normal$parameters, 
                         combos_normal$random_slopes, combos_normal$n.indiv, 
                         combos_normal$n.measurements, combos_normal$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_normal_3, "simulations_normal_3.rds")

set.seed(0904)
simulations_normal_4 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_normal$dist, combos_normal$parameters, 
                         combos_normal$random_slopes, combos_normal$n.indiv, 
                         combos_normal$n.measurements, combos_normal$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_normal_4, "simulations_normal_4.rds")

set.seed(1111)
simulations_normal_5 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_normal$dist, combos_normal$parameters, 
                         combos_normal$random_slopes, combos_normal$n.indiv, 
                         combos_normal$n.measurements, combos_normal$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_normal_5, "simulations_normal_5.rds")


######### Exponential Distribution: 400 simulations

set.seed(9998)
simulations_exp_1 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_exp$dist, combos_exp$parameters, 
                         combos_exp$random_slopes, combos_exp$n.indiv, 
                         combos_exp$n.measurements, combos_exp$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_exp_1, "simulations_exp_1.rds")


set.seed(9997)
simulations_exp_2 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_exp$dist, combos_exp$parameters, 
                         combos_exp$random_slopes, combos_exp$n.indiv, 
                         combos_exp$n.measurements, combos_exp$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_exp_2, "simulations_exp_2.rds")

set.seed(4035)
simulations_exp_3 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_exp$dist, combos_exp$parameters, 
                         combos_exp$random_slopes, combos_exp$n.indiv, 
                         combos_exp$n.measurements, combos_exp$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_exp_3, "simulations_exp_3.rds")

set.seed(42847)
simulations_exp_4 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_exp$dist, combos_exp$parameters, 
                         combos_exp$random_slopes, combos_exp$n.indiv, 
                         combos_exp$n.measurements, combos_exp$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_exp_4, "simulations_exp_4.rds")

set.seed(0505)
simulations_exp_5 <- 1:numsim %>%
  map_df(~ pmap_dfr(list(combos_exp$dist, combos_exp$parameters, 
                         combos_exp$random_slopes, combos_exp$n.indiv, 
                         combos_exp$n.measurements, combos_exp$corr), 
                    .f = sim_mixedeffects))
# saveRDS(simulations_exp_5, "simulations_exp_5.rds")


