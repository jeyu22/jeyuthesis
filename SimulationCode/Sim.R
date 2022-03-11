
library(tidyverse)
library(lmerTest)
library(SimMultiCorrData)
library(purrr)
library(broom)


## Function that returns p values from Kenward-Roger and Satterthwaite from linear mixed model
sim_mixedeffects <- function(dist, parameters, random_slopes = TRUE,
                             n.indiv, n.measurements, corr) {
   
  corr.data <- c(1, corr, corr, 1)
  corr1 <- matrix(corr.data, nrow = 2)

   #save kurtosis and skew values
  m <- round(calc_theory(Dist = dist, params = parameters), 8)

  if (random_slopes == TRUE) {
    r_effects <- rcorrvar(
      n = n.indiv, k_cont = 2, method = "Fleishman", means = c(m[1], m[1]), vars = c(m[2]^2, m[2]^2), skews = c(m[3], m[3]),
      skurts = c(m[4], m[4]),
      errorloop = TRUE, rho = corr1, seed = sample(1000:9999, 1)
    )

    r_effects1 <- r_effects$continuous_variables$V1 - m[1]
    r_effects2 <- r_effects$continuous_variables$V2 - m[1]
  } else {
    r_effects1 <- nonnormvar1("Fleishman",
      means = m[1], vars = m[2]^2, skews = m[3],
      skurts = m[4], n = n.indiv, seed = sample(1000:9999, 1)
    )$continuous_variable$V1

    r_effects2 <- rep(0, n.indiv)
  }

  data <- expand.grid(indiv = as.factor(1:n.indiv), time = 1:n.measurements)
  matrix_fixedeffects <- model.matrix(~time, data) # model matrix for fixed effects

  betas <- c(3.1, 0) # fixed effects coefficient vector

  # model matrix for random intercepts
  matrix_intercept <- model.matrix(~ 0 + indiv, data) 

  # model matrix for random slopes
  matrix_slope <- model.matrix(~ 0 + indiv, data) * data$time 

  
  errors <- rnorm(nrow(data), 0, .2)


  data$Y_val <- matrix_fixedeffects %*% betas + 
     matrix_intercept %*% r_effects1 + 
     matrix_slope %*% r_effects2 + errors

  # construct model
  if (random_slopes == FALSE) {
    model <- lmerTest::lmer(Y_val ~ time + (1 | indiv), data)
    var_ranef_df <- data.frame(VarCorr(model))
    var_ranef <- data.frame(
      var_rand_intercept = var_ranef_df[1, 4],
      var_rand_residual = var_ranef_df[2, 4]
    )
  } else {
    model <- lmerTest::lmer(Y_val ~ time + (time | indiv), data)
    var_ranef_df <- data.frame(VarCorr(model))
    var_ranef <- data.frame(
      var_rand_intercept = var_ranef_df[1, 4],
      var_rand_slope = var_ranef_df[2, 4],
      var_rand_residual = var_ranef_df[4, 4]
    )
  }


  a_KR <- tidy(anova(model, ddf = "Kenward-Roger"))
  a_S <- tidy(anova(model))
  p <- data.frame(
    distribution = dist, number_individuals = n.indiv,
    params = paste0(parameters, collapse = ","),
    number_measurements = n.measurements,
    rslope = random_slopes, mean = m[1], sd = m[2], skew = m[3], kurtosis = m[4],
    fifth = m[5], sixth = m[6]
  )


  colnames(a_KR) <- paste("KR", colnames(a_KR), sep = "_")
  colnames(a_S) <- paste("S", colnames(a_S), sep = "_")


  return(cbind(a_KR, a_S, p, var_ranef))
}

sim_mixedeffects(
  dist = "Gaussian", parameters = c(0, 2), random_slopes = TRUE,
  n.indiv = 10, n.measurements = 4,
  corr = -.38
)




## Plotting the exponential and lognormal distributions


#### lognormal
curve(dlnorm(x, meanlog = 0, sdlog = .25), from = 0, to = 5)
curve(dlnorm(x, meanlog = .5, sdlog = .1), from = 0, to = 5)
curve(dlnorm(x, meanlog = 1, sdlog = .5), from = 0, to = 5)
curve(dlnorm(x, meanlog = 0, sdlog = .9), from = 0, to = 5)

## Exponential
x <- seq(0, 8, 0.1)
plot(x, dexp(x, .2),
  type = "l",
  ylab = "", lwd = 2, col = "red"
)

plot(x, dexp(x, 4),
  type = "l",
  ylab = "", lwd = 2, col = "red"
)

plot(x, dexp(x, .9),
  type = "l",
  ylab = "", lwd = 2, col = "red"
)


## Setting up factors for lognormal condition
dist <- c("Lognormal")
corr <- c(-.38)
parameters <- list(c(0, .25), c(.5, .1), c(1, .5), c(0, .9))
random_slopes <- c(TRUE, FALSE)
n.indiv <- c(10, 18, 26)
n.measurements <- c(4, 8)

# all possible combinations
combos_lnormal <- crossing(dist, parameters, random_slopes, 
                           n.indiv, n.measurements, corr)


## Setting up factors for exponential condition
dist <- c("Exponential")
corr <- c(-.38)
parameters <- c(4, .2, .9)
random_slopes <- c(TRUE, FALSE)
n.indiv <- c(10, 18, 26)
n.measurements <- c(4, 8)

# all possible combinations
combos_exp <- crossing(dist, parameters, random_slopes, 
                       n.indiv, n.measurements, corr)


numsim <- 10

######
set.seed(2021)

simulations_lnormal_1 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))

#saveRDS(simulations_lnormal_1, "simulations_lnormal_1.rds")
set.seed(2005)

simulations_lnormal_2 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_2, "simulations_lnormal_2.rds")

set.seed(2011)

simulations_lnormal_3 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_3, "simulations_lnormal_3.rds")

set.seed(1012)
simulations_lnormal_4 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_4, "simulations_lnormal_4.rds")

set.seed(3434)
simulations_lnormal_5 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_5, "simulations_lnormal_5.rds")

set.seed(8190)
simulations_lnormal_6 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_6, "simulations_lnormal_6.rds")


set.seed(7788)
simulations_lnormal_7 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_7, "simulations_lnormal_7.rds")

set.seed(0317)
simulations_lnormal_8 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_8, "simulations_lnormal_8.rds")

set.seed(0618)
simulations_lnormal_9 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_9, "simulations_lnormal_9.rds")

set.seed(2222)
simulations_lnormal_10 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,
                          combos_lnormal$random_slopes,combos_lnormal$n.indiv,
                          combos_lnormal$n.measurements,combos_lnormal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_lnormal_10, "simulations_lnormal_10.rds")

## 

set.seed(7006)
simulations_normal_1 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,
                          combos_normal$random_slopes,combos_normal$n.indiv,
                          combos_normal$n.measurements,combos_normal$corr),
                     .f = sim_mixedeffects))
#saveRDS(simulations_normal_1, "simulations_normal_1.rds")



set.seed(3692)
simulations_normal_2 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_2, "simulations_normal_2.rds")

set.seed(0814)
simulations_normal_3 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_3, "simulations_normal_3.rds")

set.seed(0904)
simulations_normal_4 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_4, "simulations_normal_4.rds")

set.seed(1111)
simulations_normal_5 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_5, "simulations_normal_5.rds")

set.seed(3333)
simulations_normal_6 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_6, "simulations_normal_6.rds")

set.seed(4444)
simulations_normal_7 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_7, "simulations_normal_7.rds")

set.seed(6666)
simulations_normal_8 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_8, "simulations_normal_8.rds")

set.seed(7777)
simulations_normal_9 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_9, "simulations_normal_9.rds")

set.seed(8888)
simulations_normal_10 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_normal$dist,combos_normal$parameters,combos_normal$random_slopes,combos_normal$n.indiv,combos_normal$n.measurements,combos_normal$corr),.f = sim_mixedeffects))
#saveRDS(simulations_normal_10, "simulations_normal_10.rds")

#########

set.seed(9998)
simulations_exp_1 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_1, "simulations_exp_1.rds")


set.seed(9997)
simulations_exp_2 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_2, "simulations_exp_2.rds")

set.seed(4035)
simulations_exp_3 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_3, "simulations_exp_3.rds")

set.seed(42847)
simulations_exp_4 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_4, "simulations_exp_4.rds")

set.seed(0505)
simulations_exp_5 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_5, "simulations_exp_5.rds")


set.seed(9994)
simulations_exp_6 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_6, "simulations_exp_6.rds")

set.seed(9993)
simulations_exp_7 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_7, "simulations_exp_7.rds")

set.seed(3009)
simulations_exp_8 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_8, "simulations_exp_8.rds")

set.seed(2203)
simulations_exp_9 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_9, "simulations_exp_9.rds")

set.seed(9991)
simulations_exp_10 <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))
#saveRDS(simulations_exp_10, "simulations_exp_10.rds")


