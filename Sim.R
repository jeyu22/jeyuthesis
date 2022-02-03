
library(tidyverse)
library(lmerTest)
library(SimMultiCorrData)
library(purrr)
m <-calc_theory(Dist = "Exponential", params = .5)

r_effects <-rcorrvar(n = 12, k_cont = 2, method = 'Fleishman', means = c(m[1],m[1]), vars = c(m[2],m[2]), skews = c(m[3],m[3]),
                     skurts =  c(m[4],m[4]),
                     errorloop = TRUE, rho = corr1)


sim_mixedeffects <- function(dist, parameters, random_slopes=TRUE,
                             n.indiv, n.measurements, corr){
   
   corr.data <- c(1,corr,corr,1)
   corr1 <- matrix(corr.data,nrow=2)
 
   m<-round(calc_theory(Dist = dist, params = parameters),8)
  
   if(random_slopes == TRUE){
     r_effects <-rcorrvar(n = n.indiv, k_cont = 2, method = "Fleishman", means = c(m[1],m[1]), vars = c(m[2],m[2]), skews = c(m[3],m[3]),
                          skurts =  c(m[4],m[4]),
                          errorloop = TRUE, rho = corr1,seed = sample(1000:9999,1))
     
     r_effects1 <- r_effects$continuous_variables$V1
     r_effects2 <- r_effects$continuous_variables$V2
     
   } else {

    r_effects1 <-nonnormvar1("Fleishman", means = m[1], vars = m[2]^2, skews =  m[3],
                             skurts = m[4], n = n.indiv,seed = sample(1000:9999,1))$continuous_variable$V1
    
    r_effects2 <- rep(0,n.indiv)
    
   }
    
   data <- expand.grid(indiv = LETTERS[1:n.indiv], time = 1:n.measurements)
   matrix_fixedeffects <- model.matrix(~ time, data)   # model matrix for fixed effects
   
   betas <- c(3.1, 0)   # fixed effects coefficient vector
   
   matrix_intercept <- model.matrix(~ 0 + indiv, data)   # model matrix for random intercepts
   
   matrix_slope <-  model.matrix(~ 0 + indiv, data) * data$time   # model matrix for random slopes
   
   errors <- rnorm(nrow(data), 0, 2)  
   
   data$Y_val <- matrix_fixedeffects %*% betas + matrix_intercept %*% r_effects1 + matrix_slope %*% r_effects2 + errors
   
   # construct model
   model <- lmerTest::lmer(Y_val ~ time + (time|indiv) , data) 
   

  
  return(data.frame(distribution = dist, number_individuals = n.indiv,
                    params = parameters,
                       number_measurements = n.measurements,
                       rslope = random_slopes, p_value_KR = anova(model,ddf = "Kenward-Roger")$`Pr(>F)`,
         p_value_S =anova(model,ddf = "Satterthwaite")$`Pr(>F)`,
         f_value_KR = anova(model, ddf = "Kenward-Roger")$`F value`,
         f_value_S = anova(model)$`F value`,
         df_S = anova(model)$`DenDF`,
         df_KR = anova(model, ddf = "Kenward-Roger")$`DenDF`,
         msq_KR = anova(model, ddf = "Kenward-Roger")$`Mean Sq`,
         msq_S =anova(model)$`Mean Sq` ))
}


curve(dlnorm(x, meanlog=0, sdlog=.25), from=0, to=5)
curve(dlnorm(x, meanlog=.5, sdlog=.1), from=0, to=5)
curve(dlnorm(x, meanlog=.1, sdlog=.5), from=0, to=5)
x <- seq(0, 8, 0.1)
plot(x, dexp(x, .2), type = "l",
     ylab = "", lwd = 2, col = "red")

dist <- c("Lognormal")
corr <- c(-.38)
parameters <- list(c(0,.25), c(.5,.1))
random_slopes <- c(TRUE,FALSE)
n.indiv <- c(10,18,26)
n.measurements <- c(4,8)

combos_lnormal<-crossing(dist,parameters,random_slopes,n.indiv,n.measurements,corr)

dist <- c("Exponential")
corr <- c(-.38)
parameters <- c(4,.2)
random_slopes <- c(TRUE,FALSE)
n.indiv <- c(10,18,26)
n.measurements <- c(4,8)

combos_exp <- crossing(dist,parameters,random_slopes,n.indiv,n.measurements,corr)



numsim <-10
simulations_exp <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects))

simulations_lnormal <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,combos_lnormal$random_slopes,combos_lnormal$n.indiv,combos_lnormal$n.measurements,combos_lnormal$corr),.f = sim_mixedeffects))

## WOrking params
# log normal (10,6), (0,.25)


saveRDS(simulations_exp,"exponential.rds")
saveRDS(simulations_lnormal,"lnormal.rds")

sim_mixedeffects(dist = "Exponential", parameters = .2,random_slopes = TRUE,
                 n.indiv = 24, n.measurements = 4, corr = -.38
                 )
