
library(tidyverse)
library(lmerTest)
library(SimMultiCorrData)
library(purrr)
library(broom)

## Function that returns p values from Kenward-Roger and Satterthwaite from linear mixed model 
sim_mixedeffects <- function(dist, parameters, random_slopes=TRUE,
                             n.indiv, n.measurements, corr){
   
   corr.data <- c(1,corr,corr,1)
   corr1 <- matrix(corr.data,nrow=2)
   
   # save M 
 
   m<-round(calc_theory(Dist = dist, params = parameters),8)
  
   if(random_slopes == TRUE){
     r_effects <-rcorrvar(n = n.indiv, k_cont = 2, method = "Fleishman", means = c(m[1],m[1]), vars = c(m[2]^2,m[2]^2), skews = c(m[3],m[3]),
                          skurts =  c(m[4],m[4]),
                          errorloop = TRUE, rho = corr1,seed = sample(1000:9999,1))
     
     r_effects1 <- r_effects$continuous_variables$V1 - m[1]
     r_effects2 <- r_effects$continuous_variables$V2 - m[1]
     
   } else {

    r_effects1 <-nonnormvar1("Fleishman", means = m[1], vars = m[2]^2, skews =  m[3],
                             skurts = m[4], n = n.indiv,seed = sample(1000:9999,1))$continuous_variable$V1
    
    r_effects2 <- rep(0,n.indiv)
    
   }
    
   data <- expand.grid(indiv = as.factor(1:n.indiv), time = 1:n.measurements)
   matrix_fixedeffects <- model.matrix(~ time, data)   # model matrix for fixed effects
   
   betas <- c(3.1, 0)   # fixed effects coefficient vector
   
   matrix_intercept <- model.matrix(~ 0 + indiv, data)   # model matrix for random intercepts
   
   matrix_slope <-  model.matrix(~ 0 + indiv, data) * data$time   # model matrix for random slopes
   
   errors <- rnorm(nrow(data), 0, .2)  
   
   data$fixed_effects <-  matrix_fixedeffects %*% betas
   data$intercept_wrandom <- matrix_intercept %*% r_effects1
   data$slope_wrandom <- matrix_slope %*% r_effects2
   data$errors <- errors
   data$r_effects1 <- r_effects1
   data$r_effects2 <- r_effects2
   
   data$Y_val <- matrix_fixedeffects %*% betas + matrix_intercept %*% r_effects1 + matrix_slope %*% r_effects2 + errors
   
   # construct model
   if(random_slopes == FALSE){
      model <- lmerTest::lmer(Y_val ~ time + (1|indiv) , data) 
   } else{
      model <- lmerTest::lmer(Y_val ~ time + (time|indiv) , data) 
   }
   
   # add random effects value 
   # ADD NORMAL
   a_KR <- tidy(anova(model, ddf = "Kenward-Roger"))
   a_S <- tidy(anova(model))
   p <- data.frame(distribution = dist, number_individuals = n.indiv,
                   params = paste0(parameters, collapse = ","),
                   number_measurements = n.measurements,
                   rslope = random_slopes)
   moms <- data.frame(mean = m[1], sd = m[2], skew = m[3], kurtosis = m[4],
                      fifth = m[5], sixth = m[6])
   
   colnames(a_KR) <- paste("KR", colnames(a_KR), sep = "_")
   colnames(a_S) <- paste("S", colnames(a_S), sep = "_")

  
  return(cbind(a_KR,a_S,moms,p,data)
        
  )
}

w<-sim_mixedeffects(dist = 'Exponential', parameters = .2, random_slopes = TRUE,
                 n.indiv = 10, n.measurements = 4,
                 corr = -.38)




## Plotting the exponential and lognormal distributions


#### lognormal
curve(dlnorm(x, meanlog=0, sdlog=.25), from=0, to=5)
curve(dlnorm(x, meanlog=.5, sdlog=.1), from=0, to=5)
curve(dlnorm(x, meanlog=1, sdlog=.5), from=0, to=5)
curve(dlnorm(x, meanlog=0, sdlog=.9), from=0, to=5)

## Exponential 
x <- seq(0, 8, 0.1)
plot(x, dexp(x, .2), type = "l",
     ylab = "", lwd = 2, col = "red")

plot(x, dexp(x, 4), type = "l",
     ylab = "", lwd = 2, col = "red")

plot(x, dexp(x, .9), type = "l",
     ylab = "", lwd = 2, col = "red")


## Setting up factors for lognormal condition
dist <- c("Lognormal")
corr <- c(-.38)
parameters <- list(c(0,.25), c(.5,.1),c(1,.5),c(0,.9))
random_slopes <- c(TRUE,FALSE)
n.indiv <- c(10,18,26)
n.measurements <- c(4,8)

# all possible combinations
combos_lnormal<-crossing(dist,parameters,random_slopes,n.indiv,n.measurements,corr)


## Setting up factors for exponential condition
dist <- c("Exponential")
corr <- c(-.38)
parameters <- c(4,.2,.9)
random_slopes <- c(TRUE,FALSE)
n.indiv <- c(10,18,26)
n.measurements <- c(4,8)

# all possible combinations
combos_exp <- crossing(dist,parameters,random_slopes,n.indiv,n.measurements,corr)


numsim <-10000
set.seed(2021)
## Using all possible combiniations, simulate
simulations_exp <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_exp$dist,combos_exp$parameters,combos_exp$random_slopes,combos_exp$n.indiv,combos_exp$n.measurements,combos_exp$corr),.f = sim_mixedeffects, .id = "sim_num"))

simulations_lnormal <- 1:numsim %>%
   map_df(~ pmap_dfr(list(combos_lnormal$dist,combos_lnormal$parameters,combos_lnormal$random_slopes,combos_lnormal$n.indiv,combos_lnormal$n.measurements,combos_lnormal$corr),.f = sim_mixedeffects))


#saveRDS(simulations_exp,"~/jeyuthesis/SimulationCode/exponential.rds")
#saveRDS(simulations_lnormal,"~/jeyuthesis/SimulationCode/lnormal.rds")

######

## Type 1 
exponential %>%
   filter(S_p.value < .05
   ) %>%
   nrow()

4/240


exponential %>%
   filter(KR_p.value < .05,
          params ==4) %>%
   nrow()

exponential %>%
   filter(KR_p.value < .05,
         params ==.9) %>%
   nrow()

#########


lnormal %>%
   filter(S_p.value < .05
   ) %>%
   nrow()
3/240


#### testing params


m<-calc_theory(Dist = c("Lognormal"), params = c(0,.9))
nonnormvar1("Fleishman", means = m[1], vars = m[2]^2, skews =  m[3],
            skurts = m[4], n = n.indiv,seed = sample(1000:9999,1))$continuous_variable$V1
