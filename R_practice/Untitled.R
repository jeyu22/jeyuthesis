library(tidyverse)

load("/Users/jessicayu/Downloads/ICPSR_21600/DS0008/21600-0008-Data.rda")
ds3 <-da21600.0008

asian3 <- ds3 %>%
  filter(H3OD4D == "(1) (1) Marked", BIO_SEX3 == "(2) (2) Female",
         H3OD5A == "(1) (1) Marked")

lbls <- sort(levels(asian3$H3HGT_F))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian3$H3HGT_F <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian3$H3HGT_F))
asian3$H3HGT_F <- add.value.labels(asian3$H3HGT_F, lbls)

lbls <- sort(levels(asian3$H3HGT_I))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian3$H3HGT_I <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian3$H3HGT_I))
asian3$H3HGT_I <- add.value.labels(asian3$H3HGT_I, lbls)

asian3 %>%
  select(H3WGT)

# weight, ID 
asian3_new <- asian3 %>%
  select(H3WGT,AID,H3HGT_F,H3HGT_I,H3SP21,CALCAGE3,H3WP10Y) %>%
  mutate(kg = H3WGT/2.2046, m = (H3HGT_I + (H3HGT_F*12))*0.0254) %>%
  mutate(bmi = (kg/(m^2)))



gf_density(data = asian3 , ~bmi)

load("/Users/jessicayu/Downloads/ICPSR_21600 WAVE4/DS0022/21600-0022-Data.rda")
ds4 <- da21600.0022

asian4 <- ds4 %>%
  filter(AID %in% asian3$AID) %>%
  select(AID, H4BMI)

gf_density(data = asian4 , ~H4BMI)

################ WAVE 2
load("/Users/jessicayu/Downloads/ICPSR_21600 3/DS0005/21600-0005-Data.rda")
ds2 <-da21600.0005

asian2 <- ds2 %>%
  filter(AID %in% asian1$AID
         ) 



lbls <- sort(levels(asian2$H2GH52F))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian2$H2GH52F <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian2$H2GH52F))
asian2$H2GH52F <- add.value.labels(asian2$H2GH52F, lbls)

lbls <- sort(levels(asian2$H2GH52I))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian2$H2GH52I <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian2$H2GH52I))
asian2$H2GH52I <- add.value.labels(asian2$H2GH52I, lbls)



# weight, ID 
asian2_new <- asian2 %>%
  select(H2GH53,AID,H2GH52F,H2GH52I,H2WP7,H2PF24,H2DA8) %>%
  mutate(kg = H2GH53/2.2046, m = (H2GH52I + (H2GH52F*12))*0.0254) %>%
  mutate(bmi = (kg/(m^2)))

gf_density(data = asian2_new , ~bmi)




##################### WAVE 1 

load("/Users/jessicayu/Downloads/ICPSR_21600 5/DS0001/21600-0001-Data.rda")


ds1 <-da21600.0001

asian1 <- ds1 %>%
  filter(H1GI7A == "(1) (1) Marked",BIO_SEX == "(2) (2) Female",
         AID != "90577229", AID != "94577129", AID != "94576044",
         AID != "96578138", AID != "99575118", AID != "99576347",
         AID %in% 
           asian3_new_glue$AID)

# 23 > 16 > 15


lbls <- sort(levels(asian1$H1GH59A))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1GH59A <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1GH59A))
asian1$H1GH59A <- add.value.labels(asian1$H1GH59A, lbls)

lbls <- sort(levels(asian1$H1GH59B))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1GH59B <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1GH59B))
asian1$H1GH59B <- add.value.labels(asian1$H1GH59B, lbls)

lbls <- sort(levels(asian1$H1WP7))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1WP7 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1WP7))
asian1$H1WP7 <- add.value.labels(asian1$H1WP7, lbls)

lbls <- sort(levels(asian1$H1WP7))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1WP7 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1WP7))
asian1$H1WP7 <- add.value.labels(asian1$H1WP7, lbls)


lbls <- sort(levels(asian1$H1DS7))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1DS7 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1DS7))
asian1$H1DS7 <- add.value.labels(asian1$H1DS7, lbls)


lbls <- sort(levels(asian1$H1TO1))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1TO1 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1TO1))
asian1$H1TO1 <- add.value.labels(asian1$H1TO1, lbls)

lbls <- sort(levels(asian1$H1FV5))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1FV5 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1FV5))
asian1$H1FV5 <- add.value.labels(asian1$H1FV5, lbls)

lbls <- sort(levels(asian1$H1TO12))
lbls <- (sub("^\\([0-9]+\\) +(.+$)", "\\1", lbls))
asian1$H1TO12 <- as.numeric(sub("^\\(0*([0-9]+)\\).+$", "\\1", asian1$H1TO12))
asian1$H1TO12 <- add.value.labels(asian1$H1TO12, lbls)

# weight, ID 
asian1_new <- asian1 %>%
  select(H1GH60,AID,H1GH59A,H1GH59B, H1WP7,H1DS7,
         H1TO1,H1FV5,H1TO12) %>%
  mutate(kg = H1GH60/2.2046, m = (H1GH59B + (H1GH59A*12))*0.0254) %>%
  mutate(bmi = (kg/(m^2))) %>%
  select(-m) %>% select(-kg) %>% select(-H1GH60) %>% select(-H1GH59A) %>%
  select(-H1GH59B)

asian1_new <- asian1_new %>%
  rename(decisions_eat = H1WP7, 
         run_away = H1DS7,
         alcohol = H1TO12,
         physical_fight = H1FV5,
         cig = H1TO1
         ) %>%
  mutate(wave = 1,
         decisions_eat = factor(decisions_eat),
         run_away = factor(run_away),
         alcohol = factor(alcohol),
         physical_fight = factor(physical_fight),
         cig = factor(cig),
         AID = factor(AID))


ggformula::gf_density(data = asian1_new , ~bmi)



# HOW CLOSE TO BIO MOM  W1,W2,W3  

# DECISIONS ABOUT WHAT TO EAT W1,W2 H1WP7, 

# How many hours a week do you watch television? W1 W2 W3 W4 H1DA8,H2DA8
#H3DA7, H4DA1


# H1DS7 run away from home?


# H1TO4 age cigs

# H1TO1 drink 

#H1FV1 you saw someone stab

# H1FV3 someone SHOT YOU
#H1FV4 SOMEONE STABBED YOU


###### COMBING ALL OF THEM

asian1_copy <- asian1_new %>% select(-bmi)

asian2_new_glue <- asian1_copy %>%
  left_join(asian2_new %>% select(AID, bmi), by = c("AID" = "AID")) %>%
  mutate(wave = 2)


asian3_new_glue <- asian1_copy %>%
  left_join(asian3_new %>% select(AID, bmi), by = c("AID" = "AID")) %>%
  mutate(wave = 3) 

asian2 %>% select(CALCAGE2)
age<-asian3 %>% dplyr::select(H3OD1Y,AID)

 
asian_women_ds <- rbind(asian1_new, asian2_new_glue, asian3_new_glue)
asian_women_ds <- asian_women_ds %>%
  left_join(age, by = "AID") %>%
  mutate(study_year = case_when(wave == 1 ~ 1995,
                                wave == 2 ~ 1996,
                               wave == 3 ~ 2002 ),
         real_age = (study_year - H3OD1Y)) 

asian_women_ds <- asian_women_ds %>%
  dplyr::select(-study_year) %>%dplyr::select(-H3OD1Y)
  

saveRDS(asian_women_ds,"child_data.Rds")
model <- lmer(bmi ~  decisions_eat *  run_away * cig * alcohol * physical_fight * (1|AID) * (1|wave), asian_women_ds) 
summary(model)
v<-as.data.frame(VarCorr(model,comp=c("Variance","Std.Dev.")))
## 


anova(model, ddf = "Kenward-Roger")
anova(model, ddf = "Satterthwaite")

mat2.data <- c(w$constants,w$constants)
mat2 <- matrix(mat2.data,nrow=2)
mat2
cor1<-cov2cor(mat2)
cor1

mat1.data <- c(1,.2,.2,1)
mat1 <- matrix(mat1.data,nrow=2)
mat1
x = rlnorm(500,1,.6)
grid = seq(0,25,.1)

plot(grid,dlnorm(grid,1,.6),type="l",xlab="x",ylab="f(x)")
lines(density(x),col="red")

w<-find_constants(method = 'Fleishman', skews= 0.9527364, skurts = 0.3887563, seed = 15)
findintercorr(n = 12, k_cont = 2,method = 'Fleishman', rho = mat1 ,
              constants =mat2)
rcorrvar(n = 100, k_cont = 2, method = 'Polynomial', means = c(0,0), vars = c(1,1), skew = 0.952,
         skurts = 1.75, fifths = c(-3.4346798,-3.4346798), sixths = c(-13.8229145 ,-13.8229145 ),
          errorloop = TRUE, rho = mat1
         )
valid_corr(n = 12, k_cont = 2, method = "Fleishman", means = c(0,0), vars = c(1,1), skews = 2.1,
           skurts =  6.1, rho = mat1)


round(calc_theory(Dist = "Lognormal", params = c(3.02092234,0.16299311 )),8)
rlno
ncol(mat2)

library(MASS)
fitdistr(asian_women_ds$bmi,densfun = "log-normal")

calc_moments(asian_women_ds$bmi)

H_exp <- nonnormvar1("Fleishman", means = 0, vars = 1, skews =  0.9527364,
                     skurts = 0.3887563, , n = 12, seed = 1234)
H_exp2 <- nonnormvar1("Fleishman", means = 0, vars = 1, skews =  0.9527364,
                     skurts = 0.3887563, , n = 12, seed = 3000)
data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)

calc_lower_skurt(method = "Fleishman", skews = 0.9527364)
