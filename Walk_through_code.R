##########################################################################################################
##########################################################################################################
#########################                                                        #########################
#########################            Age, Growth, and Reproduction of            #########################
#########################                   Centrophorus uyato                   #########################
#########################                                                        #########################
#########################                                                        #########################
#########################             Brian J Moe, Charles F Cotton,             #########################
#########################                   and Joseph Travis                    #########################
#########################                                                        #########################
##########################################################################################################
##########################################################################################################
#########################                                                        #########################
#########################                      Results Code                      #########################
#########################                                                        #########################
##########################################################################################################
##########################################################################################################
rm(list = ls())

## Required packages

# Statistical ~~~~
## "stats"
## "tidyverse"
## "mvtnorm"
## "tmvtnorm"
## "mgcv"
## "dglm"
## "msm"

# Graphical ~~~~
## "ggplot2"
## "patchwork"
## "wesanderson"


### load data provided ###

# "length_age.Rdata"
# "embryo_lengths.Rdata"
# "fecundity.Rdata"

### source functions ###

# Functions 1. Length structure.R
# Functions 2. Length-Weight regressions.R
# Functions 3. Maturity.R
# Functions 4. Frequentist Growth.R
# Functions 5. Bayesian Growth.R

##### Length Structure #####
## ~~~~~~~~~~~~~~~~~~~~~~ ##

length_dist <- length_structure(FL, Sex, length_age.dat)

#### Length-Length Regressions ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

FL_TL <- Length_Length(FL, TL, Sex, length_age.dat)
FL_PCL <- Length_Length(FL, PCL, Sex, length_age.dat)
TL_PCL <- Length_Length(TL, PCL, Sex, length_age.dat)


#### Length-Weight Regressions ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
length_weight <- Length_Weight(FL, BW, Sex, length_age.dat)


#### Fecundity ####
### ~~~~~~~~~~~ ###
ovarian <- fecundity.dat[, c("Left.ovary", "Right.ovary", "Ovarian.f")] %>% replace(is.na(.), 0)
uterine <- fecundity.dat[, c("Left.uterus", "Right.uterus", "Uterine.f")] %>% replace(is.na(.), 0)

repro <- fecundity_analysis(ovarian, uterine)


#### Birth Size ####
### ~~~~~~~~~~~~ ###

birth_logit <- MC_logit_birth(embryo.dat$FL, length_age.dat$FL, step = 0.095)


#### Length/Age at Maturity ####
### ~~~~~~~~~~~~~~~~~~~~~~~~ ###

L50_logit <- MC_logit_x50(FL, Mat, Sex, length_age.dat, step = c(0.093, 0.059))

A50_logit_Age1 <- MC_logit_x50(Age1, Mat, Sex, length_age.dat, step = c(0.092, 0.059))

A50_logit_Age2 <- MC_logit_x50(Age2, Mat, Sex, length_age.dat, step = c(0.080, 0.045))


#### Identify Bayesian growth priors ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

fem.L50 <- L50_logit$Females$Bootstrap.Chain$Raw.Parameters$x50
male.L50 <- L50_logit$Males$Bootstrap.Chain$Raw.Parameters$x50

fem.t50_Age1 <- A50_logit_Age1$Females$Bootstrap.Chain$Raw.Parameters$x50
male.t50_Age1 <- A50_logit_Age1$Males$Bootstrap.Chain$Raw.Parameters$x50

fem.t50_Age2 <- A50_logit_Age2$Females$Bootstrap.Chain$Raw.Parameters$x50
male.t50_Age2 <- A50_logit_Age2$Males$Bootstrap.Chain$Raw.Parameters$x50

L0.prior <- birth_logit$Bootstrap$Raw.Parameters$x50


#### Fit 4-parameter Bayesian Growth models ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

MH_chains = 9999
n_independence = 999
n_bootstraps = 999
LOO_chains = 999

## Females -----------

### Age1 Band count ----
Female_Age1_4par <- 
  fit_Bayes_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 1), 
    L0.prior, 
    fem.L50, 
    fem.t50_Age1,
    step = c(0.031, 0.029, 0.027),
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
    )

### Age2 band count ----
Female_Age2_4par <- 
  fit_Bayes_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 1), 
    L0.prior, 
    fem.L50, 
    fem.t50_Age2,
    step = c(0.031, 0.030, 0.029),
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )

## Males -------------

### Age1 Band count ----
Male_Age1_4par <- 
  fit_Bayes_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 0), 
    L0.prior, 
    male.L50, 
    male.t50_Age1,
    step = c(0.024, 0.021, 0.016),    
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )

### Age2 band count ----
Male_Age2_4par <- 
  fit_Bayes_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 0), 
    L0.prior, 
    male.L50, 
    male.t50_Age2,
    step = c(0.0275, 0.022, 0.019),   
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )


#### Fit 3-parameter Bayesian Growth models ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

## Females -----------

### Age1 Band count ----
Female_Age1_3par <- 
  fit_Bayes_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 1), 
    L0.prior, 
    fem.L50, 
    fem.t50_Age1,
    step = c(0.074, 0.073, 0.066),
    k_based = TRUE,
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )

### Age2 band count ----
Female_Age2_3par <- 
  fit_Bayes_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 1), 
    L0.prior, 
    fem.L50, 
    fem.t50_Age2,
    step = c(0.069, 0.077, 0.073),
    k_based = TRUE,   
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )

## Males -------------

### Age1 Band count ----
Male_Age1_3par <- 
  fit_Bayes_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 0), 
    L0.prior, 
    male.L50, 
    male.t50_Age1,
    step = c(0.047, 0.034, 0.025),
    k_based = TRUE,    
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )

### Age2 band count ----
Male_Age2_3par <- 
  fit_Bayes_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 0), 
    L0.prior, 
    male.L50, 
    male.t50_Age2,
    step = c(0.056, 0.041, 0.030),
    k_based = TRUE,
    L_infty.inform = FALSE,
    t_mat.inform = TRUE,
    LOOCV = TRUE,
    MH_chains = MH_chains,
    n_independence = n_independence,
    n_bootstraps = n_bootstraps,
    LOO_chains = LOO_chains
  )



#### Fit Frequentist Growth models ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

## Females -----------

### Age1 Band count ----
Female_frequentist_Age1 <-
  fit_frequentists_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 1),
    mean(L0.prior), 
    mean(fem.L50), 
    mean(fem.t50_Age1)
  )

### Age2 Band count ----
Female_frequentist_Age2 <-
  fit_frequentists_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 1),
    mean(L0.prior), 
    mean(fem.L50), 
    mean(fem.t50_Age2)
  )


## Males -------------

### Age1 Band count ----
Male_frequentist_Age1 <-
  fit_frequentists_growth(
    Age1, FL, 
    length_age.dat %>% 
      filter(Sex == 0),
    mean(L0.prior), 
    mean(male.L50), 
    mean(male.t50_Age1)
  )

### Age2 Band count ----
Male_frequentist_Age2 <-
  fit_frequentists_growth(
    Age2, FL, 
    length_age.dat %>% 
      filter(Sex == 0),
    mean(L0.prior), 
    mean(male.L50), 
    mean(male.t50_Age2)
  )

