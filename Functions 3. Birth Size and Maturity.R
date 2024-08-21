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
#########################           Functions 3. Birth Size and Maturity         #########################
#########################                                                        #########################
##########################################################################################################
##########################################################################################################

stats.packs <- 
  c("stats", "tidyverse", "mvtnorm", "tmvtnorm", "mgcv", "dglm")  

graphical.packs <- 
  c("ggplot2", "patchwork", "wesanderson")

pack.list <- c(stats.packs, graphical.packs)
invisible(lapply(pack.list, require, character.only = TRUE))

## Metropolis-Hastings within Gibbs
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


get_x50_Metropolis.Hastings <- function(x0, x1, acceptance_rate = 0.24, iterations = 9999, chains = 6) {
  
  mlgamma <- function(x) {
    fit <- dglm(
      x ~ 1, 
      family = Gamma(link = "log"), 
      mustart = mean(x)
    )
    mu <- exp(fit$coefficients)
    shape <- exp(-fit$dispersion.fit$coefficients)
    scale <- mu/shape
    result <- c(shape, scale)
    names(result) <- c("shape", "scale")
    result
  }
  
  rtnorm <- 
    function(n, mu, sd, lower = -Inf, upper = Inf) {
      
      p_lower <- pnorm(lower, mu, sd)
      p_upper <- pnorm(upper, mu, sd)
      
      qnorm(runif(n, p_lower, p_upper), mu, sd)
    }
  
  get_Metropolis.Hastings_priors <- function(x, y, data) {
    x <- substitute(x)
    y <- substitute(y)
    
    input.data <- data.frame(x = data[[x]],
                             y = data[[y]]) %>% na.omit()
    
    x0 <- input.data$x[input.data$y == 0]
    x1 <- input.data$x[input.data$y == 1]
    
    ref_x <- mean(c(max(x0), min(x1)) )
    
    beta_prior <- 1 - (ref_x/max(x1))
    alpha_prior <- -ref_x * beta_prior
    priors <- c(alpha_prior, beta_prior)
    
    lower <- c(-Inf, 0)
    upper <- c(0, 1)
    
    prior_data.frame <- 
      data.frame(Prior = priors, Lower.lim = lower, Upper.lim = upper)
    rownames(prior_data.frame) <- c("alpha", "beta")
    
    return(prior_data.frame)
  }
  
  logit_log_likelihood <- function(x, y, data, beta, gamma.shape = NA, gamma.scale = NA,
                                   informative = FALSE) {
    x <- substitute(x)
    y <- substitute(y)
    
    temp.dat <- data.frame(x = data[[x]],
                           y = data[[y]]) %>% na.omit()
    
    x <- cbind(1, temp.dat$x)
    y <- temp.dat$y
    
    p <- 1 / (1 + exp(-x %*% beta))
    
    if (informative != FALSE) {
      
      l <- sum(dbinom(y, 1, p, log = TRUE))  + 
        dgamma(-beta[1], gamma.shape[1], scale = gamma.scale[1], log = T) +
        dgamma(beta[2], gamma.shape[2], scale = gamma.scale[2], log = T)
    } else {
      l <- sum(dbinom(y, 1, p, log = TRUE))
    }
    return(l)
  }
  
  get_R_x50 <- function(x50.MetHast_object) {
    
    L <- nrow(x50.MetHast_object %>% filter(chain == 1))
    
    j <- x50.MetHast_object %>% pull(chain) %>% unique() %>% length()
    
    x.j <- 
      foreach(i = icount(j), .combine = rbind) %do% {
        test.diag <- x50.MetHast_object %>% filter(chain == i)
        
        beta_accepted <- test.diag %>% filter(alpha == 1)
        
        beta0.mu <- beta_accepted %>% pull(beta0) %>% mean()
        beta0.var <- beta_accepted %>% pull(beta0) %>% var()
        
        beta1.mu <- beta_accepted %>% pull(beta1) %>% mean()
        beta1.var <- beta_accepted %>% pull(beta1) %>% var()
        
        
        data.frame(beta0_mu = beta0.mu,
                   beta0_var = beta0.var,
                   beta1_mu = beta1.mu,
                   beta1_var = beta1.var)
        
      }
    
    x.star_beta0 <- x.j[, 1] %>% mean()
    
    B_beta0 <- (L / (j - 1)) * sum((x.j[, 1] - x.star_beta0)^2)
    
    W_beta0 <- sum(x.j[, 2])/j
    
    beta0_R <- 
      (((L - 1) / L * W_beta0) + ((1 / L) * B_beta0)) / W_beta0
    
    x.star_beta1 <- x.j[, 3] %>% mean()
    
    B_beta1 <- (L / (j - 1)) * sum((x.j[, 3] - x.star_beta1)^2)
    
    W_beta1 <- sum(x.j[, 4])/j
    
    beta1_R <- 
      (((L - 1) / L * W_beta1) + ((1 / L) * B_beta1)) / W_beta1
    
    r.dat <- data.frame(beta0 = beta0_R, beta1 = beta1_R)
    rownames(r.dat) <- "R"
    
    return(r.dat %>% round(4))
  }
  
  summary_x50.MetHast <- function(x50.MetHast_object) {
    
    x.trim <- x50.MetHast_object %>% filter(alpha == 1)
    
    data.frame(Mean = apply(x.trim[,2:5], 2, mean),
               Variance = apply(x.trim[,2:5], 2, var)) %>%
      dplyr::mutate(
        "2.5%" = Mean - sqrt(Variance) * qnorm(0.975),
        "97.5%" = Mean + sqrt(Variance) * qnorm(0.975)
      ) %>% round(4)
  }
  
  
  x0.data <- data.frame(x = x0, y = 0)
  
  x1.data <- data.frame(x = x1, y = 1)
  
  temp.dat <- rbind(x0.data, x1.data) %>% 
    as.data.frame() %>% na.omit()
  
  get.priors <- get_Metropolis.Hastings_priors(x, y, temp.dat)
  
  prior <- get.priors[, 1] %>% as.numeric()
  lower <- get.priors[, 2] %>% as.numeric()
  upper <- get.priors[, 3] %>% as.numeric()
  
  library(progress) %>% suppressWarnings()
  
  pb <- progress_bar$new(
    format = "Chains complete = :letter [:bar] :elapsed | eta: :eta",
    total = chains,
    width = 120
  )
  
  progress_count <- 1:chains
  
  progress <- function(n) {
    pb$tick(tokens = list(letter = progress_count[n])) 
  }
  
  opts <- list(progress = progress)
  
  library(foreach) %>% invisible()
  library(doParallel) %>% invisible()
  library(doSNOW) %>% invisible() %>% suppressWarnings()
  
  nCores <- detectCores()
  clm <- makeCluster(nCores)
  registerDoSNOW(clm) 
  
  burn_in <- seq(1, floor(iterations * 0.1))
  
  set.seed(1234)
  
  seed.vals <- sample(1:100000, chains, replace = F)
  
  sim.x50 <- foreach(
    n = icount(chains), .combine = rbind, .options.snow = opts) %dopar% 
    {
      library(tidyverse)
      
      set.seed(seed.vals[n])
      
      sample_beta <- beta <- prior
      
      sigma <- 0.15
      
      z <- 1
      alpha_track <- NA
      
      for (i in 1:iterations) {
        
        beta_proposed <- rtnorm(2, beta, abs(prior * sigma[i]), lower, upper)
        
        # Calculate acceptance ratio
        g_beta2 <- 
          logit_log_likelihood(x, y, temp.dat, beta_proposed)
        
        g_beta1 <- 
          logit_log_likelihood(x, y, temp.dat, beta) 
        
        R <- exp(g_beta2 - g_beta1)
        
        alpha <- rbinom(1, 1, min(1, R))
        
        if (alpha == 1) {
          beta <- beta_proposed
          z <- z + 1
        }
        
        sample_beta <- rbind(sample_beta, beta) %>% as.data.frame()
        
        alpha_track <- c(alpha_track, alpha)
        
        m <- z / (i + 1)
        tau <- sigma[i] * m / acceptance_rate
        
        sigma <- c(sigma, tau)
      }
      
      names(sample_beta) <- c("beta0", "beta1")
      
      x50 <- -sample_beta[, 1] / sample_beta[, 2]
      
      sample_chain <- data.frame(chain = n, x50, sample_beta, sigma = sigma, alpha = alpha_track)
      
      sample_chain[-burn_in, ]
    }
  
  stopCluster(clm)
  
  R.vals <- get_R_x50(sim.x50)
  
  summary <- summary_x50.MetHast(sim.x50)
  
  cat("\n")
  cat("Summary\n")
  print(summary)
  cat("\n\n")
  cat("Gelman-Rubin R\n")
  print(R.vals)
  cat("\n\n")
  
  return(list(
    "Summary" = summary,
    "Gelman-Rubin R" = R.vals,
    "Raw Results" = sim.x50)
  )
  
}


## Metropolis-Hastings Independence Bootstrap
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_x50_bootstrap <- function(x0, x1, x50.MetHast_object,
                              bootstraps = 999, 
                              indep_chains = 9999) {
  
  logit_log_likelihood <- function(x, y, data, beta, gamma.shape = NA, gamma.scale = NA,
                                   informative = FALSE) {
    x <- substitute(x)
    y <- substitute(y)
    
    temp.dat <- data.frame(x = data[[x]],
                           y = data[[y]]) %>% na.omit()
    
    x <- cbind(1, temp.dat$x)
    y <- temp.dat$y
    
    p <- 1 / (1 + exp(-x %*% beta))
    
    if (informative != FALSE) {
      
      l <- sum(dbinom(y, 1, p, log = TRUE))  + 
        dgamma(-beta[1], gamma.shape[1], scale = gamma.scale[1], log = T) +
        dgamma(beta[2], gamma.shape[2], scale = gamma.scale[2], log = T)
    } else {
      l <- sum(dbinom(y, 1, p, log = TRUE))
    }
    return(l)
  }
  
  summary_x50.bootstrap <- function(x50.bootstrap_object) {
  
    data.frame(Mean = apply(x50.bootstrap_object, 2, mean),
               Variance = apply(x50.bootstrap_object, 2, var)) %>%
      dplyr::mutate(
        "2.5%" = Mean - sqrt(Variance) * qnorm(0.975),
        "97.5%" = Mean + sqrt(Variance) * qnorm(0.975)
      ) %>% round(4)
  }
  
  x0.data <- data.frame(x = x0, y = 0)
  
  x1.data <- data.frame(x = x1, y = 1)
  
  temp.dat <- rbind(x0.data, x1.data) %>% 
    as.data.frame() %>% na.omit()
  
  x.trim <- x50.MetHast_object[[3]] %>% filter(alpha == 1)
  
  beta0 <- x.trim %>% pull(beta0)
  
  beta1 <- x.trim %>% pull(beta1)
  
  d.beta0 <- mlgamma(-beta0)
  d.beta1 <- mlgamma(beta1)
  
  gamma.shape <- c(d.beta0[1], d.beta1[1])
  gamma.scale <- c(d.beta0[2], d.beta1[2])
  
  lower <- cbind(beta0, beta1) %>% apply(2, min)
  upper <- cbind(beta0, beta1) %>% apply(2, max)
  
  library(progress) %>% suppressWarnings()
  
  pb <- progress_bar$new(
    format = "Bootstraps complete = :letter [:bar] :elapsed | eta: :eta",
    total = bootstraps,
    width = 120
  )
  
  progress_count <- 1:bootstraps
  
  progress <- function(n) {
    pb$tick(tokens = list(letter = progress_count[n])) 
  }
  
  opts <- list(progress = progress)
  
  nCores <- detectCores()
  clb <- makeCluster(nCores)
  registerDoSNOW(clb) 
  
  set.seed(1234)
  
  seed.vals <- sample(1:1000000, bootstraps, replace = F)
  
  sim.boot <- 
    foreach(
      n = icount(bootstraps), .combine = rbind,
      .options.snow = opts
    ) %dopar% 
    {
      
      library(tidyverse)
      
      set.seed(seed.vals[n])
      
      boot.index <- sample(
        rownames(temp.dat), 
        nrow(temp.dat),
        replace = T)
      
      boot.dat <- temp.dat[boot.index, ]
      
      library(foreach) 
      library(doParallel) 
      library(doSNOW) %>% suppressWarnings()
      
      beta_proposed <- 
        foreach(icount(indep_chains), .combine = rbind) %dopar%
        {
          library(tidyverse)
          beta_test <- runif(2, lower, upper)
          
          temp.LL <- logit_log_likelihood(x, y, boot.dat, beta_test, 
                                          gamma.shape, gamma.scale,
                                          informative = TRUE)
          
          c(beta_test, temp.LL)
        } 
      
      wi <- exp(beta_proposed[, 3]) / 
        sum(exp(beta_proposed[, 3]))
      
      LL_index <- sample(1:length(wi), 1, prob = wi)
      
      beta_proposed[LL_index, -3]
    }
  
  stopCluster(clb)
  
  x50.stat <- -sim.boot[, 1] / sim.boot[, 2]
  
  raw.data <- data.frame(x50 = x50.stat, 
                         beta0 = sim.boot[, 1], 
                         beta1 = sim.boot[, 2]) %>% 
    `row.names<-`(NULL)
  
  summary <- summary_x50.bootstrap(raw.data)
  
  cat("\n")
  cat("Summary\n")
  print(summary)
  cat("\n\n")
  
  return(list(Summary = summary,
              "Raw Results" = raw.data))
}




