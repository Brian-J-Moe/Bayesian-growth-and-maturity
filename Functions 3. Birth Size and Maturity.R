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

#### Predict missing lengths of embryos ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
pred_embryo.lt <- function(length.x, length.y, data, lt.embryo) {
  length.x <- substitute(length.x)
  length.y <- substitute(length.y)
  
  temp.dat <- data.frame(x = data[[length.x]],
                         y = data[[length.y]]) %>% na.omit()
  
  m1 <- glm(y ~ x, data = temp.dat)
  
  embryo.dat <- data.frame(x = lt.embryo)
  predict(m1, newdata = embryo.dat)
}

## Metropolis-Hastings within Gibbs
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


get_x50_MH.Gibbs <- function(x0, x1, acceptance_rate = 0.24, iterations = 50000, chains = 4) {
  
  
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
  
  get_MH.Gibbs_priors <- function(x, y, data) {
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
    upper <- c(0, Inf)
    
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
  
  get_R_x50 <- function(x50_MH.Gibbs_object) {
    
    L <- nrow(x50_MH.Gibbs_object %>% filter(chain == 1))
    
    j <- x50_MH.Gibbs_object %>% pull(chain) %>% unique() %>% length()
    
    x.j <- 
      foreach(i = icount(j), .combine = rbind) %do% {
        test.diag <- x50_MH.Gibbs_object %>% filter(chain == i)
        
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
  
  summary_x50.Gibbs <- function(x50_MH.Gibbs_object) {
    
    x.trim <- x50_MH.Gibbs_object %>% filter(alpha == 1)
    
    data.frame(Mean = apply(x.trim[,3:6], 2, mean),
               Variance = apply(x.trim[,3:6], 2, var)) %>%
      dplyr::mutate(
        "2.5%" = Mean - sqrt(Variance) * qnorm(0.975),
        "97.5%" = Mean + sqrt(Variance) * qnorm(0.975)
      ) %>% round(4)
  }
  
  
  x0.data <- data.frame(x = x0, y = 0)
  
  x1.data <- data.frame(x = x1, y = 1)
  
  temp.dat <- rbind(x0.data, x1.data) %>% 
    as.data.frame() %>% na.omit()
  
  get.priors <- get_MH.Gibbs_priors(x, y, temp.dat)
  
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
      
      for (i in 1:(iterations - 1)) {
        
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
      
      sample_chain <- data.frame(chain = n, iteration = 1:length(x50), x50, 
                                 sample_beta, sigma = sigma, alpha = alpha_track)
      
      sample_chain[-burn_in, ]
    }
  
  stopCluster(clm)
  
  R.vals <- get_R_x50(sim.x50)
  
  r.t <- R.vals < 1.1
  
  r.true <- isTRUE(isTRUE(r.t[1] & r.t[2]))
  
  summary <- summary_x50.Gibbs(sim.x50)
  
  cat("\n\n")
  cat("Summary\n")
  print(summary)
  cat("\n\n\n")
  cat("Gelman-Rubin Convergence Diagnostic (R)\n\n")
  print(R.vals)
  cat("\n")
  cat(" Chains = "); cat(chains); cat("\n")
  if (r.true == TRUE) {
    cat(" Convergence Achieved")
  } else{
    cat(" No Convergence. Add more chains")
  }
  cat("\n\n\n")
  
  return(list(
    "Summary" = summary,
    "Gelman-Rubin R" = R.vals,
    "Raw Results" = sim.x50)
  )
  
}


## Metropolis-Hastings Independence Bootstrap
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_x50_MH.bootstrap <- function(x0, x1, x50_MH.Gibbs_object,
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
  
  x.trim <- x50_MH.Gibbs_object[[3]] %>% filter(alpha == 1)
  
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
  
  cat("\n\n")
  cat("Summary\n")
  print(summary)
  cat("\n\n\n")
  
  return(list(Summary = summary,
              "Raw Results" = raw.data))
}




### Maturity Ogives ####
### ~~~~~~~~~~~~~~~~~ ###
generate_ogive_data <- function(x, maturity, sex, data, MC_logit_x50.obj) {
  
  x <- substitute(x)
  maturity <- substitute(maturity)
  sex <- substitute(sex)
  
  temp.data <- data.frame(x = data[[x]],
                         y = data[[maturity]],
                         Sex = data[[sex]]) %>% 
    na.omit() %>%
    dplyr::mutate(
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males")))
  
  x50.female <- MC_logit_x50.obj$Females$Bootstrap.Chain$x50.Results[3, ]
  x50.male <- MC_logit_x50.obj$Males$Bootstrap.Chain$x50.Results[3, ]
  
  x50_data <- data.frame(
    Sex = c(1, 0), rbind(x50.female, x50.male)) %>%
    dplyr::mutate(
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males")))
  
  beta.female <- MC_logit_x50.obj$Females$Bootstrap.Chain$x50.Results[c('alpha', 'beta'), 'Mean']
  beta.male <- MC_logit_x50.obj$Males$Bootstrap.Chain$x50.Results[c('alpha', 'beta'), 'Mean']
  
  range_x <- rbind(
    range(temp.data[temp.data$Sex == "Females", "x"]),
    range(temp.data[temp.data$Sex == "Males", "x"]))
  
  female_seq <- seq(range_x[1, 1], range_x[1, 2], by = diff(range_x[1, ]) / 1000) 
  
  male_seq <- seq(range_x[2, 1], range_x[2, 2], by = diff(range_x[2, ]) / 1000) 
  
  
  beta.female <- MC_logit_x50.obj$Females$Bootstrap.Chain$Raw.Parameters[, -3] 
  
  beta.male <- MC_logit_x50.obj$Males$Bootstrap.Chain$Raw.Parameters[, -3]
  
  n <- nrow(beta.female)
  
  for (i in 1:n) {
    
    temp.beta1 <- beta.female[i, ] %>% as.numeric()
    temp.beta2 <- beta.male[i, ] %>% as.numeric()
    
    temp.pred_fem = 1 / (1 + exp(-cbind(1, female_seq) %*% temp.beta1) )
    temp.pred_male = 1 / (1 + exp(-cbind(1, male_seq) %*% temp.beta2) )
    
    if (i == 1) {
      pred.fem <- temp.pred_fem
      pred.male <- temp.pred_male
    } else {
      pred.fem <- cbind(pred.fem, temp.pred_fem)
      pred.male <- cbind(pred.male, temp.pred_male)
    }
  }
  
  pred.fem <- data.frame(Sex = 1, x = female_seq, fit = rowMeans(pred.fem))
  pred.male <- data.frame(Sex = 0, x = male_seq, fit = rowMeans(pred.male))
  
  
  pred_data <- rbind(pred.fem, pred.male) %>% as.data.frame() %>%
    dplyr::mutate(
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males")))
  
  return(list(fit = pred_data, x50 = x50_data) )
}


##### Gibbsings Posterior plot #####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

plot_MH.Gibbs <- 
  function(x50_MH.Gibbs_object) {
    
    plot.theme <- 
      theme_bw(base_size = 8) %+replace%
      theme(
        plot.background = element_rect(color = NA),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color = 'black', size = 10),
        strip.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(2.5,"line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        plot.margin = unit(c(10,10,10,10), "pt"),
        aspect.ratio = 0.5
      )
    
    plot.data <- x50_MH.Gibbs_object$`Raw Results`
    
    beta0.dens <- 
      foreach(n = icount(plot.data %>% pull(chain) %>% unique() %>% max()), 
              .combine = rbind) %do% {
                
                sub.dat <- plot.data %>% filter(chain == n & alpha == 1)
                
                beta0.dens <- density(sub.dat %>% pull(beta0))
                
                point.beta0 <- data.frame(chain = n, 
                                          x = beta0.dens$x,
                                          y = beta0.dens$y)
              } %>% as.data.frame()
    
    beta1.dens <- 
      foreach(n = icount(plot.data %>% pull(chain) %>% unique() %>% max()), 
              .combine = rbind) %do% {
                
                sub.dat <- plot.data %>% filter(chain == n & alpha == 1)
                
                beta1.dens <- density(sub.dat %>% pull(beta1))
                
                point.beta1 <- data.frame(chain = n, 
                                          x = beta1.dens$x,
                                          y = beta1.dens$y)
              } %>% as.data.frame()
    
    beta0_seq <- 
      seq(min(plot.data %>% 
                filter(alpha == 1) %>% 
                pull(beta0)), 
          max(plot.data %>% 
                filter(alpha == 1) %>% 
                pull(beta0)), 
          by = 0.01)
    
    beta1_seq <- 
      seq(min(plot.data %>% 
                filter(alpha == 1) %>% 
                pull(beta1)), 
          max(plot.data %>% 
                filter(alpha == 1) %>% 
                pull(beta1)), 
          by = 0.01)
        
    gpars.beta0 <- 
      mlgamma(-(plot.data %>% 
                  filter(alpha == 1) %>% 
                  pull(beta0))) %>% 
      t() %>% as.data.frame()
    
    gpars.beta1 <- 
      mlgamma(plot.data %>% 
                filter(alpha == 1) %>% 
                pull(beta1)) %>% 
      t() %>% as.data.frame()
    
    dens.gamma_beta0 <- 
      data.frame(
        x = -beta0_seq,
        y = dgamma(-beta0_seq,
                   shape = gpars.beta0$shape,
                   scale = gpars.beta0$scale)) %>%
      mutate(
        x = -x
      )
    dens.gamma_beta1 <- 
      data.frame(
        x = beta1_seq,
        y = dgamma(beta1_seq,
                   shape = gpars.beta1$shape,
                   scale = gpars.beta1$scale))
    
    (beta0.dens.plot <-
        plot.data %>%
        ggplot() +
        plot.theme +
        geom_density(aes(x = beta0, color = as.factor(chain) ),
                     adjust = 1, linewidth = 0.75) +
        geom_line(data = dens.gamma_beta0, 
                  aes(x = x, y = y), linewidth = 0.75) +
        ggtitle(expression(paste(beta[0], " Density"))) + 
        scale_color_discrete(guide = "none")
    )
    
    (beta1.dens.plot <-
        plot.data %>%
        ggplot() +
        plot.theme +
        geom_density(aes(x = beta1, color = as.factor(chain) ), 
                     adjust = 1, linewidth = 0.75) +
        geom_line(data = dens.gamma_beta1, 
                  aes(x = x, y = y), linewidth = 0.75) +
        ggtitle(expression(paste(beta[1], " Density"))) +
        scale_color_discrete(guide = "none")
    )
    
    (beta0.track.plot <-
        plot.data %>%
        ggplot(aes(y = beta0, x = iteration, color = as.factor(chain) )) +
        plot.theme +
        geom_line(linewidth = 0.75) +
        ggtitle(expression(paste(beta[0], " Trace"))) +
        scale_color_discrete(guide = "none")
    )
    
    (beta1.track.plot <-
        plot.data %>%
        ggplot(aes(y = beta1, x = iteration, color = as.factor(chain) )) +
        plot.theme +
        geom_line(linewidth = 0.75) +
        ggtitle(expression(paste(beta[1], " Trace"))) +
        scale_color_discrete(guide = "none")
    )
    
    beta0.track.plot + beta0.dens.plot + beta1.track.plot + beta1.dens.plot + plot_layout(ncol = 2) 
  }

