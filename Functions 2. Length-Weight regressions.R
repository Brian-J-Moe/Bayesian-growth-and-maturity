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
#########################       Functions 2. Length and Weight Regressions       #########################
#########################                                                        #########################
##########################################################################################################
##########################################################################################################

##### Required Packages #####
## ~~~~~~~~~~~~~~~~~~~~~~~ ##

stats.packs <- 
  c("stats", "tidyverse") 

graphical.packs <- 
  c("ggplot2", "patchwork", "wesanderson")
  

pack.list <- c(stats.packs, graphical.packs)

invisible(lapply(pack.list, require, character.only = TRUE))


###########################################################################
###############             Auxilary functions              ###############
###########################################################################

### Significant Figures with correct decimals ###
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
sigfigs <- 
  function(x, digits = 2) {
    x <- x %>% as.matrix()
    new.x <- x
    dims <- dim(x)
    for (i in 1:dims[1]) {
      for (ii in 1:dims[2]) {
        if (is.na(x[i, ii])) {
          new.x[i, ii] <- NA
        } else {
          if (abs(x[i, ii]) > 1) {
            if (x[i,ii] > 0) {
              dig <- x[i, ii]
              new.x[i,ii] <- floor(dig) + signif(dig - floor(dig), digits)
            } else {
              dig <- x[i, ii]
              new.x[i,ii] <- ceiling(dig) + signif(dig - ceiling(dig), digits)
            }
          } else {
            new.x[i,ii] <- signif(x[i, ii], digits)
          } } } }
    return(new.x)
  }

##### Set plot theme #####
## ~~~~~~~~~~~~~~~~~~~~ ##

plot.theme <- 
  theme_bw(base_size = 8) %+replace%
  theme(
    plot.background = element_rect(color = NA),
    legend.background = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = 'black', size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(2.5,"line"),
    legend.key = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    plot.margin = unit(c(10,10,10,10), "pt"),
    aspect.ratio = 3/4
  )

Fem_col <- wesanderson::wes_palettes$FantasticFox1[3] 
Male_col <- wesanderson::wes_palettes$FantasticFox1[4]


##########################################################################
###############             Primary functions              ###############
##########################################################################

#### Length-Length Functions ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~ ###

Length_Length <- function(length.1, length.2, sex, data) {
  
  length.1 <- substitute(length.1)
  length.2 <- substitute(length.2)
  sex <- substitute(sex)
  
  temp.dat <- data.frame(Sex = data[[sex]], 
                         Length.1 = data[[length.1]],
                         Length.2 = data[[length.2]]) %>% na.omit()
  
  m1.fix <- glm(Length.2 ~ Length.1 + Sex, data = temp.dat)
  m1.null <- glm(Length.2 ~ Length.1, data = temp.dat)
  
  sum1.fix <- summary(m1.fix)
  sum1.null <- summary(m1.null)
  
  fixed = sum1.fix$coefficients[3,4] <= 0.05
  
  if (fixed != FALSE)
    m1 <- m1.fixed
  else
    m1 <- m1.null
  
  pred <- predict(m1, newdata = temp.dat, se.fit = TRUE)[1:2] %>%
    as.data.frame()
  
  pred <- cbind(temp.dat, pred) %>%
    dplyr::mutate(
      lwr = fit - se.fit * qnorm(0.975),
      upr = fit + se.fit * qnorm(0.975) )
  
  plot1 <- 
    pred %>%
    dplyr::mutate(
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males"))) %>%
    ggplot2::ggplot() +
    plot.theme +
    ggplot2::geom_point(aes(x = Length.1, y = Length.2, color = Sex, fill = Sex, shape = Sex),
                        size = 1.5, alpha = 0.5, stroke = 0.5) +
    ggplot2::geom_ribbon(aes(x = Length.1, y = fit, 
                             ymin = lwr, ymax = upr), alpha = 0.5) +
    ggplot2::geom_line(aes(x = Length.1, y = fit), 
                       linewidth = 1.15) +
    ggplot2::scale_shape_manual(values = c(21,24), name = NULL) +
    ggplot2::scale_fill_manual(values = c(Fem_col, Male_col), name = NULL) +
    ggplot2::scale_color_manual(values = c(Fem_col, Male_col), name = NULL) +
    ggplot2::guides(fill = guide_legend(override.aes = list(size = 5)))
  
  if (fixed != FALSE) {
    plot2 <- 
      plot1 + 
      ggplot2::guides(fill = 'none', shape = 'none',
                      color = 'none') +
      facet_grid(. ~ Sex)
    return(list(Models = list(Fixed.effects = m1.fix, 
                                   Null = m1.null),
                     Summaries = list(Fixed.effects = sum1.fix, 
                                      Null = sum1.null),
                Plot = plot2) )
  } else {
    return(list(Models = list(Fixed.effects = m1.fix, 
                                   Null = m1.null),
                     Summaries = list(Fixed.effects = sum1.fix, 
                                      Null = sum1.null),
                Plot = plot1) )
  }
} 

plot_LxL <- function(PCL, FL, TL, sex, data) {
  
  PCL <- substitute(PCL)
  FL <- substitute(FL)
  TL <- substitute(TL)
  sex <- substitute(sex)
  
  temp.dat <- data.frame(Sex = data[[sex]], 
                         PCL = data[[PCL]],
                         FL = data[[FL]],
                         TL = data[[TL]])
  
  axis.lims <- range(temp.dat[, 2:4], na.rm = TRUE)
  
  p1 <- Length_Length(FL, TL, Sex, temp.dat)$Plot +
    ggplot2::scale_x_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025))) + 
    theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = margin(r = 0, b = 0)) +
    ggplot2::scale_y_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025)),
                                name = "Total Length (cm)") 
  
  p2 <- Length_Length(FL, PCL, Sex, temp.dat)$Plot +
    ggplot2::scale_x_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025)), 
                                name = "Fork Length (cm)") +
    ggplot2::scale_y_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025)),
                                name = "Pre-Caudal Length (cm)") +
    theme(plot.margin = margin(r = 0, t = 0))
  
  p3 <- Length_Length(TL, PCL, Sex, temp.dat)$Plot +
    ggplot2::scale_x_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025)), 
                                name = "Total Length (cm)") +
    ggplot2::scale_y_continuous(limits = axis.lims,
                                expand = expansion(mult = c(0.025, 0.025))) +
    theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(l = 0, t = 0))
  
  p1 + guide_area() + p2 + p3 + 
    plot_layout(nrow = 2, guides = "collect") 
}


#### Length-Weight Functions ####
### ~~~~~~~~~~~~~~~~~~~~~~~~~ ###

Length_Weight <- function(lt, wt, sex, data) {
  
  lt <- substitute(lt)
  wt <- substitute(wt)
  sex <- substitute(sex)
  
  temp.dat <- data.frame(Sex = data[[sex]], 
                         Length = data[[lt]],
                         Weight = data[[wt]]) %>% na.omit()
  
  glm.fix <- glm(Weight ~ Length + Sex, data = temp.dat)
  glm.fix_quad <- glm(Weight ~ Length + I(Length^2) + Sex, data = temp.dat)
 
  glm.null <- glm(Weight ~ Length, data = temp.dat)
  glm.null_quad <- glm(Weight ~ Length + I(Length^2), data = temp.dat)
  
  AIC.fixed <- c(AIC(glm.fix), AIC(glm.fix_quad))
  AIC.null <- c(AIC(glm.null), AIC(glm.null_quad))
  
  AIC.results0 <- rbind(AIC.fixed, AIC.null)
  AIC.results <- AIC.results0 - min(AIC.results0)
  
  colnames(AIC.results) <- c("GLM", "Quad GLM")
  rownames(AIC.results) <- c("Fixed", "Null")
  
  model.index <- which(c(AIC.fixed, AIC.null) == min(c(AIC.fixed, AIC.null)) )
  
  mod.list <- list(glm.fix, glm.fix_quad, 
                   glm.null, glm.null_quad)
  
  opt.modl <- mod.list[[model.index]]
  
  if (model.index < 3) {
    fixed = TRUE
  } else { 
    fixed = FALSE
  }
  
  x.range0 <- range(temp.dat$Length[temp.dat$Sex == 0])
  x.range1 <- range(temp.dat$Length[temp.dat$Sex == 1])
  
  dat_fit0 <- data.frame(Length = seq(x.range0[1], x.range0[2], 
                                 by = diff(x.range0) / 1000),
                         Sex = 0)
  dat_fit1 <- data.frame(Length = seq(x.range1[1], x.range1[2], 
                                 by = diff(x.range1) / 1000),
                         Sex = 1)
  dat_fit <- rbind(dat_fit0, dat_fit1) %>% as.data.frame()
  
  pred <- predict(opt.modl, newdata = dat_fit, se.fit = TRUE)[1:2] %>%
    as.data.frame()
  
  pred <- cbind(dat_fit, pred) %>%
    dplyr::mutate(
      lwr = fit - se.fit * qnorm(0.975),
      upr = fit + se.fit * qnorm(0.975),
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males")))
  
  plot1 <- 
    temp.dat %>%
    dplyr::mutate(
      Sex = factor(Sex, 
                   levels = c(1, 0),
                   labels = c("Females", "Males"))) %>%
    ggplot2::ggplot() +
    plot.theme +
    ggplot2::geom_point(aes(x = Length, y = Weight, 
                            color = Sex, fill = Sex, shape = Sex),
                        size = 1.5, alpha = 0.5, stroke = 0.5) +
    ggplot2::geom_ribbon(data = pred, aes(x = Length, y = fit, 
                             ymin = lwr, ymax = upr), alpha = 0.5) +
    ggplot2::geom_line(data = pred, aes(x = Length, y = fit), 
                       linewidth = 1.15) +
    ggplot2::scale_shape_manual(values = c(21,24), name = NULL) +
    ggplot2::scale_fill_manual(values = c(Fem_col, Male_col), name = NULL) +
    ggplot2::scale_color_manual(values = c(Fem_col, Male_col), name = NULL) +
    ggplot2::guides(fill = guide_legend(override.aes = list(size = 3))) +
    ggplot2::xlab("Fork Length (cm)") +
    ggplot2::ylab("Body Weight (kg)") 
  
  if (fixed != FALSE) {
    plot2 <- 
      plot1 + 
      ggplot2::guides(fill = 'none', shape = 'none',
                      color = 'none') +
      facet_grid(. ~ Sex) +
      theme(
        strip.text.x = element_text(margin = margin(1, 0, 1, 0, "mm"))
        )
    return(list(Model = opt.modl,
                AIC.table = AIC.results,
                Plot = plot2) )
  } else {
    return(list(Model = opt.modl,
                AIC.table = AIC.results,
                Plot = plot1) )
  }
} 

