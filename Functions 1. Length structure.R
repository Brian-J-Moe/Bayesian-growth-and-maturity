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
#########################             Functions 1. Length Structure              #########################
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

##### Length Structure #####
## ~~~~~~~~~~~~~~~~~~~~~~ ##

length_structure <- 
  function(lt, sex, data) {
    
    lt <- substitute(lt)
    sex <- substitute(sex)
    new.dat <- cbind(data[[lt]], data[[sex]]) %>% na.omit()
    
    f <- new.dat[new.dat[,2] == 1, 1]
    m <- new.dat[new.dat[,2] == 0, 1]
    
    n.catch <- c(length(f), length(m)) 
    
    lt.range <- rbind(range(f), range(m))
    
    mean <- round(c(mean(f), mean(m)), 2)
    st.d <- round(c(sd(f), sd(m)), 2)
    
    shapiro <-  
      c(sigfigs(stats::shapiro.test(f)$p.value, 3),
        sigfigs(stats::shapiro.test(m)$p.value, 3) )
    
    lt_stats <- 
      cbind(n.catch, mean, st.d, lt.range, shapiro)
    
    colnames(lt_stats) <- c("n", "Mean", "SD", "Min", "Max", "Shapiro-Wilk")
    rownames(lt_stats) <- c("Females", "Males")
    
    return(lt_stats)
  }

##### Density and QQ plots #####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

density_plot <- 
  function(lt, sex, data) {
    
    lt <- substitute(lt)
    sex <- substitute(sex)
    
    temp.dat <- data.frame(Sex = data[[sex]], 
                           Length = data[[lt]],
                           type = "Length Distribution ") %>% 
      na.omit() %>%
      dplyr::mutate(
        Sex = factor(Sex, 
                     levels = c(1, 0),
                     labels = c("Females", "Males")))
    
    temp.dat %>%
      ggplot2::ggplot(aes(x = Length)) +
      plot.theme +
      facet_grid(type ~ Sex, switch = "x") +
      
      ggplot2::geom_density(aes(color = Sex, fill = Sex), 
                            linewidth = 1.15, alpha = 0.25) +
      
      ggplot2::scale_color_manual(
        values = c(Fem_col, Male_col), name = NULL) +
      ggplot2::scale_fill_manual(
        values = c(Fem_col, Male_col), name = NULL) +
      
      ggplot2::scale_x_continuous(
        expand = expansion(mult = c(0.025, 0.025)),
        name = "Fork Length (cm)", position = "top") +
      ggplot2::scale_y_continuous(breaks = seq(0.025, .1, by = 0.025),
        expand = expansion(mult = c(0.025, 0.025)),
        name = "Density") +
      guides(color = "none", fill = "none") +
      theme(strip.text.x = element_text(margin = margin(1, 0, 1, 0, "mm")),
        strip.text.y = element_text(margin = margin(0, 1, 0, 1, "mm")))
  }

qq_plot <- 
  function(lt, sex, data) {
    
    lt <- substitute(lt)
    sex <- substitute(sex)
    
    temp.dat <- data.frame(Sex = data[[sex]], 
                           Length = data[[lt]],
                           type = "Q-Q Plot") %>% 
      na.omit() %>%
      dplyr::mutate(
        Sex = factor(Sex, 
                     levels = c(1, 0),
                     labels = c("Females", "Males")))
    
    temp.dat %>%
       ggplot(aes(sample = Length)) +
  plot.theme +
  
  facet_grid(type ~ Sex) +
  geom_qq(aes(fill = Sex, color = Sex, shape = Sex), 
          size = 1.5, stroke = 0.5, alpha = 0.5) +
      
      geom_qq_line(linewidth = 1.15) +
      
      ggplot2::scale_color_manual(
        values = c(Fem_col, Male_col), name = NULL) +
      
      ggplot2::scale_fill_manual(
        values = c(Fem_col, Male_col), name = NULL) +
      
      ggplot2::scale_shape_manual(
        values = c(21, 24), name = NULL) +
      
      ggplot2::scale_x_continuous(
        expand = expansion(mult = c(0.025, 0.025)), 
        name = "Theoretical Quantiles") +
      
      ggplot2::scale_y_continuous(breaks = seq(20,120, by = 30),
        expand = expansion(mult = c(0.025, 0.025)), 
        name = "Fork Length (cm)") +
      
      guides(color = "none", fill = "none", shape = "none") +
      theme(strip.text.y = element_text(margin = margin(0, 1, 0, 1, "mm")),
            strip.text.x = element_text(margin = margin(1, 0, 1, 0, "mm")))
  } 

density_qq_plot <- 
  function(lt, sex, data) {
    
    lt <- substitute(lt)
    sex <- substitute(sex)
    
    temp.dat <- data.frame(Sex = data[[sex]], 
                           Length = data[[lt]]) %>% 
      na.omit() 
    
    annotation1 <- data.frame(text = c("A", "B"), Sex = c(0, 1)) %>%
      dplyr::mutate(
        Sex = factor(Sex,
                     levels = c(1, 0),
                     labels = c("Females", "Males"))
      )
    
    annotation2 <- data.frame(text = c("C", "D"), Sex = c(0, 1)) %>%
      dplyr::mutate(
        Sex = factor(Sex,
                     levels = c(1, 0),
                     labels = c("Females", "Males"))
      )
    
    gg_dense <- density_plot(Length, Sex, temp.dat) 
    
    #  geom_text(data = annotation1, aes(x = -Inf, y = Inf, 
     #                                   hjust = -0.5, vjust = 1.5,
      #                                  label = text),
       #         size = 8)
    
    gg_qq <- qq_plot(Length, Sex, temp.dat) 
    
    new_dense <- gg_dense + theme(plot.margin = margin(b = 0))   
    new_qq <- gg_qq + theme(plot.margin = margin(t = 0), 
                                    strip.background.x = element_blank(),
                                    strip.text.x = element_blank() )
    
    new_dense / new_qq
#    gg_qq / plot_spacer() / gg_density +
#      plot_layout(heights = c(5, -1, 5))
    
    
    
  }

