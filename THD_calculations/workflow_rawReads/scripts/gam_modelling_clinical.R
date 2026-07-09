library(data.table)
library(mgcv)
library(tidyverse)

library(marginaleffects)




setwd("~/Desktop/THD_calculations_local/nextclade/")

files = list.files("./", 
           pattern = "index_time_parameters_ww_time_full.tsv",
           recursive = TRUE,
           full.names = TRUE)


dts = files %>% map_dfr(~ fread(.x)) %>% mutate(date = date_decimal(time))


dt_byAmp = dts %>% 
  group_by(subclade,`legacy.clade`,timescale,time_window,amplicon,date) 

dt_acrossAmp = dts %>% 
  filter(time_window == 0.019 & timescale == 0.4) 





palette_colors <- as.vector(pals::cols25(14))
subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
               "C.1.9.4","D","D.1","D.3.1","D.5",
               "C.1","C.1.1","C.1.8","D.2","D.3")
custom_palette <- setNames(palette_colors, subclades)


#######ACROSS AMPLICON ##########

dt_plot = dts %>% 
  filter(time_window == 0.019 & timescale == 0.4) %>% 
  filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE)



seq_included = dt_plot %>%  pull(seqName) %>% unique()

########### OVERALL ######
ggplot(dt_plot, 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("Per sequence fitness estimate")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(dt_plot, 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))





########### D.3.1 ######

#BY AMPLICON

ggplot(dt_plot %>% filter(subclade == "D.3.1"), 
           aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))
#Overall

ggplot(dt_plot %>% filter(subclade == "D.3.1"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))




########### C.1.9 ######

#BY AMPLICON

ggplot(dt_plot %>% filter(subclade == "C.1.9"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Overall
ggplot(dt_plot %>% filter(subclade == "C.1.9"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))


########### C.1.9.3 ######

#BY AMPLICON

ggplot(dt_plot %>% filter(subclade == "C.1.9.3"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Overall

ggplot(dt_plot %>% filter(subclade == "C.1.9.3"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))






#####GAM MODELLING #####
amp3_weight = 0.3

dt_plot = dt_plot %>%
  mutate(amp_weighted = case_when(
  # amplicon == "amplicon3" ~ amp3_weight,
  # .default = (1-amp3_weight)/4))
    amplicon == "amplicon3" ~ 1,
    .default = 1))

gam_D.3 <- bam(indicies ~  
                 # s(time, bs = 'cr')+
                 # s(time, by = as.factor(amplicon), bs="cr", k = -1) +
                 # as.factor(amplicon),

                 s(time, bs = 'cs') + s(time, amplicon, bs="fs", k = 5),
               
               gamma = 5,
               
               weights = amp_weighted,
               data = dt_plot %>% filter(subclade == "D.3") %>% mutate(amplicon = as.factor(amplicon)),
               family = gaussian(link = identity),
               method = "fREML"
               
)


gam_D.3.1 <- bam(indicies ~  
                # s(time, bs = 'cr')+
                # s(time, by = as.factor(amplicon), bs="cr", k = -1) +
                # as.factor(amplicon),

                s(time, bs = 'cs') + s(time, amplicon, bs="fs", k = 5),
                
                gamma = 5,
                
                weights = amp_weighted,
                data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                family = gaussian(link = identity),
                method = "fREML"
                
)




gam_C.1 <- bam(indicies ~  
                   # s(time, bs = 'cs')+
                   # s(time, by = as.factor(amplicon), bs="cs", k = -1) +
                   # as.factor(amplicon),
               

                  s(time, bs = 'cs') + s(time, amplicon, bs="fs", k = 5),
                 
                  gamma = 5,
                 
                 weights = amp_weighted,
                 data = dt_plot %>% filter(subclade == "C.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)


gam_C.1.9 <- bam(indicies ~  
                # s(time, bs = 'cs')+
                # s(time, by = as.factor(amplicon), bs="cs", k = -1) +
                # as.factor(amplicon),
                
                s(time, bs = 'cs') + s(time, amplicon, bs="fs", k = 5),
                
                gamma = 5,
                
                weights = amp_weighted,
                data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                family = gaussian(link = identity),
                method = "fREML"
                
                )




gam_C.1.9.3 <- bam(indicies ~  
                  # s(time, bs = 'cr')+
                  # s(time, by = as.factor(amplicon), bs="cr", k = -1) +
                  # as.factor(amplicon),

                  s(time, bs = 'cs') + s(time, amplicon, bs="fs", k = 5),
                  
                  gamma = 5,
                
                 weights = amp_weighted,
                 data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
)







gam <- gam_C.1.9

AIC(gam)



summary(gam)



gam.check(gam)






plot_predictions(gam_C.1.9,
                 condition = list("time", "amplicon"),
                 points = 0.5,
                 newdata = dt_plot %>% filter(subclade == "C.1.9"))+
  coord_cartesian(ylim = c(-1, 1))+
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
    rows = vars(subclade),
    labeller = label_both) 



plot_predictions(gam_C.1.9,
                 condition = list("time"),
                 points = 0.5,
                 newdata = dt_plot %>% filter(subclade == "C.1.9"))+
  # coord_cartesian(ylim = c(-2, 3))+
  coord_cartesian(ylim = c(-1, 1))+
  scale_fill_manual(values = custom_palette) +
  facet_grid(
            # cols = vars(amplicon),
             rows = vars(subclade),
             labeller = label_both) 









#####COMBINED APPROACH (no weighted averaging needed) #########

#### Generating Overall Predictions (Replaces Weighted Average) #######

# By creating a grid of only 'time' values, we ask marginaleffects to
# produce the "population-level" prediction, averaged over the random 
# effects of 'amplicon'. This directly gives us the smoothed overall trend.

# 1. Define a time grid for each subclade
time_grid_D.3 <- dt_plot %>% filter(subclade == "D.3") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_D.3.1 <- dt_plot %>% filter(subclade == "D.3.1") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1<- dt_plot %>% filter(subclade == "C.1") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9 <- dt_plot %>% filter(subclade == "C.1.9") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9.3 <- dt_plot %>% filter(subclade == "C.1.9.3") %>% with(seq(min(time), max(time), length.out = 200))



overall_pred_D.3 <- predictions(
  gam_D.3,
  newdata = datagrid(
    time     = time_grid_D.3,
    amplicon = levels(dt_plot %>% filter(subclade == "D.3") %>% 
                        mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_D.3
  ),
  by = "time"
) %>% mutate(subclade = "D.3")


overall_pred_D.3.1 <- predictions(
  gam_D.3.1,
  newdata = datagrid(
    time     = time_grid_D.3.1,
    amplicon = levels(dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_D.3.1
  ),
  by = "time"
) %>% mutate(subclade = "D.3.1")


overall_pred_C.1 <- predictions(
  gam_C.1,
  newdata = datagrid(
    time     = time_grid_C.1,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1
  ),
  by = "time"
) %>% mutate(subclade = "C.1")


overall_pred_C.1.9 <- predictions(
  gam_C.1.9,
  newdata = datagrid(
    time     = time_grid_C.1.9,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9")


overall_pred_C.1.9.3 <- predictions(
  gam_C.1.9.3,
  newdata = datagrid(
    time     = time_grid_C.1.9.3,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9.3
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9.3")





# 3. Combine into a single data frame
overall_predictions <- bind_rows(
  overall_pred_D.3,
  overall_pred_D.3.1,
  overall_pred_C.1,
  overall_pred_C.1.9,
  overall_pred_C.1.9.3
) 


# 4. Plot the results
results <- ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  # geom_point(data = dt_plot %>% filter(subclade %in% c("D.3.1","D.3","C.1.9.3","C.1.9","C.1"))
  #            , aes(x = time,  y = indicies, colour = "grey"),alpha = 0.1)+
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(subclade), labeller = label_both) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()

# Display the plot
print(results)


ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()


######FULL LENGTH SEQUENCES MODELLING#### 

file = "results_gisaid_2526/index_time_parameters_ww_time_fullLength.tsv"


dt = file %>% map_dfr(~ fread(.x)) %>% mutate(date = date_decimal(time))

dt_plot = dt %>% filter(seqName %in% seq_included)
  # filter(time_window == 0.019 & timescale == 0.4) %>% 
  # filter(date <= as.Date("2025-03-15")) %>%
  # filter(is.na(subclade) == FALSE)



#no amplicon term now




gam_D.3.1 <- bam(indicies ~  
                   
                  s(time, bs = 'cs') ,
                 
  
                 data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)





gam_C.1.9 <- bam(indicies ~  
                   
                   s(time, bs = 'cs'),
                 data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)




gam_C.1.9.3 <- bam(indicies ~  
                    
                     
                     s(time, bs = 'cs'),
                   
                   data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                   family = gaussian(link = identity),
                   method = "fREML"
)







gam <- gam_D.3.1

AIC(gam)



summary(gam)



gam.check(gam)





plot_predictions(gam,
                 condition = list("time"),
                 points = 0.5,
                 newdata = dt_plot %>% filter(subclade == "C.1.9"))+
  scale_fill_manual(values = custom_palette) +
  facet_grid(
    # cols = vars(amplicon),
    rows = vars(subclade),
    labeller = label_both) 
 




time_grid_D.3.1 <- dt_plot %>% filter(subclade == "D.3.1") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9 <- dt_plot %>% filter(subclade == "C.1.9") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9.3 <- dt_plot %>% filter(subclade == "C.1.9.3") %>% with(seq(min(time), max(time), length.out = 200))





overall_pred_D.3.1 <- predictions(
  gam_D.3.1,
  newdata = datagrid(
    time     = time_grid_D.3.1,
    amplicon = levels(dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_D.3.1
  ),
  by = "time"
) %>% mutate(subclade = "D.3.1")




overall_pred_C.1.9 <- predictions(
  gam_C.1.9,
  newdata = datagrid(
    time     = time_grid_C.1.9,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9")


overall_pred_C.1.9.3 <- predictions(
  gam_C.1.9.3,
  newdata = datagrid(
    time     = time_grid_C.1.9.3,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9.3
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9.3")





# 3. Combine into a single data frame
overall_predictions <- bind_rows(
  
  overall_pred_D.3.1,
  
  overall_pred_C.1.9,
  overall_pred_C.1.9.3
) 


# 4. Plot the results
results <- ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_point(data = dt_plot %>% filter(subclade %in% c("D.3.1","C.1.9.3","C.1.9"))
             , aes(x = time,  y = indicies, colour = "grey"),alpha = 0.1)+
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(subclade), labeller = label_both) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()

# Display the plot
print(results)


ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()



ggplot(dt_plot %>% filter(subclade %in% c("D.3","D.3.1","C.1", "C.1.9","C.1.9.3")), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))


######ALL SWISS GISAID MODELLING#### 

file = "results_gisaid_2526_Swiss/index_time_parameters_ww_time_swiss.tsv"


dt = file %>% map_dfr(~ fread(.x)) %>% mutate(date = date_decimal(time))

dt_plot = dt %>% 
filter(time_window == 0.019 & timescale == 0.4) %>%
filter(date <= as.Date("2025-03-15")) %>%
filter(is.na(subclade) == FALSE)



#no amplicon term now


# gam_D.3.1 <- NULL

gam_D.3.1 <- bam(indicies ~  
                   
                   s(time, bs = 'cr', k=9) ,
                 
                 
                 data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)





gam_C.1.9 <- bam(indicies ~  
                   
                   s(time, bs = 'cs', k=9),
                 data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)




gam_C.1.9.3 <- bam(indicies ~  
                     
                     
                     s(time, bs = 'cs', k=9),
                   
                   data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                   family = gaussian(link = identity),
                   method = "fREML"
)







gam <- gam_D.3.1

AIC(gam)



summary(gam)



gam.check(gam)





plot_predictions(gam,
                 condition = list("time"),
                 points = 0.5,
                 newdata = dt_plot %>% filter(subclade == "C.1.9"))+
  scale_fill_manual(values = custom_palette) +
  facet_grid(
    # cols = vars(amplicon),
    rows = vars(subclade),
    labeller = label_both) 





time_grid_D.3.1 <- dt_plot %>% filter(subclade == "D.3.1") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9 <- dt_plot %>% filter(subclade == "C.1.9") %>% with(seq(min(time), max(time), length.out = 200))
time_grid_C.1.9.3 <- dt_plot %>% filter(subclade == "C.1.9.3") %>% with(seq(min(time), max(time), length.out = 200))




overall_pred_D.3.1 <- predictions(
  gam_D.3.1,
  newdata = datagrid(
    time     = time_grid_D.3.1,
    amplicon = levels(dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_D.3.1
  ),
  by = "time"
) %>% mutate(subclade = "D.3.1")




overall_pred_C.1.9 <- predictions(
  gam_C.1.9,
  newdata = datagrid(
    time     = time_grid_C.1.9,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9")


overall_pred_C.1.9.3 <- predictions(
  gam_C.1.9.3,
  newdata = datagrid(
    time     = time_grid_C.1.9.3,
    amplicon = levels(dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)) %>% pull(amplicon)),
    model    = gam_C.1.9.3
  ),
  by = "time"
) %>% mutate(subclade = "C.1.9.3")





# 3. Combine into a single data frame
overall_predictions <- bind_rows(
  
  overall_pred_D.3.1,
  
  overall_pred_C.1.9,
  overall_pred_C.1.9.3
) 


# 4. Plot the results
results <- ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_point(data = dt_plot %>% filter(subclade %in% c("D.3.1","D.3","C.1.9.3","C.1.9","C.1"))
             , aes(x = time,  y = indicies, colour = "grey"),alpha = 0.1)+
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(subclade), labeller = label_both) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()

# Display the plot
print(results)


ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()



ggplot(dt_plot %>% filter(subclade %in% c("D.3","D.3.1","C.1", "C.1.9","C.1.9.3")), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

