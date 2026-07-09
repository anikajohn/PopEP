
library(data.table)
library(mgcv)
library(tidyverse)
library(marginaleffects)
library(pals)




setwd("/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_debug_haplo_H1_24_25")

files = list.files("over_all_locations/", 
           pattern = "combined_index_time_parameters_processed.tsv",
           recursive = TRUE,
           full.names = TRUE)


files_per_loc = list.files(".", 
                   pattern = "^index_time_parameters_processed.tsv",
                   recursive = TRUE,
                   full.names = TRUE)


dts = files %>% map_dfr(~ fread(.x))

dts_per_loc = files_per_loc %>% map_dfr(~ fread(.x))

# dts_abundance = files_abundance %>% map_dfr(~ fread(.x))
# dts_abundance_plot = dts_abundance %>% 
#   group_by(amplicon,location,subclade,date) %>% 
#   summarise(count_sum = sum(count)) %>% 
#   filter(date <= as.Date("2025-03-15")) %>%
#   filter(is.na(subclade) == FALSE) %>% 
#   
  

# dt_byAmp = dts %>% 
#   group_by(subclade,`legacy-clade`,timescale,time_window,amplicon,date) %>% 
#   summarise(index_average = sum(indicies*count)/sum(count))
# 
# 
# dt_acrossAmp = dts %>% 
#   filter(time_window == 0.019 & timescale == 0.4) %>% 
#   group_by(subclade,`legacy-clade`,timescale,time_window,date) %>% 
#   summarise(index_average = sum(indicies*count)/sum(count))








palette_colors <- as.vector(pals::cols25(14))
subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
               "C.1.9.4","D","D.1","D.3.1","D.5",
               "C.1","C.1.1","C.1.8","D.2","D.3")
custom_palette <- setNames(palette_colors, subclades)


#######ACROSS AMPLICON ##########

dt_plot = dts %>% 
  filter(time_window == 0.019 & timescale == 0.4) %>% 
  filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE) %>% 
  mutate(naive_weights = case_when(amplicon == "amplicon3" ~ count * 2,
                                   .default = count))
dt_plot_per_loc = dts_per_loc %>% 
  filter(time_window == 0.019 & timescale == 0.4) %>% 
  filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE)


dt_day <- dt_plot %>%
  group_by(date, subclade, location, amplicon) %>%
  summarise(count_sum = sum(count, na.rm = TRUE), .groups = "drop")

dt_day %>%
  filter(subclade %in% c("C.1.9.3", "C.1.9", "D.3.1")) %>%
  count(date,location, amplicon, subclade) %>% view(
    
  )


by_day_amp_loc = ggplot(dt_day %>% filter(subclade %in% c("C.1.9.3","C.1.9","D.3.1")), 
                                  aes(color=subclade, x=date, y=count_sum, fill = subclade)) + 
  # geom_col(position="stack")+
  geom_bar(stat = "identity")+
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette)+
  facet_grid(
    cols = vars(amplicon),
    rows = vars(location)
             ) +
  ylab("sequence count")+
  scale_y_log10() +
  theme(
        axis.text.x = element_text(angle = 45, hjust = 1))

by_day_amp = ggplot(dt_day %>% filter(subclade %in% c("C.1.9.3","C.1.9","D.3.1")), 
                                  aes(color=subclade, x=date, y=count_sum, fill = subclade)) + 
  # geom_col(position="stack")+
  geom_bar(stat = "identity")+
  scale_color_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette)+
  facet_grid(
    cols = vars(amplicon)
    # rows = vars(location)
             ) +
  ylab("sequence count")+
  scale_y_log10() +
  theme(
        axis.text.x = element_text(angle = 45, hjust = 1))




ggsave("by_day_amp.png", plot = by_day_amp)
ggsave("by_day_amp_loc.png", plot = by_day_amp_loc)


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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
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
  geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  ggtitle("count weights")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))


########### BY AMPLICON & LOCATION ######



ggplot(dt_plot_per_loc %>% filter(subclade %in% c("C.1.9.3","C.1.9","D.3.1")), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  # geom_smooth(aes(x = date, y = indicies, colour = subclade, weight = count),
  #             method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(
    # cols = vars(amplicon),
     rows = vars(location),
     labeller = label_both
    ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(0, 1))



#####GAM MODELLING #####
amp3_weight = 0.2

dt_plot = dt_plot %>%
  mutate(amp_weighted = case_when(
  amplicon == "amplicon3" ~ count * amp3_weight,
  .default = count *(1-amp3_weight)/4))

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

# # 2. Get overall predictions for each model
# overall_pred_D.3 <- predictions(gam_D.3, newdata = datagrid(time = time_grid_D.3)) %>% mutate(subclade = "D.3")
# overall_pred_D.3.1 <- predictions(gam_D.3.1, newdata = datagrid(time = time_grid_D.3.1)) %>% mutate(subclade = "D.3.1")
# overall_pred_C.1 <- predictions(gam_C.1, newdata = datagrid(time = time_grid_C.1)) %>% mutate(subclade = "C.1")
# overall_pred_C.1.9 <- predictions(gam_C.1.9, newdata = datagrid(time = time_grid_C.1.9)) %>% mutate(subclade = "C.1.9")
# overall_pred_C.1.9.3 <- predictions(gam_C.1.9.3, newdata = datagrid(time = time_grid_C.1.9.3)) %>% mutate(subclade = "C.1.9.3")


#2.1 weighted average over balanced prediction
# User-defined amplicon weights (edit these values as needed)
# amplicon_weights <- c(
#   "amplicon2" = 0.2,
#   "amplicon3" = 1,
#   "amplicon4" = 1,
#   "amplicon5" = 1,
#   "amplicon6" = 1
# )
# 
# # Helper using user-supplied weights
# # Helper using user-supplied weights
# make_weighted_grid <- function(model, subclade_name, time_grid) {
#   subdata <- dt_plot %>%
#     filter(subclade == subclade_name) %>%
#     mutate(amplicon = as.factor(amplicon))
#   
#   datagrid(
#     time     = time_grid,
#     amplicon = levels(subdata$amplicon),
#     model    = model
#   ) %>%
#     mutate(wts = amplicon_weights[as.character(amplicon)])
# }
# 
# # 2.1 Weighted average over amplicons at each time point
# overall_pred_D.3 <- predictions(
#   gam_D.3,
#   newdata = make_weighted_grid(gam_D.3, "D.3", time_grid_D.3),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "D.3")
# 
# overall_pred_D.3.1 <- predictions(
#   gam_D.3.1,
#   newdata = make_weighted_grid(gam_D.3.1, "D.3.1", time_grid_D.3.1),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "D.3.1")
# 
# overall_pred_C.1 <- predictions(
#   gam_C.1,
#   newdata = make_weighted_grid(gam_C.1, "C.1", time_grid_C.1),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1")
# 
# overall_pred_C.1.9 <- predictions(
#   gam_C.1.9,
#   newdata = make_weighted_grid(gam_C.1.9, "C.1.9", time_grid_C.1.9),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1.9")
# 
# overall_pred_C.1.9.3 <- predictions(
#   gam_C.1.9.3,
#   newdata = make_weighted_grid(gam_C.1.9.3, "C.1.9.3", time_grid_C.1.9.3),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1.9.3")
# 


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







# ####PRECISION WEIGHTS #######
# 
# make_precision_weighted_grid <- function(model, subclade_name, time_grid) {
#   subdata <- dt_plot %>%
#     filter(subclade == subclade_name) %>%
#     mutate(amplicon = as.factor(amplicon))
#   
#   grid <- datagrid(
#     time     = time_grid,
#     amplicon = levels(subdata$amplicon),
#     model    = model
#   )
#   
#   # Step 1: get cell-level predictions to obtain per-cell SE
#   cell_preds <- predictions(model, newdata = grid)
#   
#   # Step 2: combined weight = user weight × precision (1 / SE²)
#   cell_preds %>%
#     mutate(
#       user_wt = amplicon_weights[as.character(amplicon)],
#       prec_wt = 1 / std.error^2,
#       wts     = user_wt * prec_wt
#     )
# }
# 
# 
# 
# overall_pred_D.3 <- predictions(
#   gam_D.3,
#   newdata = make_precision_weighted_grid(gam_D.3, "D.3", time_grid_D.3),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "D.3")
# 
# overall_pred_D.3.1 <- predictions(
#   gam_D.3.1,
#   newdata = make_precision_weighted_grid(gam_D.3.1, "D.3.1", time_grid_D.3.1),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "D.3.1")
# 
# overall_pred_C.1 <- predictions(
#   gam_C.1,
#   newdata = make_precision_weighted_grid(gam_C.1, "C.1", time_grid_C.1),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1")
# 
# overall_pred_C.1.9 <- predictions(
#   gam_C.1.9,
#   newdata = make_precision_weighted_grid(gam_C.1.9, "C.1.9", time_grid_C.1.9),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1.9")
# 
# overall_pred_C.1.9.3 <- predictions(
#   gam_C.1.9.3,
#   newdata = make_precision_weighted_grid(gam_C.1.9.3, "C.1.9.3", time_grid_C.1.9.3),
#   by = "time", wts = "wts"
# ) %>% mutate(subclade = "C.1.9.3")


# 3. Combine into a single data frame
overall_predictions <- bind_rows(
  # overall_pred_D.3,
  overall_pred_D.3.1,
  # overall_pred_C.1,
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
  # facet_grid(cols = vars(subclade), labeller = label_both) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()

# Display the plot
print(results)


ggplot(overall_predictions, aes(x = time, y = estimate, colour = subclade)) +
  geom_point(data = dt_plot %>% filter(subclade %in% c("D.3.1","C.1.9.3","C.1.9"))
             , aes(x = time,  y = indicies, colour = "lightgrey"), alpha = 0.1)+
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
  coord_cartesian(ylim = c(-0.5, 1))+
  theme_bw()









 

dt_abundance <- dt_plot %>%
  filter(subclade %in% c("C.1","C.1.9","C.1.9.3", "D.3", "D.3.1")) %>% 
  mutate(week_date = floor_date(as.Date(date), unit = "week")) %>% 
  group_by(week_date, subclade, amplicon) %>%
  summarise(count_sum = sum(count, na.rm = TRUE), 
            count_sum_amp = sum(amp_weighted, na.rm = TRUE),
            .groups = "drop")


abundance <- 
  ggplot(dt_abundance, aes(fill=subclade, x=week_date, y=count_sum_amp)) + 
     geom_col(position="stack")+
     scale_fill_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             rows = vars(subclade),
             labeller = label_both) +
     ylab("sequence count 
          amp weighted")+
     theme(text = element_text(size = 20),
                     axis.text.x = element_text(angle = 45, hjust = 1))

abundance



