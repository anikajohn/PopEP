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

window_thr = 0.019
scale_thr = 0.4

palette_colors <- as.vector(pals::cols25(14))
subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
               "C.1.9.4","D","D.1","D.3.1","D.5",
               "C.1","C.1.1","C.1.8","D.2","D.3")
custom_palette <- setNames(palette_colors, subclades)

target_subclades = c("C.1.9.3","C.1.9","D.3.1","D.3")


#######ACROSS AMPLICON ##########

dt_plot = dts %>% 
  filter(time_window == window_thr & timescale == scale_thr) %>% 
  filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE) %>% 
  mutate(week_date = floor_date(as.Date(date), unit = "week"))



seq_included = dt_plot %>%  pull(seqName) %>% unique()

########### OVERALL ABUNDANCE ######



ggplot(dt_plot,
         # filter(date >=  date_min & date <= date_max), 
       aes(x=week_date, fill = subclade))+
  geom_bar()+
  scale_fill_manual(values = custom_palette)+
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


########### OVERALL THD ######
ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("Per sequence fitness estimate")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")




ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE, formula = y ~ s(x, k = 3)) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")





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
                 s(time, bs = 'cs', k = 2) + 
                 s(time, amplicon, bs="fs",k = 2),
               
               # gamma = 5,
               
               
               # weights = amp_weighted,
               data = dt_plot %>% filter(subclade == "D.3") %>% mutate(amplicon = as.factor(amplicon)),
               family = betar(link = "logit"),
               method = "fREML"
               
)


gam_D.3.1 <- bam(indicies ~  
                # s(time, bs = 'cr')+
                # s(time, by = as.factor(amplicon), bs="cr", k = -1) +
                # as.factor(amplicon),

                  s(time, bs = 'cs', k = 2) + 
                  s(time, amplicon, bs="fs",k = 2),
                
                # gamma = 5,
                
                
                # weights = amp_weighted,
                data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                family = betar(link = "logit"),
                method = "fREML"
                
)




gam_C.1 <- bam(indicies ~  
                 s(time, bs = 'cs', k = 2) + 
                 s(time, amplicon, bs="fs",k = 2),
               
               # gamma = 5,
               
                 
                 # weights = amp_weighted,
                 data = dt_plot %>% filter(subclade == "C.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = betar(link = "logit"),
                 method = "fREML"
                 
)


gam_C.1.9 <- bam(indicies ~  
                   s(time, bs = 'cs', k = 2) + 
                   s(time, amplicon, bs="fs",k = 2),
                 
                 # gamma = 5,
                 
                
                # weights = amp_weighted,
                data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                family = betar(link = "logit"),
                method = "fREML"
                
                )




gam_C.1.9.3 <- bam(indicies ~  
                     s(time, bs = 'cs', k = 2) + 
                     s(time, amplicon, bs="fs",k = 2),
                   
                   # gamma = 5,
                   
                
                 # weights = amp_weighted,
                 data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                 family = betar(link = "logit"),
                 method = "fREML"
)







# gam <- gam_C.1.9
# 
# AIC(gam)
# 
# 
# 
# summary(gam)
# 
# 
# 
# gam.check(gam)
# 
# 
# 
# 
# 
# 
# plot_predictions(gam_C.1.9,
#                  condition = list("time", "amplicon"),
#                  points = 0.5,
#                  newdata = dt_plot %>% filter(subclade == "C.1.9"))+
#   coord_cartesian(ylim = c(-1, 1))+
#   scale_fill_manual(values = custom_palette) +
#   facet_grid(cols = vars(amplicon),
#     rows = vars(subclade),
#     labeller = label_both) 
# 
# 
# 
# plot_predictions(gam_C.1.9,
#                  condition = list("time"),
#                  points = 0.5,
#                  newdata = dt_plot %>% filter(subclade == "C.1.9"))+
#   # coord_cartesian(ylim = c(-2, 3))+
#   coord_cartesian(ylim = c(-1, 1))+
#   scale_fill_manual(values = custom_palette) +
#   facet_grid(
#             # cols = vars(amplicon),
#              rows = vars(subclade),
#              labeller = label_both) 









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
)  %>% 
  mutate(date = lubridate::date_decimal(time)) %>% 
  filter(subclade %in% target_subclades)


# 4. Plot the results
results <- ggplot(overall_predictions, aes(x = date, y = estimate, colour = subclade)) +
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


ggplot(overall_predictions, aes(x = date, y = estimate, colour = subclade)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade), alpha = 0.2) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(-0.5, 1))
  theme(legend.position = "bottom",
  text = element_text(size = 20))

#########################################
######FULL LENGTH SEQUENCES MODELLING#### 
# 1Europe
#########################################


file = "results_gisaid_2526/index_time_parameters_ww_time_fullLength.tsv"


dt = file %>% 
  map_dfr(~ fread(.x)) %>% mutate(date = date_decimal(time), amplicon = "full_length")



###Little excursion to check how much wrongly classified ####


dt_ground_truth = dt %>% filter(seqName %in% seq_included) %>% 
  filter(time_window == window_thr & timescale == scale_thr) %>% 
  # filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE) %>% 
  mutate(week_date = floor_date(as.Date(date), unit = "week"))


dt_ground_truth = dt_ground_truth[,c("seqName","subclade","legacy.clade","date","amplicon")]

setnames(dt_ground_truth, c("subclade", "legacy.clade","amplicon"), 
         c("subclade_gt", "legacy.clade_gt","amplicon_gt"))


dt_verification = merge(dt_plot,dt_ground_truth) %>% 
  filter(subclade == subclade_gt)
  

D.3.1_sequences_gt = dt_verification %>% 
  filter(subclade == "D.3.1") %>% 
  pull(seqName) %>% unique()

C.1.9.3_sequences_gt = dt_verification %>% 
  filter(subclade == "C.1.9.3") %>% 
  pull(seqName) %>% unique()

C.1.9_sequences_gt = dt_verification %>% 
  filter(subclade == "C.1.9") %>% 
  pull(seqName) %>% unique()

write



ggplot(dt_verification,
       # filter(date >=  date_min & date <= date_max), 
       aes(x=week_date, fill = subclade))+
  geom_bar()+
  scale_fill_manual(values = custom_palette)+
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


########### OVERALL THD ######
ggplot(dt_verification %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("Per sequence fitness estimate")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")




ggplot(dt_verification %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE, formula = y ~ s(x, k = 3)) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


ggplot(dt %>% filter(seqName %in% correctly_classified) %>% 
         filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE, formula = y ~ s(x, k = 3)) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")





########### D.3.1 ######

#BY AMPLICON

ggplot(dt_verification %>% filter(subclade == "D.3.1"), 
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

ggplot(dt_verification %>% filter(subclade == "D.3.1"), 
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

ggplot(dt_verification %>% filter(subclade == "C.1.9"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Overall
ggplot(dt_verification %>% filter(subclade == "C.1.9"), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
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











#####processing of full length clinical data continuous ######


dt_plot = dt %>% filter(seqName %in% seq_included) %>% 
  filter(time_window == window_thr & timescale == scale_thr) %>% 
  # filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE) %>% 
  mutate(week_date = floor_date(as.Date(date), unit = "week"))


D.3.1_sequences = dt_plot %>% filter(subclade == "D.3.1") %>% pull(seqName) %>% unique()
C.1.9.3_sequences = dt_plot %>% filter(subclade == "C.1.9.3") %>% pull(seqName) %>% unique()
C.1.9_sequences = dt_plot %>% filter(subclade == "C.1.9") %>% pull(seqName) %>% unique()

write(D.3.1_sequences, 
          "~/Desktop/THD_calculations_local/nextclade/results_gisaid_2526/D.3.1_sequence_names.tsv",
          sep = "\t")

write(C.1.9_sequences, 
      "~/Desktop/THD_calculations_local/nextclade/results_gisaid_2526/C.1.9_sequence_names.tsv",
      sep = "\t")

write(C.1.9.3_sequences, 
      "~/Desktop/THD_calculations_local/nextclade/results_gisaid_2526/C.1.9.3_sequence_names.tsv",
      sep = "\t")



########### OVERALL ABUNDANCE ######

ggplot(dt_plot,
       # filter(date >=  date_min & date <= date_max), 
       aes(x=week_date, fill = subclade))+
  geom_bar()+
  scale_fill_manual(values = custom_palette)+
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


########### OVERALL THD ######
ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("Per sequence fitness estimate")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")




ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE, formula = y ~ s(x, k = 2)) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


  




gam_D.3.1 <- bam(indicies ~  
                   
                   s(time, bs = 'cs', k = 2) + 
                   s(time, amplicon, bs="fs",k = 2),
                 
                
                 data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)





gam_C.1.9 <- bam(indicies ~  
                   
                   s(time, bs = 'cs', k = 2) + 
                   s(time, amplicon, bs="fs",k = 2),
                 
                 data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)




gam_C.1.9.3 <- bam(indicies ~  
                     s(time, bs = 'cs', k = 2) + 
                     s(time, amplicon, bs="fs",k = 2),
                   
                   
                   data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                   family = gaussian(link = identity),
                   method = "fREML"
)















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
) %>%  mutate(date = lubridate::date_decimal(time)) %>% 
  filter(subclade %in% target_subclades)

# 4. Plot the results
results <- ggplot(overall_predictions, aes(x = date, y = estimate, colour = subclade)) +
  geom_point(data = dt_plot %>% filter(subclade %in% c("D.3.1","C.1.9.3","C.1.9"))
             , aes(x = date,  y = indicies, colour = "grey"),alpha = 0.1)+
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
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Display the plot
print(results)


ggplot(overall_predictions, aes(x = date, y = estimate, colour = subclade)) +
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
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")



ggplot(overall_predictions, aes(x = date, y = estimate, colour = subclade_simple)) +
  geom_point(data = dt_plot_ww %>% filter(subclade_simple %in% target_subclades)
             , aes(x = date,  y = indicies), alpha = 0.1,colour = "lightgrey")+
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = subclade_simple), alpha = 0.2) +
  facet_grid(cols = vars(subclade_simple), labeller = label_both) +
  scale_colour_manual(values = custom_palette) +
  scale_fill_manual(values = custom_palette) +
  labs(
    # title = "Overall GAM Prediction (Averaged over Amplicons)",
    x = "Time",
    y = "Predicted Fitness Index"
  ) +
  # coord_cartesian(ylim = c(0, 1))+
  theme_bw()



######ALL SWISS GISAID MODELLING#### 

file = "results_gisaid_2526_Swiss/index_time_parameters_ww_time_swiss.tsv"


dt = file %>% 
  map_dfr(~ fread(.x)) %>% mutate(date = date_decimal(time), amplicon = "full_length")

dt_plot = dt %>% 
  filter(time_window == window_thr & timescale == scale_thr) %>% 
  filter(date <= as.Date("2025-03-15")) %>%
  filter(is.na(subclade) == FALSE) %>% 
  mutate(week_date = floor_date(as.Date(date), unit = "week"))



########### OVERALL ABUNDANCE ######

ggplot(dt_plot,
       # filter(date >=  date_min & date <= date_max), 
       aes(x=week_date, fill = subclade))+
  geom_bar()+
  scale_fill_manual(values = custom_palette)+
  facet_grid(cols = vars(amplicon))+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


########### OVERALL THD ######
ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  ggtitle("Per sequence fitness estimate")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")




ggplot(dt_plot %>% filter(subclade %in% target_subclades), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE, formula = y ~ s(x, k = 3)) +
  scale_colour_manual(values = custom_palette) +
  facet_grid(cols = vars(amplicon),
             labeller = label_both) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


# gam_D.3.1 <- NULL

gam_D.3.1 <- bam(indicies ~  
                   
                   s(time, bs = 'cs', k = 2) + 
                   s(time, amplicon, bs="fs",k = 2),
                 
                 
                 data = dt_plot %>% filter(subclade == "D.3.1") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)





gam_C.1.9 <- bam(indicies ~  
                   
                   s(time, bs = 'cs', k = 2) + 
                   s(time, amplicon, bs="fs",k = 2),
                 
                 data = dt_plot %>% filter(subclade == "C.1.9") %>% mutate(amplicon = as.factor(amplicon)),
                 family = gaussian(link = identity),
                 method = "fREML"
                 
)




gam_C.1.9.3 <- bam(indicies ~  
                     s(time, bs = 'cs', k = 2) + 
                     s(time, amplicon, bs="fs",k = 2),
                   
                   
                   data = dt_plot %>% filter(subclade == "C.1.9.3") %>% mutate(amplicon = as.factor(amplicon)),
                   family = gaussian(link = identity),
                   method = "fREML"
)














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
) %>% mutate(date = lubridate::date_decimal(time)) %>% 
  filter(subclade %in% target_subclades)



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
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")




ggplot(dt_plot %>% filter(subclade %in% c("D.3","D.3.1","C.1", "C.1.9","C.1.9.3")), 
       aes(x = date, y = indicies, colour = subclade)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(x = date, y = indicies, colour = subclade),
              method = "gam", se = TRUE) +
  scale_colour_manual(values = custom_palette) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

