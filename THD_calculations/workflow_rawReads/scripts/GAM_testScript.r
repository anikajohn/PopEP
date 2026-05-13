# library(tidyr)
library(ggplot2)
# library(stringr)
library(data.table)
library(pals)
library(mgcv)
library(marginaleffects)
library(tidyverse)
library(car)

setwd("/Users/anjohn/Desktop/euler/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/")

#H1N1 data looks best, so doing it with thsi first 

###DATA PREPERATION ### 

amp2 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon2/combined_index_time_parameters.tsv")
amp3 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon3/combined_index_time_parameters.tsv")
amp4 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon4/2000n_downsampled/THD_calculation/index_time_parameters.tsv")
amp5 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon5/combined_index_time_parameters.tsv")
amp6 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon6/combined_index_time_parameters.tsv")

palette_colors <- as.vector(pals::cols25(14))
subclades <- c("C.1.9", "C.1.9.1", "C.1.9.2", "C.1.9.3",
               "C.1.9.4","D","D.1","D.3.1","D.5",
               "C.1","C.1.1","C.1.8","D.2","D.3")
custom_palette <- setNames(palette_colors, subclades)


dt_all = rbind(amp2,amp3,amp4,amp5,amp6)

dt = dt_all[time_window == 0.038 & timescale == 0.35]
dt[, amp_label := str_extract(seqName, "amplicon\\d+")]
dt[, location_code := sub(".*_([^_]+)_\\d{4}_\\d{2}_\\d{2}/.*", "\\1", seqName)]


dt_amp3 = amp3[time_window == 0.038 & timescale == 0.35]
dt_amp4 = amp4[time_window == 0.038 & timescale == 0.35]

dt_amp3[, location_code := sub(".*_([^_]+)_\\d{4}_\\d{2}_\\d{2}/.*", "\\1", seqName)]
dt_amp3[, sample := sub(".*_([A-Z][0-9]_.*)/.*", "\\1", seqName)]
dt_amp3[, batch := sub(".*/(.*)_.*", "\\1", seqName)]


dt_amp4[, location_code := sub(".*_([^_]+)_\\d{4}_\\d{2}_\\d{2}/.*", "\\1", seqName)]

cov_amp3 = fread("../results_H1_24_25_updated_ref_dataset/over_all_locations/amplicon3/amplicon_coverage/combined_coverage_processed.tsv")

dt_test = merge(dt_amp3,cov_amp3)

####UNDERSTANDING DATA####

ggplot(dt_test, aes(x = indicies, fill = subclade))+
  geom_histogram() +
  scale_fill_manual(values = custom_palette) + 
  facet_wrap(vars(location_code))

ggplot(dt_test, aes(x = indicies, fill = subclade))+
  geom_histogram() +
  scale_fill_manual(values = custom_palette)

ggplot(dt_test, aes(x = log(indicies), fill = subclade))+
  geom_histogram() +
  scale_fill_manual(values = custom_palette)



ggplot(dt_test, 
       aes(x = time, y = indicies, colour = location_code)) +
  geom_point()+
  geom_hline(yintercept=0.39,linetype='dashed', alpha = 0.5)








ggplot(dt_test, aes(x = time, y = indicies, colour = subclade))+
  geom_point() +
  geom_smooth(method = "gam")+
  scale_color_manual(values = custom_palette)





ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x = numreads, y = log(indicies), colour = subclade)) +
  geom_point()+
  scale_color_manual(values = custom_palette)

ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x= log(indicies), fill = subclade)) +
  geom_histogram()+
  scale_fill_manual(values = custom_palette)

ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x= indicies, fill = subclade)) +
  geom_histogram()+
  scale_fill_manual(values = custom_palette)



ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")] %>% filter(log(indicies) > -5), 
       aes(x= log(indicies), fill = subclade)) +
  geom_histogram()+
  scale_fill_manual(values = custom_palette)


dt_test$logit <- log(dt_test$indicies / (1 - dt_test$indicies))


ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x= logit, fill = subclade)) +
  geom_histogram()+
  scale_fill_manual(values = custom_palette)

###Checking zero-inflation#####

ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x = indicies)) +
  geom_histogram(fill = "steelblue", binwidth = 0.005) +
  geom_vline(xintercept=0.39,linetype='dashed', alpha = 0.5)

ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x = logit(indicies))) +
  geom_histogram(fill = "steelblue", binwidth = 0.025)+
  geom_vline(xintercept=logit(0.39),linetype='dashed', alpha = 0.5)

#Have a lot of highly abundant indicies >= 0.4


ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x = numreads, y = indicies, colour = subclade)) +
  geom_point()+
  scale_color_manual(values = custom_palette)

#highly abundant indicies (>= 0.4) mainly from C.1.9.3


ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x = time, y = indicies, colour = numreads)) +
  geom_point()+
  geom_hline(yintercept=0.39,linetype='dashed', alpha = 0.5)





ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x = time, y = logit(indicies), colour = numreads)) +
  geom_point()+
  geom_hline(yintercept=logit(0.39),linetype='dashed', alpha = 0.5)


ggplot(dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
       aes(x = time, y = logit(indicies), colour = subclade)) +
  geom_point()+
  scale_color_manual(values = custom_palette) +
  geom_hline(yintercept=logit(0.39),linetype='dashed', alpha = 0.5)












# library(performance)
# library(glmmTMB)
# 
# 
# model <- glmmTMB(indicies ~ time, 
#                  data = dt_test[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], 
#                  family = ziBeta_family(),  
#                  ziformula = ~1)
# 
# # Run the check
# check_zeroinflation(model)
# 




########MODELLING######




#FIT WITH GAUSSIAN (similar to run_gam_model())
gam_amp3 <- bam(indicies ~  s(time, bs = 'cr') +
                  #by = seperate smooth function for each factor
                  s(time, by = as.factor(subclade), bs="cr", k = 15) +  as.factor(subclade) ,
                  # denoising data
                  data = dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")],
                  family = gaussian(link = log),
                  method = "fREML"
                )

gam_amp3

# gam_amp4 <- bam(indicies ~ 
#                   s(time, bs = 'cr', k = 15) +
#                   # s(location_code, bs = 're') +
#                   #by = seperate smooth function for each factor
#                   s(time, by = as.factor(subclade), bs="cr", k = 15) +  as.factor(subclade) ,
#                 # data = dt_amp3,
#                 data = dt_amp4[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1", "C.1","D.3")],
#                 family = betar(link = logit),
#                 method = "fREML"
# )
# 
# gam_amp4
# ggplot(dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x=car::logit(indicies))) + 
#   geom_histogram()+
#   facet_wrap(vars(subclade))
# 
# 
# ggplot(dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x=log(indicies))) + 
#   geom_histogram() +
#   facet_wrap(vars(subclade))
# 
# 
# 
# ggplot(dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x=log(indicies)))+
#   geom_histogram() +
#   facet_wrap(vars(subclade))
# 
# ggplot(dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x=indicies))+
#   geom_histogram() +
#   facet_wrap(vars(subclade))
# 

AIC(gam_amp3)


 
summary(gam_amp3)


 
gam.check(gam_amp3)




plot_predictions(gam_amp3,
                 condition = list("time","subclade"),
                 points = 0.5,
                 newdata = dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")])+
  scale_color_manual(values = custom_palette) 








#ENCODE subclades same as Phylowave

#Compare data distribution from Phylowave pattern

dt_phylowaye = fread("/Users/anjohn/Downloads/List_all_isolates_H3N2_global_20230414.csv")


ggplot(dt_phylowaye %>% filter(Collection_time > 2010), aes(x=Collection_time))+
  geom_histogram()

ggplot(dt_amp3[subclade %in% c("C.1.9", "C.1.9.3", "D.3.1")], aes(x=time))+
  geom_histogram()





clade_levels <- levels(as.factor(dt_amp3$subclade))
clade_data <- sapply(clade_levels, function(cl) as.integer(as.factor(dt_amp3$subclade)== cl))
colnames(clade_data) <- clade_levels
clade_data <- data.frame(index = dt_amp3$indicies,time = dt_amp3$time, clade_data)



t = bam(log(index) ~ s(time, bs = 'cr') +
    s(D.3.1, bs = 'cr', k = 2, by = D.3.1)+ 
    s(C.1.9, bs = 'cr', by = C.1.9)+
    s(C.1.9.3, bs = 'cr', by = C.1.9.3),
    family = gaussian(link = identity),
    data = clade_data,
    method = 'GCV.Cp')
      


