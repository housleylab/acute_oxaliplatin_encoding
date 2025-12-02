###### STA paper

########################### ReadME begin ########################### 
# Version of Nov 6, 2025
# Stephen N. Housley 
# nickhousley@gatech.edu

# This work is licensed under the licenses 
# Paper: Creative Commons Attribution 3.0 Unported License 
# Code: GPL-3 License 
# Depends: R (>= 3.5.0)
# Version: 0.1
# Description: code to run analytics and graphic functions associated with:
#         STA paper
#         
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 

## version Hx
# v0.1- original Nov 6, 2025

# paths must be changed to accommodate end user file structure (e.g. line 46)
# run on MacOS 14.6.1 

########################### ReadME end ########################### 

##### begin TEMPLATE##### 
########################### Figure ZZZZ ###########################
########################### description
########################### load dependencies
########################### custom functions
########################### load data
########################### data wrangling
########################### quick visualization
########################### analyses/modeling
########################### saving data
########################### saving figures
########################### clean up
##### end TEMPLATE #####


########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")

########################### load general dependencies ########################### 
source("code/load_gen_dependencies.R")



########################### Figure 4i ###########################
########################### description

########################### Load Data
set.seed(0123456789)
############ Generate synthetic dataset for conduction velocities based on literature

## A alpha, 30-55 m/s.  https://physoc.onlinelibrary.wiley.com/doi/abs/10.1113/jphysiol.1985.sp015574
## myelinated mechanoreceptors (mode 25–30 m/s, https://www.sciencedirect.com/science/article/abs/pii/0006899382907685
## mean 32.6–46.6 m/s for various classes of mechanoreceptors, https://journals.physiology.org/doi/abs/10.1152/jn.1993.69.5.1684 
numSim = 1000
synD_cv_muscle <- rnorm(numSim, mean = 55, sd = 2)
synD_cv_cutan <- rnorm(numSim, mean = 45, sd = 2)



############## Scenario 1 if it comes from intraxonal site 
### intra-cellular record
synD_dist_intra_axonal_site <- rnorm(numSim, mean = 0, sd = 1)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_muscle

### extra-cellular record
synD_dist_pn_recording_site_m <- rnorm(numSim, mean = 85, sd = 1)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_muscle

### STA result 
intraxonal_cd <- muscle_cd + intraAxon_cd
hist(intraxonal_cd)


############## Scenario 2 if it comes from receptor site 
### intra-cellular record
synD_dist_intra_axonal_site <- rnorm(numSim, mean = 0, sd = 1)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_muscle

### extra-cellular record
synD_dist_pn_recording_site_m <- rnorm(numSim, mean = -85, sd = 1)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_muscle

### STA result 
receptorsite_musc_cd <- muscle_cd + intraAxon_cd
hist(receptorsite_musc_cd)


############## Scenario 3 if it comes from middle of peripheral nerve

### intra-cellular record
synD_dist_intra_axonal_site <- rnorm(numSim, mean = 42.5, sd = 1)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_muscle
hist(intraAxon_cd)

### extra-cellular record
synD_dist_pn_recording_site_m <- rnorm(numSim, mean = -42.5, sd = 1)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_muscle

### STA result 
peripheralsite_mid_cd <- muscle_cd + intraAxon_cd
hist(peripheralsite_mid_cd)


############## Scenario 4 if it comes from DRG

### intra-cellular record
synD_dist_intra_axonal_site <- rnorm(numSim, mean = -25, sd = 1)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_muscle
hist(intraAxon_cd)

### extra-cellular record
synD_dist_pn_recording_site_m <- rnorm(numSim, mean = 60, sd = 1)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_muscle
hist(muscle_cd)

### STA result 
drg_cd <- muscle_cd + intraAxon_cd
hist(drg_cd)


############## Scenario 5 if it comes from somewhere non-specific in peripheral nerve

### intra-cellular record
synD_dist_intra_axonal_site <- runif(numSim, min = 25, max = 85)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_muscle
hist(intraAxon_cd)

### extra-cellular record
synD_dist_pn_recording_site_m <- runif(numSim, min = -85, max = -25)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_muscle
hist(muscle_cd)

### STA result 
peripheralsite_random_cd <- muscle_cd + intraAxon_cd
hist(peripheralsite_random_cd)


############## Scenario 6 if it comes from cutaneous receptor site 
### intra-cellular record
synD_dist_intra_axonal_site <- rnorm(numSim, mean = 0, sd = 1)
intraAxon_cd = synD_dist_intra_axonal_site/synD_cv_cutan

### extra-cellular record
synD_dist_pn_recording_site_m <- rnorm(numSim, mean = -105, sd = 2)
muscle_cd = synD_dist_pn_recording_site_m/synD_cv_cutan

### STA result 
receptorsite_cut_cd <- muscle_cd + intraAxon_cd
hist(receptorsite_cut_cd)





########################### Data Wrangling: bind and convert to df
synD_df<- as.data.frame(cbind(receptorsite_musc_cd,
                              receptorsite_cut_cd,
                              peripheralsite_random_cd,
                              peripheralsite_mid_cd,
                              drg_cd,
                              intraxonal_cd
))

## wide to long format for graphing purposes
synD_df_long <- gather(synD_df, location, cd_ms, receptorsite_musc_cd:intraxonal_cd, factor_key=TRUE)
mu <- synD_df_long %>% 
  group_by(location) %>%
  summarise(loc.mean = mean(cd_ms))

########################### get observed data
obs_mus_df <- read_excel("data/acute ox muscle affs_A1.xls")
obs_cut_df <- read_excel("data/acute ox cutaneous affs_A1.xls")

########################### Data Wrangling
obs_mus_df_cd<- obs_mus_df %>% 
  filter(`spont? 0=control,1=no, 2=yes` == 0) %>%
  # filter(`aff type 1=ia, 2=II, 3=Ib` == 1) %>%
  dplyr::select(`conduction delay`) %>% 
  na.omit() %>%
  rename(cd_ms = `conduction delay`) %>%
  mutate(cd_ms = cd_ms*-1) %>%
  mutate(location=c('muscle_cd_obs'),
         .before=cd_ms)

obs_cut_df_cd<- obs_cut_df %>% 
  filter(`spont 0=control, 1=non spont, 2=spont` == 0) %>%
  # filter(`aff type 1=ra, 2=sa` == 2) %>%
  dplyr::select(`conduction delay`) %>% 
  na.omit() %>%
  rename(cd_ms = `conduction delay`) %>%
  mutate(cd_ms = cd_ms*-1) %>%
  mutate(location=c('cutan_cd_obs'),
         .before=cd_ms)

############ bind observed to synthetic data
synD_df_long<-rbind(synD_df_long, obs_mus_df_cd, obs_cut_df_cd)

mu <- synD_df_long %>% 
  group_by(location) %>%
  summarise(loc.mean = mean(cd_ms),
            loc.med = median(cd_ms)
  )

############ numerical summaries
synD_df_long %>% 
  group_by(location) %>% 
  summarise(mean = mean(cd_ms), 
            sd = sd(cd_ms),
            n = n()
  )

########################### quick visualization
############ simulated densities and observed 1D density plots
fig_4i<-synD_df_long %>%
  mutate(location = factor(location, levels = c('receptorsite_musc_cd',
                                                'receptorsite_cut_cd',
                                                'peripheralsite_random_cd',
                                                'peripheralsite_mid_cd',
                                                'drg_cd',
                                                'intraxonal_cd',
                                                'cutan_cd_obs',
                                                'muscle_cd_obs'))) %>%
  ggplot(aes(cd_ms, colour = location, fill=location))+
  geom_density(alpha=0.7) +
  # geom_vline(data=mu, aes(xintercept=loc.mean, color=location),
  #            linetype="dashed")+
  scale_fill_manual(values=c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#006D2C","#d95f02")) +
  scale_color_manual(values=c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#006D2C","#d95f02")) +
  theme_classic()+
  theme(legend.position = c(0.9, 0.8))+
  xlim(-4,4)

########################### saving data
########################### saving figures
ggsave(fig_4i, file = "fig_4i.pdf", width = 15, height = 15, units = "cm", path = "figures/")

########################### Clean up









########################### Figure 7a ###########################
########################### description
########################### load dependencies
########################### custom functions
sourceDirectory("code/functions/")
########################### load data
setwd("data/muscle_spike_times/")

########################### data wrangling
start.time <- Sys.time()

filelist = list.files(pattern = ".*.txt")
data = lapply(filelist, read.table, header = TRUE, sep = "\t") ## read in all files

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
gc()
# names(data) <- gsub("\\.txt$", "", fileNames)
names(data) <- stringr::str_replace(filelist, pattern = ".txt", replacement = "")   ## rename the list elements to reflect the Neuron names

lapply(seq_along(1:length(data)),function(x) names(data[[x]])<<-names(vars_rename(names(data[[x]]), ## rename the var list
                                                                                  recTime = contains("Time"), 
                                                                                  Spikes = starts_with("X1.")
                                                                                  
))) 

# t(apply(simplify2array(data), 1:2, mean))
# test<-lapply(data, function (x) x$Spikes)

df_tidy<-bind_rows(data, .id = "source_df")


df_tidy<-df_tidy %>%
  dplyr::mutate(
    sf_category = case_when(
      str_detect(source_df, "_SF-") ~ "sf_neg",
      str_detect(source_df, "_SF+") ~ "sf_pos",
      str_detect(source_df, "cont") ~ "control",
      TRUE ~ "noneDetected" # Default case if neither is found
    )
  )


df_tidy_mean <- df_tidy %>%
  filter(!is.na(recTime)) %>%
  group_by(recTime,sf_category) %>%
  summarise(n = n(),
            mean = mean(Spikes),
            median = median(Spikes),
            sd = sd(Spikes)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)
########################### quick visualization


muscleFig<-ggplot(df_tidy_mean, aes(x=recTime, y=mean, color = sf_category)) +
  geom_line(aes(x=recTime, y=mean,color = sf_category)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, fill=sf_category),color="grey70",alpha=0.4)+
  scale_fill_manual(values = c("#3A53A4","#A3A3CE", "#A4A4A5"))+
  scale_color_manual(values = c("#3A53A4","#A3A3CE", "#A4A4A5"))+
  theme_classic()  +
  facet_wrap(~ sf_category, nrow = 3)
########################### analyses/modeling
########################### saving data
########################### saving figures
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")
ggsave(muscleFig, file = "muscle_population_code_Fig_split.pdf", width = 18, height = 18, units = "cm", path = "figures")
########################### clean up
########################### Figure 7b ###########################
########################### description
########################### load dependencies
########################### custom functions
sourceDirectory("code/functions/")
########################### load data
setwd("data/cutaneous_spike_times/")

########################### data wrangling
start.time <- Sys.time()

filelist = list.files(pattern = ".*.txt")
data = lapply(filelist, read.table, header = TRUE, sep = "\t") ## read in all files

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
gc()
# names(data) <- gsub("\\.txt$", "", fileNames)
names(data) <- stringr::str_replace(filelist, pattern = ".txt", replacement = "")   ## rename the list elements to reflect the Neuron names

lapply(seq_along(1:length(data)),function(x) names(data[[x]])<<-names(vars_rename(names(data[[x]]), ## rename the var list
                                                                                  recTime = contains("Time"), 
                                                                                  Spikes = starts_with("X1.")
                                                                                  
))) 
### convert list of data frames to single data frame
df_tidy<-bind_rows(data, .id = "source_df")
### generate index column for detecting SF+ and SF-
df_tidy<-df_tidy %>%
  mutate(
    sf_category = case_when(
      str_detect(source_df, "SF-") ~ "sf_neg",
      str_detect(source_df, "SF+") ~ "sf_pos",
      str_detect(source_df, "cont") ~ "control",
      TRUE ~ "noneDetected" # Default case if neither is found
    )
  )
### need some smoothing for cutaneous data
df_tidy<-arrange(df_tidy,source_df,recTime) %>%
  mutate(ma2=rollapply(Spikes,2,mean,align='right',fill=NA))
### summaries
df_tidy_mean <- df_tidy %>%
  filter(!is.na(recTime)) %>%
  group_by(recTime,sf_category) %>%
  summarise(n = n(),
            mean = mean(ma2),
            median = median(ma2),
            sd = sd(ma2)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)
########################### quick visualization
cutaneousFig<-df_tidy_mean %>% filter( sf_category == "control" & n == 45 | sf_category == "sf_pos" | sf_category == "sf_neg") %>%
  ggplot( aes(x=recTime, y=mean, color = sf_category)) +
  geom_line(aes(x=recTime, y=mean,color = sf_category)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, fill=sf_category),color="grey70",alpha=0.4)+
  scale_fill_manual(values = c("#3A53A4","#A3A3CE", "#A4A4A5"))+
  scale_color_manual(values = c("#3A53A4","#A3A3CE", "#A4A4A5"))+
  theme_classic() +
  facet_wrap(~ sf_category, nrow = 3)
########################### analyses/modeling
########################### saving data
########################### saving figures
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")
ggsave(cutaneousFig, file = "cutaneous_population_code_Fig_split.pdf", width = 18, height = 18, units = "cm", path = "figures")
########################### clean up


########################### Figure 8b ###########################
########################### description
# behavior
########################### load dependencies

########################### custom functions
########################### load data
df <- read_excel("data/Copy of Acute ladder rung test.xlsx",
                 na = "NA")
########################### data wrangling
df$treat<-as.factor(df$treat)
df$rat<-as.factor(df$rat)
df_sum<-df %>% 
  group_by(treat,rat) %>%
  # filter(cutMusc == "cut") %>%  ### change from 'cut' to 'muscle'
  # filter(spontType == 1 | spontType == 2) %>%
  summarise(mean = mean(slips_total/n_steps, na.rm=T),
            sd = sd(slips_total/n_steps, na.rm=T),
            
  ) 
########################### quick visualization
stepping_Fig<-df %>% 
  ggplot(aes(x=treat, y=slips_per_step, colour = treat, fill=treat))+
  geom_jitter(alpha=1, width = .1)+
  geom_boxplot(alpha=0.8,outlier.shape = NA)+
  geom_boxplot(color = "black", fill = NA, fatten = 4, outlier.shape = NA)+
  geom_line(aes(group=rat_run))+
  scale_fill_manual(values = c("#A4A4A5", "#3A53A4"))+
  scale_color_manual(values = c("#A4A4A5", "#3A53A4"))+
  labs(title="acute effects on step accuracy",x="rat", y = "errors/step")+
  stat_n_text(
    y.pos = 0, #we can specify where in y axis the samle size should be denoted
    color = "black", #choose any color
    text.box = F, #draws a box outside the n
    size = 3
  )+
  theme_classic()

########################### analyses/modeling
########################### saving data
########################### saving figures
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")
ggsave(stepping_Fig, file = "stepping_Fig.pdf", width = 8, height = 10, units = "cm", path = "figures")
########################### clean up



########################### Figure 8c ###########################
########################### description
# behavior
########################### load dependencies

########################### custom functions
########################### load data
df <- read_excel("data/temp_behavior.xlsx",
                 na = "NA")
########################### data wrangling
df$`Cont/OX`<-as.factor(df$`Cont/OX`)
df$temp<-as.factor(df$temp)
df$file_name_date<-as.factor(df$file_name_date)
df$fil_ID_number<-as.factor(df$fil_ID_number)
########################### quick visualization
temp_Fig<-df %>% 
  ggplot(aes(x=`Cont/OX`, y=`total time`, colour = `Cont/OX`, fill=`Cont/OX`))+
  geom_jitter(alpha=1, width = .2)+
  geom_boxplot(alpha=0.8,outlier.shape = NA)+
  geom_line(aes(group=animalID))+
  scale_fill_manual(values = c("#A4A4A5", "#3A53A4"))+
  scale_color_manual(values = c("#A4A4A5", "#3A53A4"))+
  facet_wrap( ~ temp, ncol = 2)+
  geom_boxplot(color = "black", fill = NA, fatten = 4, outlier.shape = NA)+
  
  labs(title="acute effects on allodynia",x="rat", y = "foreFootTime")+
  stat_n_text(
    y.pos = 0, #we can specify where in y axis the samle size should be denoted
    color = "black", #choose any color
    text.box = F, #draws a box outside the n
    size = 3
  )+
  theme_classic()
########################### analyses/modeling
########################### saving data
########################### saving figures
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")
ggsave(temp_Fig, file = "temp_Fig.pdf", width = 12, height = 10, units = "cm", path = "figures")
########################### clean up

########################### Figure 9 ###########################
########################### description
########################### load dependencies
########################### custom functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
########################### load data
df <- read_excel("data/acute_drg_icp_ms.xlsx",
                 na = "NA")
df_short <-df %>%
  filter(calibration == "yes") %>%
  select(concentration_num, ppm)

########################### data wrangling


df_observed<-df %>% 
  group_by(concentration) %>%
  filter(calibration == "no") %>%
  dplyr::summarise(mean.ppm = mean(ppm, na.rm = TRUE),
            sd.ppm = sd(ppm, na.rm = TRUE),
            n.ppm = n()) %>%
  mutate(se.ppm = sd.ppm / sqrt(n.ppm),
         lower.ci.ppm = mean.ppm - qt(1 - (0.05 / 2), n.ppm - 1) * se.ppm,
         upper.ci.ppm = mean.ppm + qt(1 - (0.05 / 2), n.ppm - 1) * se.ppm)

df_short <- summarySE(df_short, measurevar="ppm", groupvars=c("concentration_num"))
detach("package:plyr", unload = TRUE)


########################### quick visualization
pd <- position_dodge(0.1) # move them .05 to the left and right

standardCurve_pt_ppm_fig<-ggplot(df_short, aes(x=concentration_num, y=ppm)) + 
  geom_errorbar(aes(ymin=ppm-se, ymax=ppm+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)+
  # scale_y_continuous(
  #   trans = "log2")+
  scale_y_log10(limits = c(0.001,1e5))+
  scale_x_continuous(
    trans = "log10")+
  geom_rect(data=df_observed, aes(xmin=2,
                                  xmax=16,
                                  ymin=lower.ci.ppm, 
                                  ymax=upper.ci.ppm),
            color="grey20",
            alpha=0.5,
            inherit.aes = FALSE)+
  theme_classic()
  # annotate('text', x = 1.6, y = (df_observed$mean.ppm/2),
  #          label = 'Observed (ppm)', 
  #          size = 4, 
  #          angle='90')+
  # annotate('text', x = 6, y = (df_observed$lower.ci.ppm-.15),
  #          label = 'Predicted \n Concentration \n Range (2-16uM)', 
  #          size = 4
  # )

observed_pt_ppm_fig<-df %>%
  group_by(concentration) %>%
  filter(calibration == "no") %>%
  ggplot(aes(x=as.factor(concentration), y=ppm, color=as.factor(concentration), fill=as.factor(concentration))) +
  geom_boxplot(alpha=0.4)+
  # geom_line(aes(group = experiment), size=0.2, color='black', alpha=0.6, linetype=2)+
  geom_point(aes(fill=concentration,group=sampleName),size=1,shape=21)+
  # facet_wrap(vars(affType))+
  scale_color_manual(values=c("#B6B5B5", "blue", "red", "#609B53"))+
  scale_fill_manual(values=c("#B6B5B5", "blue", "red", "#609B53"))+
  theme_classic()+
  scale_y_log10(limits = c(0.001,1e5))
########################### analyses/modeling
########################### saving data
########################### saving figures
setwd("~/Dropbox-GaTech/CoS/BioSci/BioSci-Housley_Lab/04-papers/acute_ox/acute_oxaliplatin_encoding/")
ggsave(standardCurve_pt_ppm_fig, file = "standardCurve_pt_ppm_fig.pdf", width = 10, height = 12, units = "cm", path = "figures")
ggsave(observed_pt_ppm_fig, file = "observed_pt_ppm_fig.pdf", width = 10, height = 12, units = "cm", path = "figures")

########################### clean up