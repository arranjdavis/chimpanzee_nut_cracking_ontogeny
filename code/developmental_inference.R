"
Authors: Sophie Berdugo & Arran J. Davis
Emails: sophie.berdugo@anthro.ox.ac.uk | arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
"

#load libraries
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggdist)
library(ggeffects)
library(RColorBrewer)
library(ggpubr)
library(lme4)
library(nlme)
library(sjPlot)
library(effects)
library(tidyr)
library(ggplot2)
library(extrafont)
library(performance)
library(car)
library(VGAM)
library(glmmTMB)
library(ordinal)
library(cowplot)
library(MuMIn)
library(rcompanion)
library(dotwhisker)

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load the data
peering = read.csv("../clean_data/clean_peering_data.csv")
efficiency_data = read.csv("../clean_data/cleaned_efficiency_data.csv")

### ### ###

#function for "not in"
'%!in%' = function(x,y)!('%in%'(x,y))

#subset the data
efficiency_data = subset(efficiency_data, efficiency_data$bout_outcome != "Failed" & 
                         efficiency_data$nut_species != "Coula nut" &
                         efficiency_data$bout_outcome != "None" &
                         efficiency_data$learner %!in% c("Clinging","Clinging and Peering"))

peering = subset(peering, peering$Whole.bout != "Unknown" &
                 peering$Year != "2009")

#make individual and sex a factor
efficiency_data$sex = as.factor(efficiency_data$sex)
efficiency_data$subject = as.factor(efficiency_data$subject)

################################################################################################################################################

### GRAPH THEMES ###

#fonts
quartzFonts(avenir = c("Avenir Book", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"))

#theme
header_size = 18
axis_size = 15

#theme for plots
avenir_theme = theme(text=element_text(size=header_size,family='avenir'),
                     axis.text.x = element_text(color = 'black', size = axis_size, vjust = 1),
                     axis.text.y = element_text(color = 'black', size = axis_size),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), face = "bold"),
                     axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 0, l = 0), face = "bold"),
                     panel.background = element_blank(),
                     panel.grid.major.x = element_line(color = '#e7e7e7'),
                     panel.grid.major.y = element_line(color = '#e7e7e7'),
                     legend.key = element_blank(),
                     plot.title = element_text(hjust = 0.5, face = "bold"))

#set the font (chose one)
par(family = 'avenir')

################################################################################################################################################

### VARIATION IN DEVELOPMENTAL PERIOD ###

### PEERING ###

#check the peering measures for each individual (on average)
peering_summary = peering %>% group_by(Subject) %>% dplyr::summarise(freq_peering = n())

#add in proportion of whole bouts peered at
table(peering$Subject, peering$Whole.bout)

Subject = c("Fanle", "Fanwaa", "Flanle", "Fotaiu", "Jeje", "Joya", "Juru", "Nto", "Peley", "Pokuru", "Poni", "Vuavua", "Yolo")
all_peering = c(62, 29, 61, 33, 54, 21, 44, 34, 15, 9, 14, 10, 50)
whole_bout = c(46, 14, 30, 20, 38, 17, 33, 18, 11, 2, 8, 5, 22)

proportion_whole_bout = data.frame(Subject, all_peering, whole_bout)
proportion_whole_bout = proportion_whole_bout %>% mutate(proportion_whole_bout = whole_bout / all_peering)

#correlation between number of peering events and proportion whole bout peered at
cor.test(proportion_whole_bout$all_peering, proportion_whole_bout$proportion_whole_bout, method = "pearson")

#plot correlation
ggscatter(proportion_whole_bout, x = "all_peering", y = "proportion_whole_bout", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Number of peering events", ylab = "Proportion of peering\nat the whole bout") + avenir_theme

#add `proportion_whole_bout` to the peering_summary data frame
peering_summary = merge(peering_summary, proportion_whole_bout[ , c("Subject", "proportion_whole_bout")], by = "Subject")

#critical t-value for alpha = 0.05 (two-tailed)
t_crit = qt(p = 0.025, df = (nrow(peering_summary) - 1))

t_values = (peering_summary$freq_peering - mean(peering_summary$freq_peering)) / (sd(peering_summary$freq_peering) / sqrt(nrow(peering_summary)))
peering_summary$freq_peering_tscore = t_values

t_values = (peering_summary$proportion_whole_bout - mean(peering_summary$proportion_whole_bout)) / (sd(peering_summary$proportion_whole_bout) / sqrt(nrow(peering_summary)))
peering_summary$proportion_whole_bout_tscore = t_values

#add new columns to data frame
peering = merge(peering, peering_summary[ , c("Subject", "freq_peering", "proportion_whole_bout")], by = "Subject")

### ### ###

#subset the efficiency data to just the individuals who have post-learning period efficiency 
efficiency_data_sub = subset(efficiency_data, efficiency_data$subject %in% unique(peering$Subject))

#add a column to the subsetted efficiency data frame with the values from the peering data frame

#make function to add column values according to individual chimpanzee
peerer_type <- function(x, n) { 
  if(x == "Fanle") y <- peering_summary[1, n]
  if(x == "Flanle") y <- peering_summary[3, n]
  if(x == "Fotaiu") y <- peering_summary[4, n]
  if(x == "Jeje") y <- peering_summary[5, n]
  if(x == "Joya") y <- peering_summary[6, n]
  if(x == "Nto") y <- peering_summary[8, n]
  if(x == "Peley") y <- peering_summary[9, n]
  if(x == "Poni") y <- peering_summary[11, n]
  if(x == "Vuavua") y <- peering_summary[12, n]
  if(x == "Yolo") y <- peering_summary[13, n]
  return(y)
}

#add the columns to the data set
efficiency_data_sub$proportion_whole_bout_tscore  = unlist(sapply(efficiency_data_sub$subject, peerer_type, n = "proportion_whole_bout_tscore"))

### OPPORTUNITY PROVISION ###

#load data
opportunity = read_csv("../clean_data/clean_opportunity_data.csv")
total_visibility = read_csv("../clean_data/clean_total_visibility_data.csv")

#sum the opportunity duration and total visibility duration for each subject
tapply(opportunity$Duration, opportunity$Subject, sum) 
tapply(total_visibility$Duration, total_visibility$Subject, sum) 

##get the proportion of time models were providing opportunities to their infants (opportunity duration/total visibility duration)
#make data frame
subjects = c("Fanle", "Fanwaa", "Flanle", "Fotaiu", "Jeje", "Joya", "Juru", "Pokuru", "Vuavua", "Yolo")
opportunity_duration = c(50202.363, 40936.223, 39746.375, 31703.803, 72592.377, 82297.464, 22684.814, 11020.994, 9052.595, 37556.967)
total_vis_duration = c(92821.10, 98604.13, 88664.35, 46993.72, 113040.15, 148101.84,  36988.16,  30078.41, 22441.43,  63370.87)

opp_prov = data.frame(subjects, opportunity_duration, total_vis_duration)

##create new proportion column in opportunity
opp_prov = opp_prov %>%
  mutate(opportunity_provision = opportunity_duration/total_vis_duration)

#check the median amount of time for opportunity provision for each individual (on average)
opp_prov_summary = opp_prov %>% group_by(subjects) %>% dplyr::summarise(mean_opp_prov = mean(opportunity_provision))

#critical t-value for alpha = 0.05 (two-tailed)
t_crit = qt(p = 0.025, df = (nrow(opp_prov_summary) - 1))
t_values = (opp_prov_summary$mean_opp_prov - mean(opp_prov_summary$mean_opp_prov)) / (sd(opp_prov_summary$mean_opp_prov) / sqrt(nrow(opp_prov_summary)))
opp_prov_summary$opp_prov_tscore = t_values

#add new columns to data frame
opp_prov = merge(opp_prov, opp_prov_summary[ , c("subjects", "opp_prov_tscore")], by = "subjects")

#add a column to the subsetted efficiency data frame with the opportunity values we have just calculated

#make function to add column values according to individual chimpanzee
opportunity_provision <- function(x, n) { 
  if(x == "Fanle") y <- opp_prov[1, n]
  if(x == "Flanle") y <- opp_prov[3, n]
  if(x == "Fotaiu") y <- opp_prov[4, n]
  if(x == "Jeje") y <- opp_prov[5, n]
  if(x == "Joya") y <- opp_prov[6, n]
  if(x == "Nto") y <- NA
  if(x == "Peley") y <- NA
  if(x == "Poni") y <- NA
  if(x == "Vuavua") y <- opp_prov[9, n]
  if(x == "Yolo") y <- opp_prov[10, n]
  return(y)
}

#add the columns to the data set
efficiency_data_sub$opportunity_provision_tscore = unlist(sapply(efficiency_data_sub$subject, opportunity_provision, n = "opp_prov_tscore"))

### MODEL TOLERANCE ###

#make function to add column values according to amount tolerance towards learners
tolerance <- function(x) { 
  if(x == "Model interaction") y <- 1
  if(x == "Model interaction,Proximity") y <- 2
  if(x == "None") y <- 0
  if(x == "Not Applicable") y <- 0
  if(x == "Proximity") y <- 1
  if(x == "Scrounge") y <- 1
  if(x == "Scrounge,Model interaction") y <- 2
  if(x == "Scrounge,Model interaction,Proximity") y <- 3
  if(x == "Scrounge,Proximity") y <- 2
  if(x == "Scrounge,Tool interaction,Proximity") y <- 3
  if(x == "Tool interaction,Model interaction") y <- 2
  if(x == "Tool interaction,Proximity") y <- 2
  return(y)
}

peering$tolerance_amount = sapply(peering$Model.tolerance, tolerance)

#check the median amount of tolerance for peering for each individual (on average)
tolerance_summary = peering %>% group_by(Subject) %>% dplyr::summarise(mean_tolerance = mean(tolerance_amount))

#critical t-value for alpha = 0.05 (two-tailed)
t_crit = qt(p = 0.025, df = (nrow(tolerance_summary) - 1))
t_values = (tolerance_summary$mean_tolerance - mean(tolerance_summary$mean_tolerance)) / (sd(tolerance_summary$mean_tolerance) / sqrt(nrow(tolerance_summary)))
tolerance_summary$mean_tolerance_tscore = t_values

#add new columns to data frame
peering = merge(peering, tolerance_summary[ , c("Subject", "mean_tolerance_tscore")], by = "Subject")

#make function to add column values according to individual chimpanzee
peering_tolerance <- function(x, n) { 
  if(x == "Fanle") y <- tolerance_summary[1, n]
  if(x == "Flanle") y <- tolerance_summary[3, n]
  if(x == "Fotaiu") y <- tolerance_summary[4, n]
  if(x == "Jeje") y <- tolerance_summary[5, n]
  if(x == "Joya") y <- tolerance_summary[6, n]
  if(x == "Nto") y <- tolerance_summary[8, n]
  if(x == "Peley") y <- tolerance_summary[9, n]
  if(x == "Poni") y <- tolerance_summary[11, n]
  if(x == "Vuavua") y <- tolerance_summary[12, n]
  if(x == "Yolo") y <- tolerance_summary[13, n]
  return(y)
}

#add the columns to the dataset
efficiency_data_sub$mean_tolerance_tscore = unlist(sapply(efficiency_data_sub$subject, peering_tolerance, n = "mean_tolerance_tscore"))

### MODEL INTOLERANCE ###

#make function to add column values according to amount intolerance towards learners
intolerance <- function(x) { 
  if(x == "None") y <- 0
  if(x == "Not Applicable") y <- 0
  if(x == "Physical") y <- 1
  if(x == "Vocal") y <- 1
  return(y)
}

peering$intolerance_amount = sapply(peering$Model.intolerance, intolerance)

#check the median amount of tolerance for peering for each individual (on average)
intolerance_summary = peering %>% group_by(Subject) %>% dplyr::summarise(mean_intolerance = mean(intolerance_amount))

#critical t-value for alpha = 0.05 (two-tailed)
t_crit = qt(p = 0.025, df = (nrow(intolerance_summary) - 1))
t_values = (intolerance_summary$mean_intolerance - mean(intolerance_summary$mean_intolerance)) / (sd(intolerance_summary$mean_intolerance) / sqrt(nrow(intolerance_summary)))
intolerance_summary$mean_intolerance_tscore = t_values

#add new columns to data frame
peering = merge(peering, intolerance_summary[ , c("Subject", "mean_intolerance_tscore")], by = "Subject")

#make function to add column values according to individual chimpanzee
peering_intolerance <- function(x, n) { 
  if(x == "Fanle") y <- intolerance_summary[1, n]
  if(x == "Flanle") y <- intolerance_summary[3, n]
  if(x == "Fotaiu") y <- intolerance_summary[4, n]
  if(x == "Jeje") y <- intolerance_summary[5, n]
  if(x == "Joya") y <- intolerance_summary[6, n]
  if(x == "Nto") y <- intolerance_summary[8, n]
  if(x == "Peley") y <- intolerance_summary[9, n]
  if(x == "Poni") y <- intolerance_summary[11, n]
  if(x == "Vuavua") y <- intolerance_summary[12, n]
  if(x == "Yolo") y <- intolerance_summary[13, n]
  return(y)
}

#add the columns to the dataset
efficiency_data_sub$mean_intolerance_tscore = unlist(sapply(efficiency_data_sub$subject, peering_intolerance, n = "mean_intolerance_tscore"))

#make function to add a 'mother' column values according to individual chimpanzee
mother <- function(x) { 
  if(x == "Fanle") y <- "Fana"
  if(x == "Flanle") y <- "Fana"
  if(x == "Fotaiu") y <- "Fana"
  if(x == "Jeje") y <- "Jire"
  if(x == "Joya") y <- "Jire"
  if(x == "Nto") y <- "Nina"
  if(x == "Peley") y <- "Pama"
  if(x == "Poni") y <- "Pama"
  if(x == "Vuavua") y <- "Velu"
  if(x == "Yolo") y <- "Yo"
  return(y)
}

efficiency_data_sub$mother <- sapply(efficiency_data_sub$subject, mother)

#drop levels of `subject` where there is no data
efficiency_data_sub$subject = droplevels(efficiency_data_sub$subject)

################################################################################################################################################

### ### PLOTTING ### ###

### T DISTRIBUTIONS ###

#peering whole bout
peering_t_distribution = ggplot(data.frame(x = c(-6, 6)), aes(x = x)) +
                              stat_function(fun = dt, args = list(df = 12), colour = "black", size = 1.1) +
                              scale_y_continuous(name = "Density") +
                              scale_x_continuous(name = "", limits = c(-8.5, 8.5), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
                              geom_vline(xintercept = c(-2.178813, 2.178813), colour ="black", size = 0.8) +
                              geom_vline(xintercept = c(3.4981587, -2.2151027, -2.0157234,  0.5029504,  2.6553822,  4.9880670,  3.6759318, 
                                                        -1.1866868,  3.3085340, -7.9583334, -0.2604737, -1.8350359, -3.1576682), 
                                         colour = "#5ab4ac", linetype = "longdash", size = 1) +
                              ggtitle(expression(paste("Proportion of peering at whole bout ", italic("t"), " distribution (df = 12)"))) +
                              avenir_theme

#opportunity provision
opportunity_t_distribution = ggplot(data.frame(x = c(-6, 6)), aes(x = x)) +
                              stat_function(fun = dt, args = list(df = 9), colour = "black", size = 1.1) +
                              scale_y_continuous(name = "") +
                              scale_x_continuous(name = "", limits = c(-8.5, 8.5), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
                              geom_vline(xintercept = c(-2.262157, 2.262157), colour ="black", size = 0.8) +
                              geom_vline(xintercept = c(-4.5943735955292, -3.52481285733698, -3.18439167287387, -2.2263895101902, 0.451116519690008, 0.880077005783734, 1.94943778520268,
                                                        2.54659016852966,  3.38198783393916,    4.320758322785), 
                                         colour = "#d8b365", linetype = "longdash", size = 1) +
                              ggtitle(expression(paste("Opportunity provision ", italic("t"), " distribution (df = 9)"))) +
                              avenir_theme

#model tolerance
tolerance_t_distribution = ggplot(data.frame(x = c(-8.5, 8.5)), aes(x = x)) +
                              stat_function(fun = dt, args = list(df = 12), colour = "black", size = 1.1) +
                              scale_y_continuous(name = "Density") +
                              scale_x_continuous(name = "t-value", limits = c(-8.5, 8.5), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
                              geom_vline(xintercept = c(-2.178813, 2.178813), colour ="black", size = 0.8) +
                              geom_vline(xintercept = c(1.8756956,  2.0853545,  0.6895761, -5.7336525,  0.5808577, -7.4471376,  0.5491265,
                                                        -0.0556330, -1.1643587,  7.2126799, -1.1643587,  0.7204750,  1.8513752), 
                                         colour = "#67a9cf", linetype = "longdash", size = 1) +
                              ggtitle(expression(paste("Model tolerance to peering ", italic("t"), " distribution (df = 12)"))) +
                              avenir_theme

#model intolerance
intolerance_t_distribution = ggplot(data.frame(x = c(-8.5, 8.5)), aes(x = x)) +
                              stat_function(fun = dt, args = list(df = 12), colour = "black", size = 1.1) +
                              scale_y_continuous(name = "") +
                              scale_x_continuous(name = "t-value", limits = c(-8.5, 8.5), breaks = c(-8, -6, -4, -2, 0, 2, 4, 6, 8)) +
                              geom_vline(xintercept = c(-2.178813, 2.178813), colour ="black", size = 0.8) +
                              geom_vline(xintercept = c(0.5167710,  2.1621343,  2.4622972, -3.2278049,  2.6465242, -4.3771301,  0.7948335, -4.3771301,
                                                        0.6799010,  8.2654477, -4.3771301, -0.5843568, -0.5843568), 
                                         colour = "#ef8a62", linetype = "longdash", size = 1) +
                              ggtitle(expression(paste("Model intolerance to peering ", italic("t"), " distribution (df = 12)"))) +
                              avenir_theme

#plot all four graphs in a grid
development_variation_grid = plot_grid(peering_t_distribution, opportunity_t_distribution, tolerance_t_distribution, intolerance_t_distribution, labels = "AUTO")
development_variation_grid

#save plot
ggsave("../plots/development_variation_grid.jpg", development_variation_grid, height = 12, width = 19)


### PROPORTION WHOLE BOUTS PEERED AT ###

##whole and partial peering by age

#get data
table(peering$Whole.bout, peering$Age)

#make data frame
age = c("0", "1", "2", "3", "4", "5")
whole_peering = c(29, 27, 58, 30, 85, 35)
partial_peering = c(29, 25, 42, 15, 45, 16)
total = c(58, 52, 100, 45, 130, 51)
whole_bout_data = data.frame(age, whole_peering, partial_peering, total)

whole_bout_data = whole_bout_data %>% mutate(proportion_whole_bout = whole_peering / (whole_peering + partial_peering))

#make data frame long
whole_bout_data = gather(whole_bout_data, whole_bout, value, whole_peering, partial_peering)

#plot data
participant_n_title = expression(paste(bold("Total "), bolditalic("n"), bold(":")))

whole_bout_by_age_barchart = ggplot(whole_bout_data, aes(fill = whole_bout, x = age, y = value)) + 
                              geom_bar(position="fill", stat="identity", alpha = 0.8) +
                              scale_fill_manual(values = c("#f16913", "#4292c6"), name = "Whole bout", labels = c("No", 'Yes')) +
                              scale_x_discrete(breaks = seq(0, 5, 1), labels = c("0–1", "1–2", "2–3", "3–4", "4–5", "5–6")) +
                              scale_y_continuous(breaks = seq(0, 1, .1)) +
                              ylab("Proportion of peering events") +
                              xlab("Learner age range (years)") +
                              geom_text(aes(label = total, y = 1.05), size = 6) +
                              annotate("text", x = 0.65, y = 1.05, label = participant_n_title, size = 6, parse = TRUE) +
                              avenir_theme
whole_bout_by_age_barchart

#save plot
ggsave("../plots/whole_bout_by_age_barchart.jpg", whole_bout_by_age_barchart, height = 9, width = 16)

################################################################################################################################################

### ### MODELLING ### ###

### BOUT DURATION ###

#full model (model is too complex; model not converging)
bout_d_full = lmer(log(bout_duration) ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + (1 | mother) + (1 | subject), data = efficiency_data_sub)
summary(bout_d_full)

#model with only `subject` as a random variable (model not converging)
bout_d_subject = lmer(log(bout_duration) ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore +  (1 | subject), data = efficiency_data_sub)
summary(bout_d_subject)

#check multicollienarity (multicollienarity detected)
vif(bout_d_subject)

#full linear model (but with sex removed because of multicollinearity)
bout_d_lm = lm(log(bout_duration) ~ 1 + age + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, data = efficiency_data_sub)
summary(bout_d_lm)
AIC(bout_d_lm)

#check multicollienarity (no multicollienarity detected)
vif(bout_d_lm)

#tabulate the model
tab_model(bout_d_lm)

################################################################################################################################################

### STRIKES PER NUT ###

spn_poisson = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + 
                        (1 | subject) + (1 | mother), data = efficiency_data_sub, family = truncated_poisson(link="log")) 
summary(spn_poisson)
MuMIn::AICc(spn_poisson)

#check over dispersion (overdispersion detected)
check_overdispersion(spn_poisson)

#remove mother term to see if that removes the over dispersion
spn_poisson_subject = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + 
                        (1 | subject), data = efficiency_data_sub, family = truncated_poisson(link="log")) 
summary(spn_poisson_subject)
MuMIn::AICc(spn_poisson_subject)

#check over dispersion (overdispersion detected)
check_overdispersion(spn_poisson_subject)

#negative binomial model to account for over dispersion
spn_nb = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + 
                (1 | subject), data = efficiency_data_sub, family = truncated_nbinom1(link="log")) 
summary(spn_nb)
MuMIn::AICc(spn_nb)

#full linear model
spn_lm = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, 
                 data = efficiency_data_sub, family = truncated_nbinom1(link="log"))
summary(spn_lm)
MuMIn::AICc(spn_lm)

anova(spn_nb, spn_lm)
tab_model(spn_lm)

################################################################################################################################################

### SUCCESS RATE ###

#load and subset the data
sr_data = read.csv('../clean_data/cleaned_efficiency_data.csv')
sr_data = subset(sr_data,
                 sr_data$nut_species != "Coula nut" &
                 sr_data$bout_outcome != "None" &
                 sr_data$learner %!in% c("Clinging","Clinging and Peering"))

#make individual and sex a factor
sr_data$sex = as.factor(sr_data$sex)
sr_data$subject = as.factor(sr_data$subject)

#subset the success rate data to just the individuals who have post-learning period efficiency 
sr_data_sub = subset(sr_data, sr_data$subject %in% unique(peering$Subject))

#add the predictor columns to the data set
sr_data_sub$proportion_whole_bout_tscore  = unlist(sapply(sr_data_sub$subject, peerer_type, n = "proportion_whole_bout_tscore"))
sr_data_sub$opportunity_provision_tscore = unlist(sapply(sr_data_sub$subject, opportunity_provision, n = "opp_prov_tscore"))
sr_data_sub$mean_tolerance_tscore = unlist(sapply(sr_data_sub$subject, peering_tolerance, n = "mean_tolerance_tscore"))
sr_data_sub$mean_intolerance_tscore = unlist(sapply(sr_data_sub$subject, peering_intolerance, n = "mean_intolerance_tscore"))

#add mother column
sr_data_sub$mother <- sapply(sr_data_sub$subject,mother)

#order the outcome variable
sr_data_sub$bout_outcome_ordered = ordered(sr_data_sub$bout_outcome, c("Failed", "Smash", "Successful"))

#run multilevel model                                                                                                                
success_rate = clmm(bout_outcome_ordered ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + 
                      (1 | subject) + (1 | mother), data = sr_data_sub) 
summary(success_rate)
MuMIn::AICc(success_rate)

#run full linear model
success_rate_lm = clm(bout_outcome_ordered ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, data = sr_data_sub)
summary(success_rate_lm)
MuMIn::AICc(success_rate_lm)

#compare model fit (use model lower AIC)
anova(success_rate, success_rate_lm)
tab_model(success_rate_lm)

################################################################################################################################################

### DISPLACEMENT RATE ###

#load and subset the data
dr_data = read.csv('../clean_data/cleaned_efficiency_data.csv')
dr_data = subset(dr_data, dr_data$nut_species != "Coula nut" &
                   dr_data$bout_outcome != "None" &
                   dr_data$learner %!in% c("Clinging","Clinging and Peering"))

#subset the displacement rate data to just the individuals who have post-learning period efficiency 
dr_data_sub = subset(dr_data, dr_data$subject %in% unique(peering$Subject))

#add the predictor columns to the data set
dr_data_sub$proportion_whole_bout_tscore  = unlist(sapply(dr_data_sub$subject, peerer_type, n = "proportion_whole_bout_tscore"))
dr_data_sub$opportunity_provision_tscore = unlist(sapply(dr_data_sub$subject, opportunity_provision, n = "opp_prov_tscore"))
dr_data_sub$mean_tolerance_tscore = unlist(sapply(dr_data_sub$subject, peering_tolerance, n = "mean_tolerance_tscore"))
dr_data_sub$mean_intolerance_tscore = unlist(sapply(dr_data_sub$subject, peering_intolerance, n = "mean_intolerance_tscore"))

#add mother column
dr_data_sub$mother = sapply(dr_data_sub$subject,mother)

#run model
displacement_rate_poisson = glmmTMB(displacement_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + (1 | subject) + (1 | mother), 
                            data = dr_data_sub, family = poisson, zi = ~ 1) 
summary(displacement_rate_poisson)
MuMIn::AICc(displacement_rate_poisson)

#check for over dispersion (overdispersion detected)
check_overdispersion(displacement_rate_poisson)

#fit negative binomial model to account for over dispersion (model did not converge)
displacement_rate = glmmTMB(displacement_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore + (1 | subject) + (1 | mother), 
                                    data = dr_data_sub, family = nbinom1, zi = ~ 1) 
summary(displacement_rate)
MuMIn::AICc(displacement_rate)

#full linear model
displacement_rate_lm = glmmTMB(displacement_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, 
                               data = dr_data_sub, family = nbinom1, zi = ~ 1) 
summary(displacement_rate_lm)
MuMIn::AICc(displacement_rate_lm)

tab_model(displacement_rate_lm)

### TABULATE ALL MODELS ###

tab_model(bout_d_lm, spn_lm, success_rate_lm, displacement_rate_lm)

################################################################################################################################################

### ### PLOT MODEL OUTPUTS ### ###

#plot coefficients of full models
dw = dwplot(list(bout_d_lm, displacement_rate_lm, spn_lm, success_rate_lm), 
            dodge_size = 0.5, 
            dot_args = list(size = 2.5),
            whisker_args = list(size = 1))


#make plot pretty
model_dotwhisker_plot = dw + 
                        geom_vline(xintercept=0, lty=2) + 
                        scale_color_manual(values = c("#238b45", "#6a51a3", "#d94801", "#2171b5"), 
                                           labels = c("Success rate", "Strikes per nut", "Displacement rate", "Bout duration"),
                                           name = "Model") +
                        scale_x_continuous(breaks = seq(-2, 1, .25), limits = c(-2, 1)) +
                        scale_y_discrete(labels = c("Sex (male)", "Model intolerance", "Model tolerance", "Opportunity provision", "Proportion whole bout\n peered at", "Age")) + 
                        xlab("Coefficient estimate") + 
                        ylab("Predictor") +
                        avenir_theme 

model_dotwhisker_plot

#save plot
ggsave("../plots/model_dotwhisker_plot.jpg", model_dotwhisker_plot, height = 9, width = 16)

################################################################################################################################################

### ### MODEL ASSUMPTIONS ### ###

#load libraries
library(DHARMa)
library(broom)
library(broom.mixed)
library(sure)
library(car)
library(rcompanion)
library(influence.ME)
library(ggfortify)

### BOUT DURATION ##

#run the model
bout_d_lm = lm(log(bout_duration) ~ 1 + age + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, data = efficiency_data_sub)

#multicollinearity
vif(bout_d_lm)

#residual normality
par(mfrow=c(1,2))
plotNormalHistogram(residuals(bout_d_lm), xlab = "Log bout duration residuals") 
qqnorm(residuals(bout_d_lm), ylab="Sample quantiles for residuals")
qqline(residuals(bout_d_lm), col="blue")
par(mfrow=c(1,1))

#homoskedasticity
plot(fitted(bout_d_lm), resid(bout_d_lm)^2,
     ylab = "Squared log bout duration residuals",
     xlab = "Fitted log bout duration values")

par(mfrow=c(2,2))
plot(bout_d_lm)
par(mfrow=c(1,1))

#influence (procedure was taken from https://towardsdatascience.com/identifying-outliers-in-linear-regression-cooks-distance-9e212e9136a)
cooksD = cooks.distance(bout_d_lm)
influential = cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential
names_of_influential = names(influential)
outliers = efficiency_data_sub[names_of_influential,]

#check the bout duration of the influential observations
View(outliers)

###  ###  ###  ###  ###  ###  ###  ### 

### SUCCESS RATE ##

#residual normality (extract surrogate residuals)
surrogate(success_rate_lm, method=c("latent"), nsim=1L)
set.seed(1225)
autoplot.clm(success_rate_lm, nsim=100,what="qq")

#homoskedasticity (residual vs fitted surrogate residuals)
set.seed(1225)
autoplot.clm(success_rate_lm, nsim=100,what="fitted",alpha=0.5) 

###  ###  ###  ###  ###  ###  ###  ### 

### STRIKES PER NUT ###

#run the model
spn_lm = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, 
                 data = efficiency_data_sub, family = truncated_nbinom1(link="log"))

#normality of residuals and homoskedasticity
spn_simres = simulateResiduals(spn_lm)
plot(spn_simres)
testDispersion(spn_simres)
testResiduals(spn_simres)
testOutliers(spn_simres, type = "bootstrap")

#top eight scaled residuals are the locations of the outliers
spn_outlier = order(spn_simres$scaledResiduals, decreasing = TRUE)
spn_outlier

#check the observation where the simulated residuals are outliers
mean(efficiency_data_sub$strike_count)
efficiency_data_sub[ c(14, 199, 239, 242, 323, 402, 461, 670), ]

#remove outliers from models
spn_lm_no_outliers = glmmTMB(strike_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, 
                             data = efficiency_data_sub[-c(14, 199, 239, 242, 323, 402, 461, 670),], family = truncated_nbinom1(link="log"))

summary(spn_lm_no_outliers)

###  ###  ###  ###  ###  ###  ###  ### 

### DISPLACEMENT RATE ###

#run the model
displacement_rate_lm = glmmTMB(displacement_count ~ 1 + age + sex + proportion_whole_bout_tscore + opportunity_provision_tscore + mean_tolerance_tscore + mean_intolerance_tscore, 
                               data = dr_data_sub, family = nbinom1, zi = ~ 1) 

#normality of residuals and homoskedasticity
dr_simres = simulateResiduals(displacement_rate_lm)
plot(dr_simres)
testDispersion(dr_simres)
testResiduals(dr_simres)

################################################################################################################################################

### PACKAGE REFERENCE LIST ###

#create function
citations <- function(includeURL = TRUE, includeRStudio = TRUE) {
  if(includeRStudio == TRUE) {
    ref.rstudio <- RStudio.Version()$citation
    if(includeURL == FALSE) {
      ref.rstudio$url <- NULL;
    }
    print(ref.rstudio, style = 'text')
    cat('\n')
  }
  
  cit.list <- c('base', names(sessionInfo()$otherPkgs))
  for(i in 1:length(cit.list)) {
    ref <- citation(cit.list[i])
    if(includeURL == FALSE) {
      ref$url <- NULL;
    }
    print(ref, style = 'text')
    cat('\n')
  }
}

#print all citations
citations()
