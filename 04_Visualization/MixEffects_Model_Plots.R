#Plotting mean RT, RSA, and mixed effects model results

Dir_flash <- "/Volumes/RESEARCH/Flanker"
# Load libraries
library(data.table)
library(dplyr)
library(broom)
library(rhdf5)
library(caret)
library(foreach)
library(doMC)
library(tidyr)
library(binhf)
library(GGally)
library(RColorBrewer)
library(RSQLite)
library(lazyeval)
library(stringr)
library(psych)
library(lme4)
library(RcppEigen)
library(rio)
library(gridExtra)

# non analysis related libraries
library(parallel)
# Source Files
setwd(Dir_R)
source('basic_lib.R')


#**PLOTS FOR RSM vector RT PREDICTABILITY ON CONGRUNET AND INCONGRUENT TRIALS**

# figpath <- tcltk::tk_choose.dir()
# if(!is.na(figpath)) setwd(figpath);getwd()
#dir.create("Plots")
filename1 <- paste0(Dir_flash,"/rIncongruent.rds")
filename2 <- paste0(Dir_flash,"/rCongruent.rds")
#results1 <- readRDS(filename1)
load(filename1)
load(filename2)

setwd("/Users/tesufuai/Documents/Flanker/Plots")

ab=c("a","b")
plottitle <- list(bquote(bold("a) ")), #~ "Congruent trial RSA fits predicting RT"),
                  bquote(bold("b) "))#~ "Inongruent trial RSA fits predicting RT")
)
#a <- get(paste0('results',as.character(g))) %>% dplyr::mutate(RSA_Model = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))
#congruent,noconflict
a <- results2 %>% dplyr::mutate(`RSA Model Vector` = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))
#incongruent, conflict
b <- results1 %>% dplyr::mutate(`RSA Model Vector` = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))#incongruent

#named list to store plots
plots <- setNames(vector("list",length=2),paste0("rsa_model_tval_",ab[1:2]))

for (g in 1:2){
  #plot tvalues for each RSA model across the span of a trial
  plots[[g]] <- ggplot(get(ab[g]), aes(x = time, y = `t value`, group = `RSA Model Vector`)) +
    #facet_wrap(~prevTrial+currentTrial, nrow=2) +
    #facet_wrap(~prevTrial, nrow=2) +
    geom_line(aes(color = `RSA Model Vector`), linewidth=1) +
    geom_ribbon(aes(ymin = `t value` - `Std. Error`, ymax = `t value` + `Std. Error`, fill = `RSA Model Vector`), alpha = .2) +
    #coord_cartesian(ylim=c(-7,4)) +
    geom_hline(yintercept=0,linetype="longdash",linewidth=0.5,color="lightgray")+
    geom_vline(xintercept=0,linetype=1,linewidth=1)+
    geom_vline(xintercept=433,linetype=1,linewidth=1)+
    geom_vline(xintercept=133,linetype="longdash",linewidth=.5, color = "darkgray")+
    scale_color_manual(values = c("goldenrod2", "dodgerblue2", "red2", "green4")) + #goldenrod2 = flanker, dodgerblue2 = target, red2 = conjunction, green4 = congruency, mediumpurple1 = TF confusion, red2 = unbound flanker, blue3 = unbound target
    scale_fill_manual(values = c("goldenrod2", "dodgerblue2", "red2", "green4")) + #goldenrod2 = flanker, dodgerblue2 = target, red2 = conjunction, green4 = congruency, mediumpurple1 = TF confusion, red2 = unbound flanker, blue3 = unbound target
    scale_x_continuous(breaks=seq(-200,1000,by=100),limits = c(-100, ceiling((433 + mean(ds_cj$subRT) + sd(ds_cj$subRT))/100) * 100), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-6,6,by=2))+
    coord_cartesian(ylim=c(-6.5,6),xlim=c(433,1000))+
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=-1.96,ymax=1.96,fill="slategray2",alpha=0.5)+
    #labs(title = paste0("RSA fits, ",freqband, " band"))
    labs(title = plottitle[[g]],
         x="Time (ms)",
         y="t-value") +
    theme_classic() +
    theme(axis.title.x = element_text(size = 16),   # x-axis label
          axis.title.y = element_text(size = 16),
          legend.box.background = element_rect(colour = NA))
}
dev.new(width=8,height=13)
grid.arrange(plots$rsa_model_tval_a,plots$rsa_model_tval_b,ncol=1)




#}

#**PLOTS FOR RT PREDICTABILITY ON n-1 CONGRUNET AND INCONGRUENT TRIALS**
#*

load("rpCongruent.rds")
load("rpIncongruent.rds")

cd=c("c","d")
plottitleprev <- list(bquote(bold("a) ") ~ "Congruent previous trial RSA fits predicting RT"),
                      bquote(bold("b) ") ~ "Incongruent previous trial RSA fits predicting RT")
)
#a <- get(paste0('results',as.character(g))) %>% dplyr::mutate(RSA_Model = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))
a <- results4 %>% dplyr::mutate(`RSA Model Vector` = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))#congruent
b <- results3 %>% dplyr::mutate(`RSA Model Vector` = factor(predictor, levels = c("FLANK_RSA", "TARG_RSA","CONJ_RSA", "CONG_RSA"), labels = c("Flanker ID", "Target ID","Conjunction","Congruency")))#incongruent

#named list to store plots
plotsprev <- setNames(vector("list",length=2),paste0("rsa_model_tval_",cd[1:2]))

for (g in 1:length(cd)){
  #plot tvalues for each RSA model across the span of a trial
  plotsprev[[g]] <- ggplot(get(ab[g]), aes(x = time, y = `t value`, group = `RSA Model Vector`)) +
    #facet_wrap(~prevTrial+currentTrial, nrow=2) +
    #facet_wrap(~prevTrial, nrow=2) +
    geom_line(aes(color = `RSA Model Vector`), linewidth=1) +
    geom_ribbon(aes(ymin = `t value` - `Std. Error`, ymax = `t value` + `Std. Error`, fill = `RSA Model Vector`), alpha = .2) +
    #coord_cartesian(ylim=c(-7,4)) +
    geom_hline(yintercept=0,linetype="longdash",linewidth=0.5,color="lightgray")+
    geom_vline(xintercept=0,linetype=1,linewidth=1)+
    geom_vline(xintercept=433,linetype=1,linewidth=1)+
    geom_vline(xintercept=133,linetype="longdash",linewidth=.5, color = "darkgray")+
    scale_color_manual(values = c("goldenrod2", "dodgerblue2", "red2", "green4")) + #goldenrod2 = flanker, dodgerblue2 = target, red2 = conjunction, green4 = congruency, mediumpurple1 = TF confusion, red2 = unbound flanker, blue3 = unbound target
    scale_fill_manual(values = c("goldenrod2", "dodgerblue2", "red2", "green4")) + #goldenrod2 = flanker, dodgerblue2 = target, red2 = conjunction, green4 = congruency, mediumpurple1 = TF confusion, red2 = unbound flanker, blue3 = unbound target
    scale_x_continuous(breaks=seq(-200,1000,by=100),limits = c(-100, ceiling((433 + mean(ds_cj$subRT) + sd(ds_cj$subRT))/100) * 100), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-6,6,by=2))+
    coord_cartesian(ylim=c(-6.5,6),xlim=c(433,1000))+
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=-1.96,ymax=1.96,fill="slategray2",alpha=0.5)+
    #labs(title = paste0("RSA fits, ",freqband, " band"))
    labs(title = plottitleprev[[g]]) +
    theme_classic() +
    theme(legend.box.background = element_rect(colour = "black"))
}
dev.new()
grid.arrange(plotsprev$rsa_model_tval_c,plotsprev$rsa_model_tval_d,ncol=1)

#Check RT for interpretation/max cut off 
ds_cj_RTfromtrialonset <- ds_cj %>% mutate(RT_from_trial_onset=433+subRT)
#calculate ratio of categorized RT windows 
nrow(ds_cj_RTfromtrialonset[ds_cj_RTfromtrialonset$RT_from_trial_onset>900])/nrow(ds_cj_RTfromtrialonset)

nrow(ds_cj_RTfromtrialonset[ds_cj_RTfromtrialonset$RT_from_trial_onset>950])/nrow(ds_cj_RTfromtrialonset)

nrow(ds_cj_RTfromtrialonset[ds_cj_RTfromtrialonset$RT_from_trial_onset>1000])/nrow(ds_cj_RTfromtrialonset)

nrow(ds_cj_RTfromtrialonset[ds_cj_RTfromtrialonset$RT_from_trial_onset>1050])/nrow(ds_cj_RTfromtrialonset)

#rewrite shit lost from the past day: Mean RSA plot without CTR vector
#new stuff to plot: line graph with RT effects for congruency and n-1 congruency condtions, 
#Bar graph with mean Flanker t-values for each congruency condition and n-1 congruency conditions










# 1) Add Congruency label once
ds_long <- dsA %>%
  mutate(Congruency = if_else(ConfTrial == 1, "Incongruent", "Congruent")) %>%
  # 2) Aggregate at the grain you actually intend (collapse duplicates!)
  group_by(SUBID,time, Congruency) %>%
  summarize(
    `Target vector`      = mean(TARG_RSA,  na.rm = TRUE),
    `Flanker vector`     = mean(FLANK_RSA, na.rm = TRUE),
    `Conjunction vector` = mean(CONJ_RSA,  na.rm = TRUE),
    `Congruency vector`  = mean(CONG_RSA,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 3) Wide -> long
  pivot_longer(
    cols = c(`Target vector`,`Flanker vector`,`Conjunction vector`,`Congruency vector`),
    names_to = "RSA_Model",
    values_to = "RSA_Value"
  )

# Group and summarize
ICRSA_summary <- ds_long %>% 
  group_by(Congruency,RSA_Model,time) %>%
  summarise(
    mean_RSA = mean(RSA_Value, na.rm = TRUE),
    se = sd(RSA_Value, na.rm = TRUE) / sqrt(n()),  # Optional: standard error
    .groups = "drop"
  ) %>% rename(`Mean t-value`=mean_RSA, Time=time
  ) %>% mutate(`t-value`=as.numeric(`Mean t-value`))





#ICRSA_summary <- ICRSA_summary %>% rename(`Mean t-value`=mean_RSA, Time=time
#                                                         ) %>% mutate(`Mean t-value`=as.numeric(`Mean t-value`))



#plot
ggplot(ICRSA_summary, aes(x = Time, y = `t-value`, group = interaction(Congruency, RSA_Model))) +
  geom_line(aes(color = RSA_Model,linetype=Congruency), linewidth = 1) +
  geom_ribbon(aes(
    ymin = `t-value` - se,
    ymax = `t-value` + se,
    fill = RSA_Model,
    group=interaction(Congruency,RSA_Model)
  ), alpha = .3, position="identity") +
  geom_hline(yintercept = 0, linetype = "longdash", linewidth = 0.5, color = "lightgray") +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 1) +
  geom_vline(xintercept = 433, linetype = 1, linewidth = 1) +
  geom_vline(xintercept = 133, linetype = "longdash", linewidth = .5, color = "darkgray") +
  scale_color_manual(values = c("goldenrod2","red2", "dodgerblue2", "green4")) +
  scale_fill_manual(values = c("goldenrod2", "red2","dodgerblue2", "green4")) +
  scale_x_continuous(
    breaks = seq(-200, 1000, by = 100),
    limits = c(0, ceiling((433 + mean(ds_cj$subRT) + sd(ds_cj$subRT)) / 100) * 100),
    expand = c(0, 0)
  ) +
  #scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  coord_cartesian(xlim = c(0, 1000)) +
  #annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96, fill = "slategray2", alpha = 0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16),   # x-axis label
        axis.title.y = element_text(size = 16),
        legend.box.background = element_rect(colour = NA))

#save(MeanRSAtvalueplot,file="MeanRSAtvalueplot.png")








#Bar graph for incongruent Flanker t- values as a function of n-1 congruency conditions
prevCong <- dsA %>% mutate(`Previous Trial Congruency`=if_else(prevConfTrial==1,"Incongruent","Congruent"),
                           `Current Trial Congruency` =if_else(ConfTrial==1,"Incongruent","Congruent")) 
FlankCong <- prevCong %>% 
  filter(time>133 &time<432) %>% #prevCong[prevCong$time>133&prevCong$time<432] %>% 
  group_by(SUBID,BLOCK,TRIAL,`Previous Trial Congruency`) %>%
  summarise(mean_flanker_t=mean(FLANK_RSA,na.rm=TRUE),
            .groups ="drop") %>% 
  #rename(`Flanker vector`     = FLANK_RSA) %>% 
  # 3) Wide -> long
  group_by(`Previous Trial Congruency`) %>% 
  summarise(`Mean t-value`=mean(mean_flanker_t),
            N            =n(),
            SE           =sd(mean_flanker_t)/sqrt(N),
            CI95         =SE*qt(0.975,df=N-1),
            .groups       ="drop")

#t-test
submean <- prevCong %>% 
  #filter(time>133 &time<432) %>% #prevCong[prevCong$time>133&prevCong$time<432] %>% 
  group_by(SUBID,BLOCK,TRIAL,`Previous Trial Congruency`) %>%
  summarise(mean_flanker_t=mean(FLANK_RSA,na.rm=TRUE),
            .groups ="drop") 

t.test(mean_flanker_t~`Previous Trial Congruency`,data=submean)



#plot
ggplot(FlankCong,aes(`Previous Trial Congruency`,`Mean t-value`,fill=`Previous Trial Congruency`)) +
  geom_col(width=0.5)+
  geom_errorbar(aes(ymin=`Mean t-value`-CI95,ymax=`Mean t-value`+CI95),linetype="longdash",width=0.1)+
  scale_fill_manual(values=c("green4","green4"))+
  labs(
    x="Previous Trial Congruency",
    y="t-value") +
  theme_classic() +
  theme(
    #plot.title   =element_text(size=18), 
    axis.title.x = element_text(size = 16),   # x-axis label
    axis.title.y = element_text(size = 16),
    axis.text.x  =element_text(size=12),
    legend.position = "none")    # y-axis label










# 
#   pivot_longer(
#     cols = `Flanker vector`,
#     names_to = "RSA_Model",
#     values_to = "RSA_Value") %>% 
#     summarize(`Mean t-value`=mean(RSA_Value,na.rm=TRUE),
#               t-val_se=sd(RSA_Value,na.rm=TRUE)/sqrt(n()),
#              .groups = "drop") 


#plot
# ggplot(subset(prevCong, `Current Trial Congruency`=="Incongruent"),aes(x=`Previous Trial Congruency`,y=`Mean t-value`,fill=`Current Trial Congruency`)) +
#          geom_col(width=0.5)+
#          geom_errorbar(aes(ymin=`Mean t-value`-se,ymax=`Mean t-value`+se),width=0.1)+
#          labs(title = "Flanker Representations for Incongruent Current Trials ") +
#          scale_fill_manual(values="green4")+
#          theme_minimal()+
#          theme(legend.position = "none") 


#RT effects for n-1 congruency conditions
#within subject means for each prev cur congruency combination
RTCong <- prevCong %>%
  group_by(SUBID, `Previous Trial Congruency`,`Current Trial Congruency`) %>% 
  summarize(`RT(ms)`=mean(RT,na.rm=TRUE),
            .groups = "drop")

#group mean RT and between subject CI for congruency conditions
dsRT <- RTCong %>% 
  group_by(`Previous Trial Congruency`,`Current Trial Congruency`) %>% 
  summarize(`Mean RT(ms)`=mean(`RT(ms)`),
            SERT     =sd(`RT(ms)`)/sqrt(n()),
            CI95RT   =SERT*qt(0.975,df=n()-1),
            .groups="drop")

ggplot(dsRT,
       aes(x=`Previous Trial Congruency`,
           y=`Mean RT(ms)`,
           color=`Current Trial Congruency`,
           group=`Current Trial Congruency`))+
  geom_line(position=position_dodge(width=0),linewidth=0.7)+
  geom_point(position = position_dodge(width = 0),size=2.5)+
  geom_errorbar(aes(ymin=`Mean RT(ms)`-CI95RT,
                    ymax=`Mean RT(ms)`+CI95RT),
                position=position_dodge(width=0),width=0.5,
                linetype="longdash")+
  labs( 
    x="Previous Trial Congruency",
    y="Mean RT(ms)",
    color="Current Trial")+
  theme_classic()+
  theme(
    #plot.title   =element_text(size=18), 
    axis.title.x = element_text(size = 16),   # x-axis label
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size=12),
    legend.title = element_text(size=15),
    legend.text = element_text(size=12))# y-axis label
#legend.position = "nonet")    
