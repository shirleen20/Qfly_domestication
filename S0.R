#______________________________________24/01/2025_______________________________

# Here we calculate fly weight weight differences and plot dry weight vs total lipids for 
# the seven strains to generate Figure S1, Table S3, and Dataset S1.

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("gtools")
library("Rmisc")
library("gridExtra")

rm(list=ls())

#_________________Outlier function to be used, when required______________________________ 
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#__________________Create a list of ambiguous lipid species from Shirleen's csv file______
ALS <- read_csv("data/Lipids to filter ploty.csv") %>% 
  dplyr::select(3) %>% 
  unique() %>% 
  c()
#______________ Load samples and remove outliers from batch02 ______________________________
# correct names of incorrectly labeled sample names
'%!in%' <- function(x,y)!('%in%'(x,y)) # a function to create opposite of %in% for filtering out

B01 <- read_csv("data/CD_results_Batch01_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))%>% 
  dplyr::filter(samples != "ct_n09")

B02 <- read_csv("data/CD_results_Batch02_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "S06_95") %>% 
  dplyr::mutate(samples = recode(samples,  syd_n87 = "syd_o87",syd_n104 = "syd_o104"))  

B03 <- read_csv("data/CD_results_Batch03_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>%
  dplyr::filter(samples != "syd_n134") %>% # remove outlier
  dplyr::mutate(samples = recode(samples, ct_n161 = "ct_o161", S06_171 = "s06_171",
                                 syd_160 = "syd_o160", syd_n138 = "syd_o138", syd_n169 = "syd_o169")) %>% 
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

B04 <- read_csv("data/CD_results_Batch04_w_QCs_library_filtered.csv") %>%
  dplyr::select(-c(pqn, sumpqn, `Area (Max.)`)) %>% 
  dplyr::filter(samples != "s06_217") %>% # remove outlier
  dplyr::filter(Name %!in% c(ALS$`Lipid Species`))

#_________________ calculate total area of all peaks in each sample and calculate their mean ___________________
F01 <- function(x){
  x %>% 
    dplyr::group_by(samples) %>% 
    dplyr::mutate(S.Area = sum(Area)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(samples, S.Area) %>% 
    unique() %>% 
    dplyr::mutate(M.S.Area = mean(S.Area)) %>% 
    dplyr::select(M.S.Area) %>% 
    slice(1)
}

MSA.B01 <- F01(x = B01)
MSA.B02 <- F01(x = B02)
MSA.B03 <- F01(x = B03)
MSA.B04 <- F01(x = B04)

NP.B01 <- B01 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B01))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B02 <- B02 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B02))) %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100)  

NP.B03 <- B03 %>% 
  dplyr::mutate(N.Area = Area)  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

NP.B04 <- B04 %>% 
  dplyr::mutate(N.Area = (Area*as.numeric(MSA.B03)/as.numeric(MSA.B04)))  %>% 
  dplyr::group_by(samples) %>% 
  dplyr::mutate(P.Area = (N.Area/sum(N.Area))*100) 

DF1 <- NP.B01 %>% 
  dplyr::bind_rows(NP.B02) %>% 
  dplyr::bind_rows(NP.B03) %>%
  dplyr::bind_rows(NP.B04) 
#______________________________________________________________________________________________________

#______________ Create main data frame ________________________________________________________________
Variables <- read_excel("data/Variables.xlsx") %>% 
  dplyr::select(ID, Weight, Time2) %>% 
  dplyr::rename(samples = ID) 

DF <- DF1 %>% 
  dplyr::left_join(., Variables, by = "samples") %>% 
  dplyr::filter(Weight < 7.5) %>% 
  dplyr::filter(Weight > 1) %>% 
  dplyr::select(Name, Weight, samples, N.Area, count, class, species, batch, DB_new, class_new) %>%
  dplyr::mutate(Bond = DB_new) %>% 
  dplyr::mutate(SubClass = class_new) %>% 
  dplyr::select(-DB_new, -class_new) %>% 
  tidyr::separate(samples, into = c("Line", "Cage"), remove = FALSE) %>% 
  dplyr::mutate(Line = replace(Line, Line == "S06", "s06")) %>% 
  dplyr::mutate(Class = class) %>% 
  dplyr::mutate(Batch = batch) %>%
  dplyr::select(-batch) %>% 
  dplyr::rename(Samples = samples) %>% 
  dplyr::mutate(Age = if_else(Batch %in% c("B01", "B02"), "19 Day", "1 Day")) %>% 
  dplyr::mutate(Time = str_extract(Cage, "[a-z]+")) %>% 
  dplyr::mutate(Time = replace_na(Time, "o")) %>% 
  dplyr::mutate(Time = if_else(Time %in% c("o"), "Old", "New")) %>% 
  dplyr::mutate(Time = replace(Time, Line == "s06", "Older")) %>% 
  dplyr::select(-count, -class) %>% 
  tidyr::separate(col = species, into = c("Carbon", NA), ":", remove = FALSE) %>% 
  mutate(Carbon = as.double(Carbon))

#___________________

# Import standards data

df.old <- read_csv("data/Table_ConcAllStandards.csv") %>% # this is data from 1st run
  dplyr::select(Name, Class, area, Conc, Batch) %>%
  dplyr::filter(Class != "Cer" & Class != "FA" & Class != "SM" & Class != "PCp") %>% 
  dplyr::group_by(Name, Class, Conc, Batch) %>% 
  dplyr::summarise(Area = mean(area)) %>% # calculates total area for any class that had more than one standard 
  dplyr::ungroup() %>%
  unique()  %>% 
  dplyr::select(Conc, Name, Class, Area)

# Add manually calculate exp 2 CL areas to DF
DF.CL <- tibble(Conc = c(0.01, 0.10, 1.00, 10.00), 
                Name = rep("CL 72:4",4), Class = rep("CL", 4),
                Area = c(1023.129827,32570.02475,579994.8638,45993575.15))

DF.Stds <- rbind(df.old, DF.CL)

Regression <- DF.Stds %>% split(~Class) %>%
  map(~lm(log10(Area) ~ log10(Conc), data = .x)$coefficient) ### Obtain all the regression coefficient

Regression_coef <- tibble(Class2 = names(Regression) ,
                          bind_rows(lapply(Regression, as.data.frame.list))) %>% ### 
  dplyr::rename(Intercept =  "X.Intercept.",  Slope = "log10.Conc.") %>% 
  mutate(Class2 = if_else(Class2 =="PEp", "PE p", Class2))

#____________create main DF_____________________

DF.main <- DF %>%
  ungroup() %>% 
  mutate(Class2 = if_else(grepl(" e", SubClass), Class, SubClass)) %>% ## All ether lipids use ester lipid standards
  left_join(Regression_coef, by = "Class2") %>%
  mutate(Conc = 10^((log10(N.Area)-Intercept)/Slope)) %>%
  mutate(LogAreaCheck = Intercept + Slope*(log10(Conc)), AreaCheck = 10^LogAreaCheck) ### Just to double check


#______________Make weight plots to identify samples with weight outliers______________________

# Histogram

DF.main %>% 
  dplyr::select(Samples, Line, Weight, Age, Time) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  unique() %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","cbr_New","cbr_Old","syd_New","syd_Old","s06_Older"))) %>%
  Rmisc::summarySE(measurevar = "Weight", groupvars = c("LineTime", "Age")) %>% 
  ggplot(aes(LineTime, Weight, fill = LineTime)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), 
                position = position_dodge(0.8),
                width = 0.5)+
  facet_wrap(~Age)+
  theme_bw()+
  ggtitle("Weight distribution in wild and domesticated Qfly strains" )+
  labs(x= "Qfly strains", y = "Weight")+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold"))+
  theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=1))+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + 
  theme(axis.title = element_text(size = 12))


#_____________label outliers in the weight plot_____________

DF.plot <- DF.main %>% 
  dplyr::select(Samples, Line, Time, Age, Weight) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(-Line,-Time) %>%
  unique() %>% 
  group_by(LineTime, Age) %>% 
  dplyr::mutate(outlier = ifelse(is_outlier(Weight), Weight, as.numeric(NA))) %>% 
  ungroup()

Wt.plot <- DF.plot %>% 
  ggplot(aes(LineTime, Weight, label = Samples))+
  geom_boxplot(outlier.color = "Red") +
  ggrepel::geom_text_repel(data = (DF.plot %>% drop_na() %>% unique()),
                           aes(label = Samples),
                           na.rm = TRUE,
                           color = "Red") +
  stat_summary(fun=mean, colour="red", geom="point",
               shape=18, size=3) +
  facet_grid(~Age) +
  theme_bw()

#__________Calculate total lipids and percent of fly weight due to lipid____________

# Here we calc total lipids in mglipid/fly and SE for each strain and we calculate
# percent of fly weight due to lipid and SE for each strain


#_______ First remake DF.new by removing weight outliers from DF.main___________

DF.new <- DF.main %>% 
  dplyr::filter(Samples != "ct_o225" & Samples != "syd_o169") %>% # remove weight outlier for Day 1
  dplyr::filter(Samples != "ct_o42" & Samples != "cbr_n43") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "cbr_n62" & Samples != "cbr_n86") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "syd_o22" & Samples != "syd_o76") # remove weight outlier for Day 19 


TL.WT <- DF.new %>% 
  ungroup %>%
  group_by(Samples) %>%
  #mutate(LipidWT = (Conc*Weight*50)) %>% #LipidWT in ng/fly
  mutate(LipidWT = (Conc*Weight*50)*1e-6) %>%  #Convert from ng to mg so LipidWt in mg/fly
  ungroup() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Total.Lipid = sum(LipidWT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Line, Age, Time, Weight, Total.Lipid) %>% 
  unique() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Percent = (as.numeric(Total.Lipid)/as.numeric(Weight)*100)) %>% 
  ungroup() %>% 
  group_by(Line, Time, Age) %>% 
  dplyr::summarise(TL=mean(Total.Lipid),TLSE = sd(Total.Lipid)/sqrt(length(Total.Lipid)),
                   PL=mean(Percent),PLSE = sd(Percent)/sqrt(length(Percent)),
                   WT=mean(Weight),WTSE = sd(Weight)/sqrt(length(Weight))) %>%  
  ungroup() %>% 
  dplyr::mutate(TL = format(round(.$TL, 3), nsmall = 3)) %>% 
  dplyr::mutate(TLSE = format(round(.$TLSE, 3), nsmall = 3)) %>% 
  dplyr::mutate(PL = format(round(.$PL, 3), nsmall = 3)) %>% 
  dplyr::mutate(PLSE = format(round(.$PLSE, 3), nsmall = 3)) %>% 
  dplyr::mutate(WT = format(round(.$WT, 3), nsmall = 3)) %>% 
  dplyr::mutate(WTSE = format(round(.$WTSE, 3), nsmall = 3)) %>% 
  dplyr::mutate(Age = replace(Age, Age == "1 Day", "D1")) %>% 
  dplyr::mutate(Age = replace(Age, Age == "19 Day", "D19")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "s06", "SD")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "cbr", "CN")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "ct", "CT")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "syd", "SD"))  %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::mutate(TL = as.numeric(TL), TLSE = as.numeric(TLSE), 
                PL = as.numeric(PL), PLSE = as.numeric(PLSE),
                WT = as.numeric(WT), WTSE = as.numeric(WTSE)) 

#_____________make a correlation plot of weight vs total lipids_________

pd = position_dodge(0.0005)

WT.TL_D1 <- TL.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  ggplot(aes(WT, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TLSE, ymax=TL+TLSE,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Total lipids (mg/fly)"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  stat_cor(aes(), color = "black", label.y.npc="top", label.x.npc = "left", geom = "label", size = 2.5) 
  
pd = position_dodge(0.4)

WT.TL_D19 <- TL.WT %>% 
  dplyr::filter(Age == "D19") %>% 
  ggplot(aes(WT, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TLSE, ymax=TL+TLSE,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Total lipids (mg/fly)"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  stat_cor(aes(), color = "black", label.y.npc="top", label.x.npc = "left", geom = "label", size = 2.5)

pd = position_dodge(0.0005)

Legend1 <- TL.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  ggplot(aes(WT, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TLSE, ymax=TL+TLSE,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size =8))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Total lipids (mg/fly)"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size = 8))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0))

#Arrange and save the correlation plots for DB

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(Legend1)

Figure1 <- ggarrange(WT.TL_D1, WT.TL_D19, ncol = 2)

Figure1_annotated <- annotate_figure(Figure1, bottom = text_grob("Dry weight (mg)", hjust = 0.5, face = "bold", size = 8),
                                    left = text_grob("Total lipids (mg/fly)", rot = 90, face = "bold", size = 8))

Figure1_arranged <- grid.arrange(Figure1_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                               widths = c(2.7, 2.7), heights = c(2.5, 0.3))


#ggsave(plot = Figure1_arranged, width = 7.0, height = 3.5, units = "in", dpi = 300,filename = "Figures/Weight_vs_total_lipid_titre_correlation_plot.jpg")              


#_________________________________________________________________________________________


#_____________Plot the percent of fly weight due to lipid_________________________

PL.WT <- TL.WT %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "CT_New", "CTnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "CT_Old", "CTold")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "SD_New", "SDnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "SD_Old", "SDold")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "SD_Older", "SDolder")) %>% 
  dplyr::mutate(LineTime= replace(LineTime, LineTime == "CN_New", "CNnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "CN_Old", "CNold")) %>% 
  dplyr::mutate(Strain = LineTime) %>% 
  dplyr::select(11,2,3,4,7,8)

pd = position_dodge(0.4)

PL.WT_D1 <- PL.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Strain, PL))+ 
  geom_point(mapping = aes(colour = Strain, shape = Strain), size = 4, position = pd) +
  theme_bw() +
  #geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=PL-PLSE, ymax=PL+PLSE,colour = Strain), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Strain")) +
  ylab(paste0("Fly weight due to lipids (%)"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

PL.WT_D19 <- PL.WT %>% 
  dplyr::filter(Age == "D19") %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Strain, PL))+ 
  geom_point(mapping = aes(colour = Strain, shape  = Strain), size = 4, position = pd) +
  theme_bw() +
  #geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=PL-PLSE, ymax=PL+PLSE,colour = Strain), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Strain")) +
  ylab(paste0("Fly weight due to lipids (%)"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

Figure2 <- ggarrange(PL.WT_D1, PL.WT_D19, ncol = 2)

Figure2_annotated <- annotate_figure(Figure2, bottom = text_grob("Strain", hjust = 0.5, face = "bold", size = 12),
                                    left = text_grob("Lipid weight (%)", rot = 90, face = "bold", size = 12))

#ggsave(plot = Figure2_arranged, width = 7.5, height = 3.5, units = "in", dpi = 300,filename = "Figures/Percent of fly weight due to lipids.jpg")              

#Save Figure2_annotated as FigureS1
ggsave(plot = Figure2_annotated, width = 8.0, height = 4.0, units = "in", dpi = 300, filename = "Figures/FigureS1.jpg")              


#_______________________________________________________________________________________________


#____________________Plot the percent of lipid weight vs fly weight_______________________________

PLWT_WT <- TL.WT %>% 
  dplyr::select(-TL, -TLSE)

pd = position_dodge(0.4)

arrows1 <- PLWT_WT  %>%
  dplyr::filter(Age == "D1") %>% 
  dplyr::select(Line, Time, WT, PL) %>%
  pivot_wider(names_from = Time, values_from = c(WT, PL)) %>%
  dplyr::rename(x_start =  WT_New, x_end = WT_Old, y_start = PL_New, y_end = PL_Old)

Extra1 <- data.frame(Line = "SD", x_start = arrows1$x_end[3], x_end = arrows1$WT_Older[3],
                     y_start = arrows1$y_end[3], y_end = arrows1$PL_Older[3])

Final_arrows1 <- rbind(arrows1 %>% dplyr::select(-WT_Older,-PL_Older), Extra1)


PLWT_WT_D1 <- PLWT_WT %>% 
  dplyr::filter(Age == "D1") %>% 
  ggplot(aes(WT, PL))+ 
  geom_point(mapping = aes(shape = Time,colour = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=PL-PLSE, ymax=PL+PLSE,colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Lipid weight (%)"))+
  ggarchery::geom_arrowsegment(data=Final_arrows1, ### Requires a separate data frame to specify start and end of lines
                               mapping=aes(x = x_start, xend = x_end, y = y_start, 
                                           yend = y_end, colour = Line),
                               arrows = arrow(angle = 30, type = "open", length = unit(0.12, "inches")), # adjust length to adjust size of arrow head 
                               arrow_positions = 0.5, # location of arrow along the line
                               size = 0.5) +
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("CT" = "red", "CN" = "green", "SD" = "blue"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
 

arrows2 <- PLWT_WT  %>%
  dplyr::filter(Age == "D19") %>% 
  dplyr::select(Line, Time, WT, PL) %>%
  pivot_wider(names_from = Time, values_from = c(WT, PL)) %>%
  dplyr::rename(x_start =  WT_New, x_end = WT_Old, y_start = PL_New, y_end = PL_Old)

Extra2 <- data.frame(Line = "SD", x_start = arrows2$x_end[3], x_end = arrows2$WT_Older[3],
                     y_start = arrows2$y_end[3], y_end = arrows2$PL_Older[3])

Final_arrows2 <- rbind(arrows2 %>% dplyr::select(-WT_Older,-PL_Older), Extra2)


PLWT_WT_D19 <- PLWT_WT %>% 
  dplyr::filter(Age == "D19") %>% 
  ggplot(aes(WT, PL))+ 
  geom_point(mapping = aes(shape = Time,colour = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=PL-PLSE, ymax=PL+PLSE,colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Lipid weight (%)"))+
  ggarchery::geom_arrowsegment(data=Final_arrows2, ### Requires a separate data frame to specify start and end of lines
                               mapping=aes(x = x_start, xend = x_end, y = y_start, 
                                           yend = y_end, colour = Line),
                               arrows = arrow(angle = 30, type = "open", length = unit(0.12, "inches")), # adjust length to adjust size of arrow head 
                               arrow_positions = 0.5, # location of arrow along the line
                               size = 0.5) +
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("CT" = "red", "CN" = "green", "SD" = "blue"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 


#Arrange and save the correlation plots for DB

Legend2 <- PLWT_WT %>% 
  dplyr::filter(Age == "D1") %>% 
  ggplot(aes(WT, PL))+ 
  geom_point(mapping = aes(shape = Time,colour = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = WT-WTSE, xmax = WT+WTSE, colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=PL-PLSE, ymax=PL+PLSE,colour = Line), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight (mg)")) +
  ylab(paste0("Lipid weight (%)"))+
  ggarchery::geom_arrowsegment(data=Final_arrows1, ### Requires a separate data frame to specify start and end of lines
                               mapping=aes(x = x_start, xend = x_end, y = y_start, 
                                           yend = y_end, colour = Line),
                               arrows = arrow(angle = 30, type = "open", length = unit(0.12, "inches")), # adjust length to adjust size of arrow head 
                               arrow_positions = 0.5, # location of arrow along the line
                               size = 0.5) +
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size = 8))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("CT" = "red", "CN" = "green", "SD" = "blue"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 


legend2 <- get_legend(Legend2)

Figure3 <- ggarrange(PLWT_WT_D1, PLWT_WT_D19, ncol = 2)

Figure3_annotated <- annotate_figure(Figure3, bottom = text_grob("Dry weight (mg)", hjust = 0.5, face = "bold", size = 8),
                                    left = text_grob("Lipid weight (%)", rot = 90, face = "bold", size = 8))

Figure3_arranged <- grid.arrange(Figure3_annotated, legend2, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                widths = c(2.7, 2.7), heights = c(2.5, 0.3))


#ggsave(plot = Figure3_arranged, width = 8.0, height = 6.0, units = "in", dpi = 300,filename = "Figures/Percent of fly weight due to lipids vs dry fly weight.jpg")              

# Arrange and figures 1,2,3 

#All_Figures <- ggarrange(Figure1_arranged, Figure2_annotated, Figure3_arranged, ncol = 1, nrow = 3,
#                      heights = c(1.6, 1, 1.6))

#________________________________________END_________________________________________

