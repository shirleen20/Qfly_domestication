#______________________________________24/01/2025_______________________________

# Here we calculate fly weight weight differences and plot dry weight vs total lipid titre for 
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

#___________ save data for the weight of each sample____________

# Find mean weight and SE for each of the seven strain

WeightsSE <- DF.main %>% 
  ungroup() %>%
  dplyr::filter(Samples != "ct_o225" & Samples != "syd_o169") %>% # remove weight outlier for Day 1
  dplyr::filter(Samples != "ct_o42" & Samples != "cbr_n43") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "cbr_n62" & Samples != "cbr_n86") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "syd_o22" & Samples != "syd_o76") %>% # remove weight outlier for Day 19
  dplyr::select(Samples, Line, Age, Time, Weight) %>% 
  unique() %>% 
  group_by(Samples) %>% 
  summarySE(groupvars = c("Line", "Time", "Age"), measurevar = "Weight") %>% 
  ungroup() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(LineTime, Age, Weight, se) %>% 
  dplyr::mutate(Wt = format(round(.$Weight, 2), nsmall = 2)) %>% 
  dplyr::mutate(WT.se = format(round(.$se, 2), nsmall = 2)) %>% 
  #tidyr::unite("WtSE", Wt:se, sep = "±") #%>% 
  #tidyr::pivot_wider(id_cols = c(Age), names_from = LineTime, values_from = WtSE)
  dplyr::select(1,2,5,6)


#__________Calculate Total lipid titre in lipid/fly and SE for each strain_____________


TLSE <- DF.main %>% 
  ungroup %>%
  group_by(Samples) %>%
  #mutate(LipidWT = (Conc*Weight*50)) %>% 
  mutate(LipidWT = (Conc*Weight*50)*1e-6) %>%  #Convert from ng to mg
  ungroup() %>% 
  dplyr::group_by(Age,Samples) %>% 
  dplyr::mutate(Total.Lipid = sum(LipidWT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Line, Age, Time, Weight, Total.Lipid) %>% 
  unique() %>% 
  group_by(Samples) %>% 
  summarySE(groupvars = c("Line", "Time", "Age"), measurevar = "Total.Lipid") %>% 
  ungroup() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>%
  dplyr::select(LineTime, Age, Total.Lipid, se) %>% 
  dplyr::mutate(TL = format(round(.$Total.Lipid, 2), nsmall = 2)) %>% 
  dplyr::mutate(TL.se = format(round(.$se, 2), nsmall = 2)) %>% 
  #tidyr::unite("TLSE", Total.Conc:se, sep = "±")
  dplyr::select(1,2,5,6)

# Create a dataframe with total lipid + SE and Weight + SE 

TL.WT <- cbind(WeightsSE, TLSE) %>% 
  dplyr::select(-5,-6) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = TRUE) %>% 
  dplyr::mutate(Wt = as.numeric(Wt), WT.se = as.numeric(WT.se), TL = as.numeric(TL), TL.se = as.numeric(TL.se)) %>% 
  dplyr::mutate(Age = replace(Age, Age == "1 Day", "D1")) %>% 
  dplyr::mutate(Age = replace(Age, Age == "19 Day", "D19")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "s06", "SD")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "cbr", "CN")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "ct", "CT")) %>% 
  dplyr::mutate(Line = replace(Line, Line == "syd", "SD")) 

#_____________make a correlation plot of weight vs total lipid titre_________

pd = position_dodge(0.4)

Wt.TL.plotD1 <- TL.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  #ggplot(aes(Wt, log10(TL)))+ 
  ggplot(aes(Wt, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight")) +
  ylab(paste0("Total lipids"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  ylim(expand = c(0, 0.50))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

Wt.TL.plotD19 <- TL.WT %>% 
  dplyr::filter(Age == "D19") %>% 
  #ggplot(aes(Wt, log10(TL)))+ 
  ggplot(aes(Wt, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight")) +
  ylab(paste0("Total lipids"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

Legend <- TL.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  #ggplot(aes(Wt, log10(TL)))+ 
  ggplot(aes(Wt, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.01, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.01, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size =10))+
  xlab(paste0("Dry weight")) +
  ylab(paste0("Total lipid titre"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size = 10))+
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

legend <- get_legend(Legend)

Figure <- ggarrange(Wt.TL.plotD1, Wt.TL.plotD19, ncol = 2)

Figure_annotated <- annotate_figure(Figure, bottom = text_grob("Dry weight", hjust = 0.5, face = "bold", size = 12),
                                    left = text_grob("Total lipid titre", rot = 90, face = "bold", size = 12))

Figure_arranged <- grid.arrange(Figure_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                widths = c(2.7, 2.7), heights = c(2.5, 0.3))


ggsave(plot = Figure_arranged, width = 7.0, height = 3.5, units = "in", dpi = 300,filename = "Figures/Weight_vs_TotalLipids_Correlation_plot.jpg")              


#_________________________________________________________________________________________


#_____________Plot the percent of fly weight due to lipid_________________________

TL.P.WT <- cbind(WeightsSE, TLSE) %>% 
  dplyr::select(1,2,3,7) %>% 
  dplyr::mutate(Wt = as.numeric(Wt), TL = as.numeric(TL)) %>% 
  dplyr::mutate(Percent = (as.numeric(TL)/as.numeric(Wt))) %>% 
  #tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = TRUE) %>% 
  dplyr::mutate(Age = replace(Age, Age == "1 Day", "D1")) %>% 
  dplyr::mutate(Age = replace(Age, Age == "19 Day", "D19")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "ct_New", "CTnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "ct_Old", "CTold")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "syd_New", "SDnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "syd_Old", "SDold")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "s06_Older", "SDolder")) %>% 
  dplyr::mutate(LineTime= replace(LineTime, LineTime == "cbr_New", "CNnew")) %>% 
  dplyr::mutate(LineTime = replace(LineTime, LineTime == "cbr_Old", "CNold")) %>% 
  dplyr::mutate(Strain = LineTime) %>% 
  dplyr::select(6,2,3,4,5)


TL.P.WT.plotD1 <- TL.P.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Strain, Percent, colour = Strain))+ 
  #geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  #geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  #geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight")) +
  ylab(paste0("% fly weight due to lipids"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  #scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  #scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) #+
  #stat_cor(aes(), color = "black", label.y.npc="top", label.x.npc = "left", geom = "label", size = 2.5)

TL.P.WT.plotD19 <- TL.P.WT %>% 
  dplyr::filter(Age == "D19") %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Strain, Percent, colour = Strain))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.0001, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 10))+
  xlab(paste0("Fly weight")) +
  ylab(paste0("Total lipids"))+
  theme(legend.position="none")+
  theme(legend.text=element_text(size = 10))+
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 12))+
  scale_colour_manual(values = c("New" = "deepskyblue", "Old" = "blue", "Older" = "red"))+
  scale_shape_manual(values = c("CT" = 2, "CN" = 5, "SD" = 0)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())+
  stat_cor(aes(), color = "black", label.y.npc="top", label.x.npc = "centre", geom = "label", size = 2.5)

Legendp <- TL.P.WT %>% 
  dplyr::filter(Age == "D1") %>% 
  #ggplot(aes(Wt, log10(TL)))+ 
  ggplot(aes(Wt, TL))+ 
  geom_point(mapping = aes(colour = Time, shape = Line), size = 4, position = pd) +
  theme_bw() +
  geom_errorbar(aes(xmin = Wt-WT.se, xmax = Wt+WT.se, colour = Time), width =  0.01, size  =  0.5, position = pd) +
  geom_errorbar(aes(ymin=TL-TL.se, ymax=TL+TL.se,colour = Time), width =  0.01, size  =  0.5, position = pd) +
  facet_grid(~ Age)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size =10))+
  xlab(paste0("Dry weight")) +
  ylab(paste0("Total lipid titre"))+
  theme(legend.position="bottom")+
  theme(legend.text=element_text(size = 10))+
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

legend <- get_legend(Legend)

Figure_p <- ggarrange(Wt.TL.plotD1, Wt.TL.plotD19, ncol = 2)

Figurep_annotated <- annotate_figure(Figure_p, bottom = text_grob("Dry weight", hjust = 0.5, face = "bold", size = 12),
                                    left = text_grob("Total lipid titre", rot = 90, face = "bold", size = 12))

Figurep_arranged <- grid.arrange(Figurep_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                widths = c(2.7, 2.7), heights = c(2.5, 0.3))

ggsave(plot = Figurep_arranged, width = 7.0, height = 3.5, units = "in", dpi = 300,filename = "Figures/Percent of l=fly weight due to lipids.jpg")              


#________________________________________END_________________________________________
