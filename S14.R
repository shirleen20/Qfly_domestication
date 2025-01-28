#___________________________________________25/01/2025________________________________________

# Here we generate plots showing significant main and/or interaction effect of strain for 
# Ether lipids category for variables abundance, CC and DB
# Strain-specific means and 95% confidence intervals are shown. 
# Where the interaction with weight was not sig, single plot is given showing estimates for 
# all the strains, where the interaction was significant three plots are given, one for 
# lightest 10% of flies one for mean weight, and one for heaviest 10% of flies.
# We use this analysis to generate Figures 4 and 5, Tables 2 and 3.

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("broom")
library("gtools")
library("emmeans")
library("multcomp")
library("multcompView")


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
  mutate(Carbon = as.double(Carbon))%>% 
  dplyr::filter(Samples != "ct_o225" & Samples != "syd_o169") %>% # remove weight outlier for Day 1
  dplyr::filter(Samples != "ct_o42" & Samples != "cbr_n43") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "cbr_n62" & Samples != "cbr_n86") %>% # remove weight outlier for Day 19 
  dplyr::filter(Samples != "syd_o22" & Samples != "syd_o76") # remove weight outlier for Day 19 

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


#___________________________Variable : Abundance each category _________________________

# complete fxn fills in zeros for samples that have missing lipids since the numbers we analyse for % abundances 
# are always using all replicates even if the lipid in question was not present in all of them 

Complete_DF <- DF.main %>% ungroup() %>% 
  complete(SubClass, nesting(Samples, Age, Line, Time, Weight), fill = list(Conc = 0)) %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Conc, Weight)


# use DF to calculate % abundance & SE for each of the 4 lipid categories for s06 only for Table 1 of lipid MS1

# add another column to DF that has the following 4 lipid categories
PLs <- c("CL","PC", "PE", "PG", "PI", "PS","LPC")
NLs <- c("DG", "TG")
ELPLs <- c("PC e", "PE e", "PE p", "PS e")
ELNLs <- c("DG e", "TG e")


#_________________Variable: Abundance each category_______________________________________________

# NOTE: where emtrends was sig we will extract and plot weight emmeans using global fly weight, top 10% and bottom 10% fly weight.

DF_glm01_1 <- Complete_DF %>% 
  ungroup %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% PLs,"PLs", 
                              ifelse(SubClass %in% ELPLs,"ELPLs", ifelse(SubClass %in% ELNLs,"ELNLs", "NLs")))) %>%  
  dplyr::group_by(Samples,Type) %>% 
  dplyr::mutate(Total.Conc.by.Type = sum(Conc)) %>% # by samples, age and type
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>%
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Weight, Age, Line, Time,  Type, Total.Conc.by.Type, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.Type/Total.Conc.by.SA)*100) %>% 
  unique() %>%
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF_glm01 %>% filter(Samples == "cbr_n04"))$Percentage)
  dplyr::select(-Total.Conc.by.Type, -Total.Conc.by.SA) %>%
  dplyr::filter(Age=="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Percentage), SD=sd(Percentage)) %>% 
  mutate(Percentage.t = (Percentage-MEAN)/SD) %>% 
  ungroup() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(Type, Age, Samples, LineTime, Percentage.t, Weight) 

# Use a mean weight for all 7 stains i.e we include S06

Weights_Category_1 <- DF.main  %>%
  dplyr::filter(Age =="1 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 


#_____Day 1 ELNL_____________

Abundance.global.1Day.ELNL <- DF_glm01_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Percentage.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(id == "ELNLs") %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

pd = position_dodge(0.4)    ### How much to jitter the points on the plot

Abundance.global.1Dayplot.ELNL <- Abundance.global.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.5, 2), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Abundance Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.2, 1.2, 1.2), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 

#ggsave(plot = Abundance.global.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/AB_ELNLglobalD1.jpg")                

#_____Day 1 ELPL_____________

Abundance.global.1Day.ELPL <- DF_glm01_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Percentage.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(id == "ELPLs") %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

Abundance.global.1Dayplot.ELPL <- Abundance.global.1Day.ELPL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1, 2), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Abundance Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.2, 1.2, 1.2), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 

#ggsave(plot = Abundance.global.1Dayplot.ELPL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure5/AB_ELPLglobalD1.jpg")                

#_____________________________Day 19 _______________________________

DF_glm01_19 <- Complete_DF %>% 
  ungroup %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% PLs,"PLs", 
                              ifelse(SubClass %in% ELPLs,"ELPLs", ifelse(SubClass %in% ELNLs,"ELNLs", "NLs")))) %>%  
  dplyr::group_by(Samples,Type) %>% 
  dplyr::mutate(Total.Conc.by.Type = sum(Conc)) %>% # by samples, age and type
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>%
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Weight, Age, Line, Time,  Type, Total.Conc.by.Type, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.Type/Total.Conc.by.SA)*100) %>% 
  unique() %>%
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF_glm01 %>% filter(Samples == "cbr_n04"))$Percentage)
  dplyr::select(-Total.Conc.by.Type, -Total.Conc.by.SA) %>%
  dplyr::filter(Age=="19 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Percentage), SD=sd(Percentage)) %>% 
  mutate(Percentage.t = (Percentage-MEAN)/SD) %>% 
  ungroup() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(Type, Age, Samples, LineTime, Percentage.t, Weight) 

# Use a mean weight for all 7 stains i.e we include S06
Weights_Category_19 <- DF.main %>%
  dplyr::filter(Age =="19 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 


#________Day 19 ELPL_______________

Abundance.global.19Day.ELPL <- DF_glm01_19 %>% 
  split(.$Type) %>% 
  map(~ lm(Percentage.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 4.59)) %>%  
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(id == "ELPLs") %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

Abundance.global.19Dayplot.ELPL <- Abundance.global.19Day.ELPL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.5, 2), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Abundance Day 19")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.4, 1.4, 1.4), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 


#ggsave(plot = Abundance.global.19Dayplot.ELPL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure5/AB_ELPLglobalD19.jpg")                

#________________________________________________________________________________________________________________


#___________________Variable: Mean carbon chain length in lipid category_________________________


# Create dataframe to calculate mean CC and mean DB in the akyl and acyl chains of ELNLs and ELPLs category

#____________

DFnew <- DF.main %>% dplyr::filter(SubClass != "PS e" & SubClass != "DG e") #remove subclasses PSe and DGe from DF because these two subclasses have all lipid species that do not have acyl chain identified

F1 <- function(X, Y){
  Subclass <- X %>% 
    filter(SubClass == Y)
  length1 <- length(Subclass)  # The starting number of columns
  Names <- colnames(Subclass) # Save column names
  Subclass_SC <- gsub(".*\\((.+)\\).*","\\1",Subclass$Name) %>%  ### keep only string within "( )"
    str_remove_all("[ep]") %>%           
    strsplit("_") ### split by "_"
  SCNumber <- Subclass_SC %>% map(length) %>% unlist ### find the number of side chains in each lipid for each individual (want to remove those with no side chains identified)
  RemoveSpecies <- which(SCNumber < max(SCNumber)) ### vector of indices of those without side chains identified. They will have length less than the maximum number of side chains
  Subclass_SC <- do.call(rbind, Subclass_SC) %>% as.data.frame() ### 
  Subclass <- cbind(Subclass, Subclass_SC) ### Join the new columns to the old data frame
  length2 <- length(Subclass) ### New number of columns
  colnames(Subclass) <- c(Names, paste("SC", 1:(length2-length1), sep = ""))  ### Rename new colunms with SC1, SC2 ... etc length2-length1 is the total number of new columns
  if (length(RemoveSpecies) == 0){
    return(Subclass)
  }
  else{
    return(Subclass[-RemoveSpecies,]) ### Output new data frame after removing rows with unindentified side chains
  }
}

F2 <- function(X, Subclass = "SubClass"){
  Subclass_names <- unique(X[[Subclass]]) ### All the unique Subclasses
  Out_list <- setNames(replicate(length(Subclass_names), data.frame()), Subclass_names) ### create a list of empty data frames with names Subclass_names with length of the number of subclasses 
  # if using vector(mode = "list", length = length(Subclass_names)) ... required renaming objects in Out_list, see next line
  # names(Out_list) <- Subclass_names ## rename each object within the list with the subclass names
  for (i in Subclass_names){ ### For each separate subclass
    Out_list[[i]] <- F1(X, i) ### in each object of that ith subclass, run F1 function on data frame X with subclass i
  }
  return(Out_list)
}

LIST1 <- F2(X = DFnew)

# add another column to LIST called type that has lipid categories
NLs <- c("TG", "DG")
ELNLs <- c("TG e", "DG e")
PLs <- c("CL","PC", "PE", "PG", "PI", "PS","LPC")
ELPLs <- c("PE p", "PE e", "PC e", "PS e")

Composition <- bind_rows(LIST1, .id = "Class") %>% 
  dplyr::select(Name, SubClass, Samples, Line, Time, Age, Weight, Conc, SC1, SC2, SC3, SC4) %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% PLs, "PLs", 
                              ifelse(SubClass %in% ELPLs, "ELPLs",
                                     ifelse(SubClass %in% ELNLs, "ELNLs",
                                            "NLs"))))


#___________________________Variable: Mean Acyl/Alkyl carbon chain length analysis_________________________

CCAlkyl_1 <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(AlkylCC.Wt = as.numeric(AlkylCC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, AlkylCC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylCC.Wt'), names_to='AlkylChain', values_to='AlkylChainCC') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.CarbonAlkyl) %>%
  unique() %>%
  dplyr::filter(Age =="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.CarbonAlkyl), SD=sd(Mean.CarbonAlkyl)) %>% 
  mutate(Mean.CarbonAlkyl.t = (Mean.CarbonAlkyl-MEAN)/SD) %>% 
  ungroup() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

# Use a mean weight for all 7 stains 
Weights_Category_1 <- DF.main  %>%
  dplyr::filter(Age =="1 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 


#___Day 1 AlkylCC ELNLs______

AlkylCC.global.1Day.ELNL <- CCAlkyl_1 %>% 
  split(.$Type) %>% 
  map(~ lm(Mean.CarbonAlkyl.t ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~LineTime, at = list(Weight = 1.96)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylCC.global.1Dayplot.ELNL <- AlkylCC.global.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Alkyl chain carbon content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.8, 0.8, 0.8), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylCC.global.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AlkylchainELNLglobalD1.jpg")

#___Day 1 AcylCC ELNLs______

CCAcyl_1 <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(Acyl1CC.Wt = as.numeric(Acyl1CC)*Conc) %>% 
  dplyr::mutate(Acyl2CC.Wt = as.numeric(Acyl2CC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Acyl1CC.Wt, Acyl2CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1CC.Wt', 'Acyl2CC.Wt'), names_to='AcylChain', values_to='AcylChainCC') %>% 
  dplyr::filter(AcylChainCC != "NA") %>% # to remove SC3 from ELPLs
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AcylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.CarbonAcyl) %>%
  unique() %>% 
  dplyr::filter(Age =="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.CarbonAcyl), SD=sd(Mean.CarbonAcyl)) %>% 
  mutate(Mean.CarbonAcyl.t = (Mean.CarbonAcyl-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

#___
AcylCC.global.1Day.ELNL <- CCAcyl_1 %>% 
  split(.$Type) %>%   
  map(~ lm(Mean.CarbonAcyl.t ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~LineTime, at = list(Weight = 1.96)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AcylCC.global.1Dayplot.ELNL <- AcylCC.global.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Acyl chain carbon content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.8, 0.8, 0.8), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AcylCC.global.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AcylchainELNLglobalD1.jpg")


#___Day 1 AcylCC ELPLs______

AcylCC.global.1Day.ELPL <- CCAcyl_1 %>% 
  split(.$Type) %>%   
  map(~ lm(Mean.CarbonAcyl.t ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~LineTime, at = list(Weight = 1.96)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELPLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AcylCC.global.1Dayplot.ELPL <- AcylCC.global.1Day.ELPL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Acyl chain carbon content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.8, 0.8, 0.8), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))


#ggsave(plot = AcylCC.global.1Dayplot.ELPL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure5/CC_AcylchainELPLglobalD1.jpg")

#________________________Day 19_______________________________

#_________Day 19 AlkylCC___________

CCAlkyl_19 <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "ALkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(AlkylCC.Wt = as.numeric(AlkylCC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, AlkylCC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylCC.Wt'), names_to='AlkylChain', values_to='AlkylChainCC') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.CarbonAlkyl) %>%
  unique() %>% 
  dplyr::filter(Age =="19 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.CarbonAlkyl), SD=sd(Mean.CarbonAlkyl)) %>% 
  mutate(Mean.CarbonAlkyl.t = (Mean.CarbonAlkyl-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

# Use a mean weight for all 7 stains 
Weights_Category_19 <- DF.main  %>%
  dplyr::filter(Age =="19 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 


#______Day 19 AlkylCC ELNLs______

AlkylCC.global.19Day.ELNL <- CCAlkyl_19 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.CarbonAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 4.59)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylCC.global.19Dayplot.ELNL <- AlkylCC.global.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-5.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.5, 1.5, 1.5), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylCC.global.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AlkylchainELNLglobalD19.jpg")

#______

AlkylCC.bottom.19Day.ELNL <- CCAlkyl_19 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.CarbonAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 3.20)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylCC.bottom.19Dayplot.ELNL <- AlkylCC.bottom.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-5.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.82, 1.82, 1.82), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylCC.bottom.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AlkylchainELNLbottomD19.jpg")

#______

AlkylCC.top.19Day.ELNL <- CCAlkyl_19 %>% 
  split(.$Type) %>%    
  map(~ lm(Mean.CarbonAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 6.14)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylCC.top.19Dayplot.ELNL <- AlkylCC.top.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-5.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(2.3, 2.3, 2.3), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylCC.top.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AlkylchainELNLtopD19.jpg")


#_________Day 19 AcylCC___________

CCAcyl_19 <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(Acyl1CC.Wt = as.numeric(Acyl1CC)*Conc) %>% 
  dplyr::mutate(Acyl2CC.Wt = as.numeric(Acyl2CC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Acyl1CC.Wt, Acyl2CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1CC.Wt', 'Acyl2CC.Wt'), names_to='AcylChain', values_to='AcylChainCC') %>% 
  dplyr::filter(AcylChainCC != "NA") %>% # to remove SC3 from ELPLs
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AcylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.CarbonAcyl) %>%
  unique() %>% 
  dplyr::filter(Age =="19 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.CarbonAcyl), SD=sd(Mean.CarbonAcyl)) %>% 
  mutate(Mean.CarbonAcyl.t = (Mean.CarbonAcyl-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

#___
AcylCC.global.19Day.ELNL <- CCAcyl_19 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.CarbonAcyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 4.59)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AcylCC.global.19Dayplot.ELNL <- AcylCC.global.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-3.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.5, 1.5, 1.5), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AcylCC.global.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AcylchainELNLglobalD19.jpg")

#______

AcylCC.bottom.19Day.ELNL <- CCAcyl_19 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.CarbonAcyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 3.20)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AcylCC.bottom.19Dayplot.ELNL <- AcylCC.bottom.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-3.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.8, 1.8, 1.8), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AcylCC.bottom.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AcylchainELNLbottomD19.jpg")

#______

AcylCC.top.19Day.ELNL <- CCAcyl_19 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.CarbonAcyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 6.14)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AcylCC.top.19Dayplot.ELNL <- AcylCC.top.19Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-5, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(2.45, 2.45, 2.45), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AcylCC.top.19Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/CC_AcylchainELNLtopD19.jpg")

#__________________________________________________________________________________________________________


#___________________Variable: Double bond content in ether lipid alkyl/acyl chains__________________________

DBAlkyl_1 <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(AlkylDB.Wt = as.numeric(AlkylDB)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, AlkylDB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylDB.Wt'), names_to='AlkylChain', values_to='AlkylChainDB') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.BondAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.BondAlkyl) %>%
  unique() %>% 
  dplyr::filter(Age =="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.BondAlkyl), SD=sd(Mean.BondAlkyl)) %>% 
  mutate(Mean.BondAlkyl.t = (Mean.BondAlkyl-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

# Use a mean weight for all 7 stains 
Weights_Category_1 <- DF.main  %>%
  dplyr::filter(Age =="1 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 


#___Day 1 AlkylDB ELNLs______

AlkylDB.global.1Day.ELNL <- DBAlkyl_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.BondAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylDB.global.1Dayplot.ELNL <- AlkylDB.global.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-3.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(2, 2, 2), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylDB.global.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/DB_AlkylchainELNLglobalD1.jpg")

#_____

AlkylDB.bottom.1Day.ELNL <- DBAlkyl_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.BondAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 1.39)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylDB.bottom.1Dayplot.ELNL <- AlkylDB.bottom.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-3.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(2, 2, 2), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylDB.bottom.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/DB_AlkylchainELNLbottomD1.jpg")

#______

AlkylDB.top.1Day.ELNL <- DBAlkyl_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.BondAlkyl.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>%
        emmeans(~interaction(LineTime), at = list(Weight = 2.54)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELNLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylDB.top.1Dayplot.ELNL <- AlkylDB.top.1Day.ELNL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-3.0, 3.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(hjust = 0)) +
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1.3, 1.3, 1.3), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylDB.top.1Dayplot.ELNL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure4/DB_AlkylchainELNLtopD1.jpg")

#______________________________________


#___Day 1 AlkylDB ELPLs______

AlkylDB.global.1Day.ELPL <- DBAlkyl_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.BondAlkyl.t ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~LineTime, at = list(Weight = 1.96)) %>%
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "ELPLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

AlkylDB.global.1Dayplot.ELPL <- AlkylDB.global.1Day.ELPL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Alkyl chain double bond content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.5, 0.5, 0.5), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))

#ggsave(plot = AlkylDB.global.1Dayplot.ELPL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure5/DB_AlkylchainELPLglobalD1.jpg")

#______________________________END_____________________________________________

