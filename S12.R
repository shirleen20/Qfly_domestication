#___________________________________________25/01/2025________________________________________

# Here we generate plots showing significant main and/or interaction effect of strain for 
# Phospholipids category for variables abundance, CC and DB
# Strain-specific means and 95% confidence intervals are shown. 
# Where the interaction with weight was not sig, single plot is given showing estimates for 
# all the strains, where the interaction was significant three plots are given, one for 
# lightest 10% of flies one for mean weight, and one for heaviest 10% of flies.
# We use this analysis to generate Figure 3, Tables 2 and 3.

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


#___________________________Variable: Abundance each category _________________________

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

# Generate marginal means 

# normalise each subclass lipid concentration to zero mean and unit variance for each age 

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

#_________________1 Day_____________

Abundance.global.1Day <- DF_glm01_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Percentage.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(id == "PLs") %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

pd = position_dodge(0.4)    ### How much to jitter the points on the plot

Abundance.global.1Dayplot.PL <- Abundance.global.1Day %>% 
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
  #ggtitle("Percent abundance Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1, 1, 1), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))  

#ggsave(plot = Abundance.global.1Dayplot.PL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure3/AB_PLglobalD1.jpg")                

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

#_____

Abundance.global.19Day <- DF_glm01_19 %>% 
  split(.$Type) %>% 
  map(~ lm(Percentage.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 4.59)) %>%  
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(id == "PLs") %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

Abundance.global.19Dayplot.PL <- Abundance.global.19Day %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.5, 1.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Percent abundance Day 19")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(1, 1, 1), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 

#ggsave(plot = Abundance.global.19Dayplot.PL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure3/AB_PLglobalD19.jpg") 

#________________________________________________________________________________________________________________


#___________________Variable: Mean carbon chain length in PL lipid category_________________________

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
                                            "NLs")))) %>% 
  filter(SubClass != "TG e" & SubClass != "PC e" & SubClass != "PE e" & SubClass != "PE p") # we remove all ELs

#__________________________

# Day 1 PL CC global mean reduced model

DF_glm02_1 <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1:SC4) %>% 
  tidyr::separate(SC1, into = c("SC1CC", "SC1DB"), remove = TRUE) %>% 
  tidyr::separate(SC2, into = c("SC2CC", "SC2DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("SC3CC", "SC3DB"), remove = TRUE) %>% 
  tidyr::separate(SC4, into = c("SC4CC", "SC4DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(SC1CC.Wt = as.numeric(SC1CC)*Conc) %>% 
  dplyr::mutate(SC2CC.Wt = as.numeric(SC2CC)*Conc) %>%
  dplyr::mutate(SC3CC.Wt = as.numeric(SC3CC)*Conc) %>% 
  dplyr::mutate(SC4CC.Wt = as.numeric(SC4CC)*Conc) %>%
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, SC1CC.Wt:SC4CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1CC.Wt', 'SC2CC.Wt', 'SC3CC.Wt', 'SC4CC.Wt'), names_to='Chain', values_to='ChainCC') %>% 
  dplyr::filter(ChainCC != "NA") %>% # to remove SC3 from DG
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(ChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Carbon = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Type, Mean.Carbon) %>% 
  unique() %>% 
  dplyr::filter(Age=="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.Carbon), SD=sd(Mean.Carbon)) %>% 
  mutate(Mean.Carbon.t = (Mean.Carbon-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) 

# Use a mean weight for all 7 stains 
Weights_Category_1 <- DF.main  %>%
  dplyr::filter(Age =="1 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 

#_____

CC.global.1Day.PL <- DF_glm02_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.Carbon.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "PLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

CC.global.1Dayplot.PL <- CC.global.1Day.PL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.0, 1.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Carbon content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.6, 0.6, 0.6), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue"))  

#ggsave(plot = CC.global.1Dayplot.PL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure3/CC_PLglobalD1.jpg") 


#______________________________________________________________________________________________


#_______________________Variable: Mean double bond in lipid categories ________________________

DF_glm03_1 <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1:SC4) %>% 
  tidyr::separate(SC1, into = c("SC1CC", "SC1DB"), remove = TRUE) %>% 
  tidyr::separate(SC2, into = c("SC2CC", "SC2DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("SC3CC", "SC3DB"), remove = TRUE) %>% 
  tidyr::separate(SC4, into = c("SC4CC", "SC4DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(SC1DB.Wt = as.numeric(SC1DB)*Conc) %>% 
  dplyr::mutate(SC2DB.Wt = as.numeric(SC2DB)*Conc) %>%
  dplyr::mutate(SC3DB.Wt = as.numeric(SC3DB)*Conc) %>% 
  dplyr::mutate(SC4DB.Wt = as.numeric(SC4DB)*Conc) %>%
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, SC1DB.Wt:SC4DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1DB.Wt', 'SC2DB.Wt', 'SC3DB.Wt', 'SC4DB.Wt'), names_to='Chain', values_to='ChainDB') %>% 
  dplyr::filter(ChainDB != "NA") %>% # to remove SC3 from DG
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(ChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Bond = Numerator/Denominator) %>% 
  ungroup() %>%
  dplyr::select(Samples, Line, SubClass, Age, Time, Type, Weight, Mean.Bond) %>% 
  unique() %>% 
  dplyr::filter(Age =="1 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.Bond), SD=sd(Mean.Bond)) %>% 
  mutate(Mean.Bond.t = (Mean.Bond-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

#_________________1 Day_____________

# Use a mean weight for all 7 stains 
Weights_Category_1 <- DF.main  %>%
  dplyr::filter(Age =="1 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 

#____

DB.global.1Day.PL <- DF_glm03_1 %>% 
  split(.$Type) %>%  
  map(~ lm(Mean.Bond.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 1.96)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "PLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)

DB.global.1Dayplot.PL <- DB.global.1Day.PL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.0, 1.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Double bond content Day 1")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.5, 0.5, 0.5), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 

#ggsave(plot = DB.global.1Dayplot.PL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure3/DB_PLglobalD1.jpg") 

#____________Day 19 ______________

DF_glm03_19 <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1:SC4) %>% 
  tidyr::separate(SC1, into = c("SC1CC", "SC1DB"), remove = TRUE) %>% 
  tidyr::separate(SC2, into = c("SC2CC", "SC2DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("SC3CC", "SC3DB"), remove = TRUE) %>% 
  tidyr::separate(SC4, into = c("SC4CC", "SC4DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(SC1DB.Wt = as.numeric(SC1DB)*Conc) %>% 
  dplyr::mutate(SC2DB.Wt = as.numeric(SC2DB)*Conc) %>%
  dplyr::mutate(SC3DB.Wt = as.numeric(SC3DB)*Conc) %>% 
  dplyr::mutate(SC4DB.Wt = as.numeric(SC4DB)*Conc) %>%
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, SC1DB.Wt:SC4DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1DB.Wt', 'SC2DB.Wt', 'SC3DB.Wt', 'SC4DB.Wt'), names_to='Chain', values_to='ChainDB') %>% 
  dplyr::filter(ChainDB != "NA") %>% # to remove SC3 from DG
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(ChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Bond = Numerator/Denominator) %>% 
  ungroup() %>%
  dplyr::select(Samples, Line, SubClass, Age, Time, Type, Weight, Mean.Bond) %>% 
  unique() %>% 
  dplyr::filter(Age =="19 Day") %>% 
  group_by(Type) %>% 
  mutate(MEAN= mean(Mean.Bond), SD=sd(Mean.Bond)) %>% 
  mutate(Mean.Bond.t = (Mean.Bond-MEAN)/SD) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

# Use a mean weight for all 7 stains i.e we include S06
Weights_Category_19 <- DF.main %>%
  dplyr::filter(Age =="19 Day") %>% 
  summarise(global_wt = mean(Weight), meanTop5pct = mean(Weight[Weight>=quantile(Weight, 0.90)]),
            meanBottom5pct = mean(Weight[Weight<=quantile(Weight, 0.10)])) 

#_____

DB.global.19Day.PL <- DF_glm03_19 %>% 
  split(.$Type) %>% 
  map(~ lm(Mean.Bond.t ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime), at = list(Weight = 4.59)) %>%  
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%  dplyr::filter(str_detect(id , "PLs")) %>% 
  tidyr::separate(col = LineTime, into = c("Line", "Time"), "_", remove = FALSE)


DB.global.19Dayplot.PL <- DB.global.19Day.PL %>% 
  dplyr::mutate(LineTime = factor(as.factor(LineTime), levels = c("ct_New","ct_Old","syd_New","syd_Old","s06_Older","cbr_New","cbr_Old"))) %>%
  ggplot(aes(x = LineTime, y = emmean, label = .group)) +
  geom_point(mapping = aes(color = Line, shape  = Time), size = 2, position = pd) +
  geom_errorbar(aes(ymin  =  lower.CL, ymax =  upper.CL, color = Line), width =  0.1, size  =  0.4, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(-1.0, 1.5), labels = scales::number_format(accuracy = 0.1)) +
  #facet_wrap(~id) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Double bond content Day 19")+ theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_text_repel(nudge_x = c(0, 0, 0), nudge_y = c(0.5, 0.5, 0.5), color = "black", segment.color = "transparent", size = 3)+
  theme(legend.position = "none")+
  theme(axis.title.x =element_blank()) + theme(axis.title.y =element_blank()) + theme(axis.text.x =element_blank()) +
  scale_colour_manual(values = c("ct" = "red", "cbr" = "green", "syd" = "blue", "s06" = "blue")) 

#ggsave(plot = DB.global.19Dayplot.PL, width = 50, height = 30, units = "mm", dpi = 300,filename = "ScriptManuscript2/Figures/Figure3/DB_PLglobalD19.jpg") 

#______________________________END_____________________________________________

