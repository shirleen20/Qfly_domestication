#_________________________________________24/01/2025____________________________________
# Here we perform t-tests on variables percent Abundance
# We test two hypotheses about the patterns of lipidome changes during domestication. 
# Hypothesis 1 is that the three old strains are more similar to one another than the 
# three new strains are to one another and 
# hypothesis 2 is that the three old strains are more similar to SDolder than the three 
# new strains are to SDolder. 
# We use this analysis to generate Table 4 and Table S11.


library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("broom")
library("gtools")
library("multcomp")
library("multcompView")
library("Rmisc")

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

df.old <- read_csv("MS1Standards/Table_ConcAllStandards.csv") %>% # this is data from 1st run
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

#__________________________________________________________________________________________________________________

#______________________________T-tests on variable Abundance_______________________________


#___________________________Variable : Abundance each category _________________________

# complete fxn fills in zeros for samples that have missing lipids since the numbers we analyse for % abundances 
# are always using all replicates even if the lipid in question was not present in all of them 

Complete_DF <- DF.main %>% ungroup() %>% 
  complete(SubClass, nesting(Samples, Age, Line, Time, Weight), fill = list(Conc = 0)) %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Conc, Weight)

# normalise each subclass lipid concentration to zero mean and unit variance for each age 

DF_glm01 <- Complete_DF %>% 
  ungroup %>% 
  dplyr::group_by(Samples,SubClass) %>% 
  dplyr::mutate(Total.Conc.by.SAsC = sum(Conc)) %>% # by samples, age and subclass 
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>% # sum area for a sample regardless of subclass (we want this code grouping by samples so we find overall average abundance over 15/16 replicates not only in samples where a particular lipid species is present)
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Weight, Total.Conc.by.SAsC, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.SAsC/Total.Conc.by.SA)*100) %>%
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF01 %>% filter(Samples == "cbr_n04"))$Percentage)
  unique() %>% 
  dplyr::select(-Total.Conc.by.SAsC, -Total.Conc.by.SA) %>%
  #tidyr::unite(SubClass.Age, SubClass, Age, remove = FALSE) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(SubClass, Age, Samples, LineTime, Percentage, Weight) %>% 
  summarySE(groupvars = c("LineTime", "SubClass", "Age"), measurevar = "Percentage") %>% 
  dplyr::select(LineTime, SubClass, Age, Percentage) 


#___________________Test hypothesis 1______________________

# Hypothesis 1: Test whether olds are closer to one another than news 
# CBR vs CT, CBR vs SYD, CT vs SYD: For each class and each day calc the t-test a) for new pop b) for old pop so calc difference between the three pairs for new and old.

DF_glm01_D1 <- DF_glm01 %>%
  dplyr::filter(Age == "1 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Percentage) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Percentage") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old")) %>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Percentage, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Percentage) 

# T-test
M1_1 <- DF_glm01_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "Abn_H1D1")
#write_excel_csv("ScriptManuscript2/T.test/PercentabundanceD1_H1.csv")

# Do FDR corrections for M1_1

FDRcorrected.M1_1 <- M1_1 %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))

DF_glm01_D19 <- DF_glm01 %>%
  dplyr::filter(Age == "19 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Percentage) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Percentage") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Percentage, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Percentage) 


# T-test
M1_19 <- DF_glm01_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>%  
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "Abn_H1D19")
#write_excel_csv("ScriptManuscript2/T.test/PercentabundanceD19_H1.csv")

FDRcorrected.M1_19 <- M1_19 %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))

#_______________________________Test hypothesis 2_____________________________________________________________

DF_glm02_D1 <- DF_glm01 %>%
  dplyr::filter(Age == "1 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Percentage) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Percentage") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Percentage, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Percentage) 

# T-test
M2_1 <- DF_glm02_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>%  
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "Abn_H2D1")
#write_excel_csv("ScriptManuscript2/T.test/PercentabundanceD1_H2.csv")

FDRcorrected.M2_1 <- M2_1 %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


DF_glm02_D19 <- DF_glm01 %>%
  dplyr::filter(Age == "19 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Percentage) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Percentage") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Percentage, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Percentage) 

# T-test

M2_19 <- DF_glm02_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>%
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "Abn_H2D19")

FDRcorrected.M2_19 <- M2_19 %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))

# Save results
T.test.Abn <- rbind(M1_1, M1_19, M2_1, M2_19) #%>% 
  #write_excel_csv("ScriptManuscript2/T.test/Percentabundance.csv")

#_________________________________END____________________________________

