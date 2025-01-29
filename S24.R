#_________________________________________24/01/2025____________________________________
# Here we calculate lipid species abundances in each class for each of the seven strains 
# We use this analysis to generate Table 5 and Table S12.

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


#__________________Calculate percentage conc of lipid species in each class and strain________________________

# NB: for the complete function used to fill in zeros we will have SubClass and Name inside the 1st nesting unlike above when calculating the % abundance of each SubClass, we only had SubClass in the 1st nesting

Abn <- DF.main %>%
  ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Name, Conc) %>% 
  #dplyr::filter(Line == "s06" & Samples != "s06_158")%>% # remove this outlier) %>% 
  complete(nesting(SubClass,Name), nesting(Samples, Age, Line, Time), fill = list(Conc = 0)) %>% 
  dplyr::group_by(Samples, SubClass, Name) %>% 
  dplyr::mutate(Total.Conc.by.SAsC = sum(Conc)) %>% # by samples, age and subclass, species 
  dplyr::ungroup() %>% 
  #dplyr::group_by(Samples,SubClass) %>% # sum area for a sample for each subclass
  dplyr::group_by(Samples) %>% # sum area for a sample regardless of subclass (we want this code grouping by samples so we find overall average abundance over 15/16 replicates not only in samples where a particular lipid species is present)
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Name, Total.Conc.by.SAsC, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.SAsC/Total.Conc.by.SA)*100) %>% # for each sample the lipid specie sum to 100% sum((S06_lipids %>% filter(Samples == "s06_02"))$Percentage) 
  unique() %>% 
  dplyr::select(-Total.Conc.by.SAsC, -Total.Conc.by.SA, -Samples) %>%
  Rmisc::summarySE(groupvars = c("Line", "Time", "Age", "SubClass", "Name"), measurevar = "Percentage") %>% # when we take mean we are summarising for whole s06 line and not for sample
  dplyr::select(Line, Time, Age, SubClass, Name, Percentage) %>% 
  tidyr::unite(LineTime, Line, Time, remove = TRUE) %>%
  tidyr::pivot_wider(id_cols = c(SubClass, Name, Age), names_from = LineTime, values_from = Percentage) %>% 
  dplyr::select(SubClass,Name,Age,ct_New,ct_Old,syd_New,syd_Old,s06_Older,cbr_New,cbr_Old) 

Abn1 <- Abn %>% 
  dplyr::filter(Age == "1 Day") #%>% 
  #write_excel_csv("ScriptManuscript2/Commonchains/Abundance1_lipidspecies_in_each_strain.csv") 

Abn19 <- Abn %>% 
  dplyr::filter(Age == "19 Day") #%>% 
  #write_excel_csv("ScriptManuscript2/Commonchains/Abundance19_lipidspecies_in_each_strain.csv")

# check to see the sum of percentages within a subclass should be same as the one for that subclass percentage abundance obtained in Table 1 of lipids MS1
#sum((Abn %>% dplyr::filter(Age == "1 Day" & LineTime == "ct_New"))$Percentage)
#sum((Abn %>% dplyr::filter(Age == "1 Day" & SubClass == "PE" & LineTime == "cbr_New"))$Percentage)


#____________________________END______________________________


