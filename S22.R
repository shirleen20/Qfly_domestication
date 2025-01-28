#_________________________________________24/01/2025____________________________________
# Here we perform t-tests on variables carbon content n alkyl and acyl chains of ether lipid classes 
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

#__________________________________________________________________________________________________________________


#___________________________Variable: Mean carbon content analysis in ELs_________________________


# Create dataframe to calculate mean CC and mean DB in the akyl and acyl chains of each ether linked classes with disambiguated chains:

DFnew <- DF.main %>% dplyr::filter(SubClass != "PS e" & SubClass != "DG e") #remove subclasses PSe and DGe from DF because these two subclasses have all lipid species that do not have acyl chain identified

#grep - finds things, gives index 1, 3, 5 not showing 1 is apple, 3 is pear uless you specify value =T) check help ?gsub or grep
#gsub is a substitution

### Function to create for each subclass data frame, all side chains separated into columns
# X = data frame to process
# Y = Subclass level as string

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


# X is the data frame to process
# Subclass is the name of the subclass variable as a string (default: "SubClass")
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

#view(LIST1[[1]])

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
                                            "NLs"))))%>% 
  dplyr::filter(SubClass %in% c("TG e", "PE p", "PE e"))  # we remove PCe from analysis because there is no variation between samples for mean 

#__________________
# Variable: Mean carbon content in alkyl chains of ether lipid classes 

CCAlkyl <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Weight, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(AlkylCC.Wt = as.numeric(AlkylCC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Weight, AlkylCC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylCC.Wt'), names_to='AlkylChain', values_to='AlkylChainCC') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Weight, Mean.CarbonAlkyl) %>% 
  unique() %>%
  tidyr::unite(SubClass.Age, SubClass, Age, remove = FALSE) %>% 
  dplyr::select(SubClass.Age, Samples, Line, SubClass, Age, Time, Mean.CarbonAlkyl, Weight) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  summarySE(groupvars = c("LineTime", "SubClass", "Age"), measurevar = "Mean.CarbonAlkyl") %>% 
  dplyr::select(LineTime, SubClass, Age, Mean.CarbonAlkyl) 

#___________________Test hypothesis 1______________________

# Hypothesis 1: Test whether olds are closer to one another than news 
# CBR vs CT, CBR vs SYD, CT vs SYD: For each class and each day calc the t-test a) for new pop b) for old pop so calc difference between the three pairs for new and old.

DF_glm01_D1 <- CCAlkyl %>%
  dplyr::filter(Age == "1 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAlkyl) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Mean.CarbonAlkyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old")) %>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAlkyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAlkyl) 


M1_1CCAlkyl <- DF_glm01_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAlkyl_H1D1")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CCAlkylD1_H1.csv")

# Do FDR corrections for M1_1

FDRcorrected.M1_1CCAlkyl <- M1_1CCAlkyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))

DF_glm01_D19 <- CCAlkyl %>%
  dplyr::filter(Age == "19 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAlkyl) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Mean.CarbonAlkyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old")) %>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAlkyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAlkyl) 

# T-test
M1_19CCAlkyl <- DF_glm01_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAlkyl_H1D19")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CCAlkylD19_H1.csv")

FDRcorrected.M1_19CCAlkyl <- M1_19CCAlkyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


#___________________Test hypothesis 2_____________________________________________________________

DF_glm02_D1 <- CCAlkyl %>%
  dplyr::filter(Age == "1 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAlkyl) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Mean.CarbonAlkyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAlkyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAlkyl) 

# T-test
M2_1CCAlkyl <- DF_glm02_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAlkyl_H2D1")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CarbonAlkylD1_H2.csv")

FDRcorrected.M2_1CCAlkyl <- M2_1CCAlkyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


DF_glm02_D19 <- CCAlkyl %>%
  dplyr::filter(Age == "19 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAlkyl) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Mean.CarbonAlkyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAlkyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAlkyl) 

# T-test
M2_19CCAlkyl <- DF_glm02_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAlkyl_H2D19")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CarbonAlkylD19_H2.csv")

FDRcorrected.M2_19CCAlkyl <- M2_19CCAlkyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


#_________________________________________________________________________________________________________________
#_________________________________________________________________________________________________________________

#___________________Variable: Mean carbon content in acyl chains of ether lipid classes____________________________ 

CCAcyl <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Conc, Weight, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(Acyl1CC.Wt = as.numeric(Acyl1CC)*Conc) %>% 
  dplyr::mutate(Acyl2CC.Wt = as.numeric(Acyl2CC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Weight, Acyl1CC.Wt, Acyl2CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1CC.Wt', 'Acyl2CC.Wt'), names_to='AcylChain', values_to='AcylChainCC') %>% 
  dplyr::filter(AcylChainCC != "NA") %>% # to remove SC3 from DG
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AcylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Weight, Mean.CarbonAcyl) %>% 
  unique() %>%
  tidyr::unite(SubClass.Age, SubClass, Age, remove = FALSE) %>% 
  dplyr::select(SubClass.Age, Samples, Line, SubClass, Age, Time, Mean.CarbonAcyl, Weight) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  summarySE(groupvars = c("LineTime", "SubClass", "Age"), measurevar = "Mean.CarbonAcyl") %>% 
  dplyr::select(LineTime, SubClass, Age, Mean.CarbonAcyl)

#___________________Test hypothesis 1______________________

# Hypothesis 1: Test whether olds are closer to one another than news 
# CBR vs CT, CBR vs SYD, CT vs SYD: For each class and each day calc the t-test a) for new pop b) for old pop so calc difference between the three pairs for new and old.

DF_glm03_D1 <- CCAcyl %>%
  dplyr::filter(Age == "1 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAcyl) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Mean.CarbonAcyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old")) %>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAcyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAcyl) 


M1_1CCAcyl <- DF_glm03_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAcyl_H1D1")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CCAcylD1_H1.csv")

FDRcorrected.M1_1CCAcyl <- M1_1CCAcyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


DF_glm03_D19 <- CCAcyl %>%
  dplyr::filter(Age == "19 Day") %>% 
  dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAcyl) %>% 
  dplyr::mutate(CN.SD_new = abs(cbr_New-syd_New), CN.CT_new = abs(cbr_New-ct_New), CT.SD_new = abs(ct_New-syd_New)) %>% 
  dplyr::mutate(CN.SD_old = abs(cbr_Old-syd_Old), CN.CT_old = abs(cbr_Old-ct_Old), CT.SD_old = abs(ct_Old-syd_Old)) %>%
  dplyr::select(-2:-8) %>% 
  pivot_longer(cols = CN.SD_new:CT.SD_old, names_to = "Comparisons", values_to="Mean.CarbonAcyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.SD_new', 'CN.CT_new', 'CT.SD_new') ~ "New",
                                 Comparisons %in% c("CN.SD_old","CN.CT_old","CT.SD_old") ~ "Old")) %>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAcyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAcyl) 

# T-test
M1_19CCAcyl <- DF_glm03_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAcyl_H1D19")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CCAcylD19_H1.csv")

FDRcorrected.M1_19CCAcyl <- M1_19CCAcyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


#___________________Test hypothesis 2_____________________________________________________________

DF_glm04_D1 <- CCAcyl %>%
  dplyr::filter(Age == "1 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAcyl) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Mean.CarbonAcyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAcyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAcyl) 

# T-test
M2_1CCAcyl <- DF_glm04_D1 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAcyl_H2D1")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CarbonAcylD1_H2.csv")

FDRcorrected.M2_1CCAcyl <- M2_1CCAcyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))


DF_glm04_D19 <- CCAcyl %>%
  dplyr::filter(Age == "19 Day") %>% 
  #dplyr::filter(LineTime!= "s06_Older") %>% 
  pivot_wider(names_from = c(LineTime), values_from = Mean.CarbonAcyl) %>% 
  dplyr::mutate(CN.S06_new = abs(cbr_New-s06_Older), CT.S06_new = abs(ct_New-s06_Older), SD.S06_new = abs(syd_New-s06_Older)) %>% 
  dplyr::mutate(CN.S06_old = abs(cbr_Old-s06_Older), CT.S06_Old = abs(ct_Old-s06_Older), SD.S06_Old = abs(syd_Old-s06_Older)) %>% 
  dplyr::select(-2:-9) %>% 
  pivot_longer(cols = CN.S06_new:SD.S06_Old, names_to = "Comparisons", values_to="Mean.CarbonAcyl") %>% 
  dplyr::mutate(Time = case_when(Comparisons %in% c('CN.S06_new', 'CT.S06_new', 'SD.S06_new') ~ "New",
                                 Comparisons %in% c("CN.S06_old","CT.S06_Old","SD.S06_Old") ~ "Old"))%>%
  tidyr::separate(col = Comparisons, into = c("Line", NA), "_", remove = TRUE) %>% 
  dplyr::select(SubClass, Line, Mean.CarbonAcyl, Time) %>% 
  pivot_wider(names_from = c(Time), values_from = Mean.CarbonAcyl) 

# T-test
M2_19CCAcyl <- DF_glm04_D19 %>% split(.$SubClass) %>% 
  map(~ t.test(.$New, .$Old, paired=TRUE, conf.level=0.95, data = .x) %>% 
        broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::select(1, 2, 4, 6, 7, 10) %>% 
  dplyr::mutate(Test = "CCAcyl_H2D19")
#write_excel_csv("ScriptManuscript2/T.test/Mean.CarbonAcylD19_H2.csv")

FDRcorrected.M2_19CCAcyl <- M2_19CCAcyl %>% 
  dplyr::select(id, 'p.value') %>% 
  dplyr::mutate(p.value = p.adjust((p.value), method = "fdr", n = length(p.value))) %>% 
  dplyr::mutate(p.value = stars.pval(p.value))

# Save results
T.test.CC <- rbind(M1_1CCAlkyl, M1_19CCAlkyl, M2_1CCAlkyl, M2_19CCAlkyl,
                   M1_1CCAcyl, M1_19CCAcyl, M2_1CCAcyl, M2_19CCAcyl) #%>% 
  #write_excel_csv("ScriptManuscript2/T.test/CC_ether.csv")

#_______________________________________END______________________________________________________


