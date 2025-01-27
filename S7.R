#____________________________________24/01/2025________________________________________

# Here we perform univariate analysis (Anova lm) on variables Carbon content and Double bond content for each ester 
# lipid class in neutral and phospholipids
# we 1st fit a full model and check if the interaction term is sig after FDR corrections, if not than we fit a reduced model 
# we do FDR corrections separately for all the variables and we extract emmeans for variables where strain was significant
# We use this analysis to generate Tables 1, S7 and S8.

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

#___________________________Variable: Mean carbon chain length analysis_________________________

#_____________create dataframe for mean carbon and double bond calculation in lipid category_____________________

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
                                            "NLs"))))

#_________________________________

DF_glm01 <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1:SC4) %>% 
  dplyr::filter(Type != 'ELNLs' & Type != "ELPLs") %>% 
  tidyr::separate(SC1, into = c("SC1CC", "SC1DB"), remove = TRUE) %>% 
  tidyr::separate(SC2, into = c("SC2CC", "SC2DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("SC3CC", "SC3DB"), remove = TRUE) %>% 
  tidyr::separate(SC4, into = c("SC4CC", "SC4DB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(SC1CC.Wt = as.numeric(SC1CC)*Conc) %>% 
  dplyr::mutate(SC2CC.Wt = as.numeric(SC2CC)*Conc) %>%
  dplyr::mutate(SC3CC.Wt = as.numeric(SC3CC)*Conc) %>% 
  dplyr::mutate(SC4CC.Wt = as.numeric(SC4CC)*Conc) %>%
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Weight, SC1CC.Wt:SC4CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1CC.Wt', 'SC2CC.Wt', 'SC3CC.Wt', 'SC4CC.Wt'), names_to='Chain', values_to='ChainCC') %>% 
  dplyr::filter(ChainCC != "NA") %>% # to remove SC3 from DG
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(ChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Carbon = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Weight, Mean.Carbon) %>% 
  unique() %>% 
  tidyr::unite(SubClass.Age, SubClass, Age, remove = FALSE) %>% 
  dplyr::select(SubClass.Age, Samples, Line, SubClass, Age, Time, Mean.Carbon, Weight) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) 

by_SubClass.Age <- DF_glm01 %>% split(.$SubClass.Age)

# use full model to test significant interaction terms
M1 <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::mutate(type = "CC")

# fit a reduced model to test significant terms
M1.R <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::mutate(type = "CC")

# Generate emmeans:
M1.R.Emmeans <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~(LineTime)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") 

# Since DG 1 Day for carbon chain length analysis has significant interaction term, we need to extract emmeans using full model

DG_1Day <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight) * LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~interaction(LineTime)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id") %>%
  filter(str_detect(id, "DG_1 Day")) #%>%
  #write_excel_csv("ScriptManuscript2/Carbon/Emmeans.DG1CC.csv")

# we also analyse emtrends when Weight:LineTime is significant and present under emmeans data in final results table

M1.emtrends <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emtrends(pairwise~LineTime, var = "Weight") %>% 
        cld(Letter = "abcdefgh") %>% 
        as.data.frame()) %>% 
  bind_rows(.id = "id") %>%  
  dplyr::mutate(type = "CC") %>% dplyr::filter(id %in% c("DG_1 Day")) #%>%
  #write_excel_csv("ScriptManuscript2/Carbon/Emtrends.DG1CC.csv")

# Do FDR corrections for FULL MODEL to see if the interaction term i.e. as.numeric(Weight):LineTime is significant, if not then we move to fitting reduced model ijn next script

DF.p.values <- M1 %>% 
  dplyr::select(id, term, 'p.value') %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value) 

FDRcorrected.CC.fullmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = p.adjust(DF.p.values$`as.numeric(Weight):LineTime`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight):LineTime`))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate( LineTime = stars.pval(LineTime)) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = stars.pval(`as.numeric(Weight):LineTime`)) %>% 
  dplyr::filter(str_detect(id, "DG_1 Day")) 

# save Anova results after FDR corrections full model
rbind(FDRcorrected.CC.fullmodel) #%>% write_excel_csv("ScriptManuscript2/Carbon/FDRcorrected.fullmodel.CC.csv")

# save Anova results for full model 

M1 <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  mutate(statistic = format(round(.$statistic, 1), nsmall = 1)) %>%
  dplyr::select(id,term,df,statistic) %>% 
  dplyr::mutate(F = statistic) %>% 
  dplyr::select(-statistic) %>% 
  pivot_wider(names_from = term, values_from = c(df, F)) %>% 
  dplyr::select(-2,-7,-11) %>%
  unite("Wt.", "df_as.numeric(Weight)" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("Str.", "df_LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("WtStr.", "df_as.numeric(Weight):LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>%
  dplyr::select(-"df_as.numeric(Weight)",-"df_LineTime",-"df_as.numeric(Weight):LineTime", -"df_Residuals") %>% 
  unite("Wt", "F_as.numeric(Weight)", "Wt.", sep = "(") %>% 
  dplyr::mutate(Wt = paste0(Wt,")")) %>%  
  unite("Str", "F_LineTime", "Str." , sep = "(") %>% 
  dplyr::mutate(Str = paste0(Str,")")) %>%
  unite("Wt.Str", "F_as.numeric(Weight):LineTime", "WtStr.", sep = "(") %>% 
  dplyr::mutate("Wt.Str" = paste0(Wt.Str,")")) %>%
  dplyr::filter(str_detect(id, "DG_1 Day")) #%>% 
  #write_excel_csv("ScriptManuscript2/Carbon/Anovaresults.fullmodel.CC.csv")

# SAVE EMMEANS RESULTS for reduced models

rbind(M1.R.Emmeans) #%>% write_excel_csv("ScriptManuscript2/Carbon/EmmeansCC.ester.csv")

#_____________________

# Do FDR corrections for REDUCED MODEL on all the terms i.e. (Intercept), as.numeric(Weight), LineTime

DF.p.values <- M1.R %>% 
  dplyr::select(id, term, 'p.value') %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value) 

FDRcorrected.CC.reducedmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate(LineTime = stars.pval(LineTime)) %>% 
  dplyr::filter(id != "DG_1 Day")

# save Anova results after FDR corrections full model
rbind(FDRcorrected.CC.reducedmodel) #%>% write_excel_csv("ScriptManuscript2/Carbon/FDRcorrected.reducedmodel.CC.csv")

# save Anova results for reduced model 

M1.R <- by_SubClass.Age %>% 
  map(~ lm(Mean.Carbon  ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  mutate(statistic = format(round(.$statistic, 1), nsmall = 1)) %>%
  dplyr::select(id,term,df,statistic) %>% 
  dplyr::mutate(F = statistic) %>% 
  dplyr::select(-statistic) %>% 
  pivot_wider(names_from = term, values_from = c(df, F)) %>% 
  dplyr::select(-2,-6,-9) %>%
  unite("Wt.", "df_as.numeric(Weight)" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("Str.", "df_LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  dplyr::select(-"df_as.numeric(Weight)",-"df_LineTime", -"df_Residuals") %>% 
  unite("Wt", "F_as.numeric(Weight)", "Wt.", sep = "(") %>% 
  dplyr::mutate(Wt = paste0(Wt,")")) %>%  
  unite("Str", "F_LineTime", "Str." , sep = "(") %>% 
  dplyr::mutate(Str = paste0(Str,")")) %>% 
  dplyr::filter(id != "DG_1 Day") #%>% 
  #write_excel_csv("ScriptManuscript2/Carbon/Anovaresults.reducedmodel.CC.csv")


#___________________________Master table for CC ester reduced model______________________________


rm(list=ls())

#_____________________
Emmeans <- readxl::read_excel("ScriptManuscript2/Carbon/EmmeansCC.ester.xlsx")

DF.p.Emmeans <- Emmeans 

A <- DF.p.Emmeans %>% 
  dplyr::select(-SE, -df) %>% 
  dplyr::mutate(emmean =  formatC(as.numeric(emmean), format = "e", digits = 2)) %>% 
  dplyr::mutate(lower.CL =  formatC(as.numeric(lower.CL), format = "e", digits = 2)) %>%
  dplyr::mutate(upper.CL =  formatC(as.numeric(upper.CL), format = "e", digits = 2)) %>%      
  dplyr::mutate(Unit = str_sub(emmean,-4,-1)) %>% 
  dplyr::mutate(emmean = str_sub(emmean, end=-5)) %>% 
  dplyr::mutate(lower.CL = str_sub(lower.CL, end=-5)) %>% 
  dplyr::mutate(upper.CL = str_sub(upper.CL, end=-5)) %>% 
  dplyr::mutate(left = "(", right = ")") %>% 
  tidyr::unite("A", "lower.CL", "upper.CL", sep = "-") %>% 
  tidyr::unite(B, "left", "A", "right", sep = "") %>%  
  tidyr::unite("C", "emmean", "B", sep = "")    %>% 
  tidyr::unite("P", "C", ".group", sep = " ")    %>%
  dplyr::filter(Unit !=  " NA") %>% 
  tidyr::pivot_wider(names_from = LineTime, values_from = P) #%>% 
  #write_excel_csv("ScriptManuscript2/Carbon/MastertableCC.ester.csv")

#_______________# Here we generate Master table for emmeans analysis full model 1 day DG CC ester lipids_____________________

rm(list=ls())

Emmeans <- readxl::read_excel("ScriptManuscript2/Carbon/Emmeans.DG1CC.xlsx")

DF.p.Emmeans <- Emmeans 

A <- DF.p.Emmeans %>% 
  dplyr::select(-SE, -df) %>% 
  dplyr::mutate(emmean =  formatC(as.numeric(emmean), format = "e", digits = 2)) %>% 
  dplyr::mutate(lower.CL =  formatC(as.numeric(lower.CL), format = "e", digits = 2)) %>%
  dplyr::mutate(upper.CL =  formatC(as.numeric(upper.CL), format = "e", digits = 2)) %>%      
  dplyr::mutate(Unit = str_sub(emmean,-4,-1)) %>% 
  dplyr::mutate(emmean = str_sub(emmean, end=-5)) %>% 
  dplyr::mutate(lower.CL = str_sub(lower.CL, end=-5)) %>% 
  dplyr::mutate(upper.CL = str_sub(upper.CL, end=-5)) %>% 
  dplyr::mutate(left = "(", right = ")") %>% 
  tidyr::unite("A", "lower.CL", "upper.CL", sep = "-") %>% 
  tidyr::unite(B, "left", "A", "right", sep = "") %>%  
  tidyr::unite("C", "emmean", "B", sep = "")    %>% 
  tidyr::unite("P", "C", ".group", sep = " ")    %>%
  dplyr::filter(Unit !=  " NA") %>% 
  tidyr::pivot_wider(names_from = LineTime, values_from = P) #%>% 
  #write_excel_csv("ScriptManuscript2/Carbon/MastertableDG1CC.csv")

#_______________

rm(list=ls())

# Here we generate Master table for emtrends analysis for 1 day DG CC 

Emtrends <- readxl::read_excel("ScriptManuscript2/Carbon/Emtrends.DG1CC.xlsx")

DF.p.Emtrends <- Emtrends 

A <- DF.p.Emtrends %>% 
  dplyr::select(-SE, -df) %>% 
  dplyr::mutate(Weight.trend =  formatC(as.numeric(Weight.trend), format = "e", digits = 2)) %>% 
  dplyr::mutate(lower.CL =  formatC(as.numeric(lower.CL), format = "e", digits = 2)) %>%
  dplyr::mutate(upper.CL =  formatC(as.numeric(upper.CL), format = "e", digits = 2)) %>%      
  dplyr::mutate(Unit = str_sub(Weight.trend,-4,-1)) %>% 
  dplyr::mutate(Weight.trend = str_sub(Weight.trend, end=-5)) %>% 
  dplyr::mutate(lower.CL = str_sub(lower.CL, end=-5)) %>% 
  dplyr::mutate(upper.CL = str_sub(upper.CL, end=-5)) %>% 
  dplyr::mutate(left = "(", right = ")") %>% 
  tidyr::unite("A", "lower.CL", "upper.CL", sep = "-") %>% 
  tidyr::unite(B, "left", "A", "right", sep = "") %>%  
  tidyr::unite("C", "Weight.trend", "B", sep = "")    %>% 
  tidyr::unite("P", "C", ".group", sep = " ")    %>%
  dplyr::filter(Unit !=  " NA") %>% 
  tidyr::pivot_wider(names_from = LineTime, values_from = P) #%>% 
  #write_excel_csv("ScriptManuscript2/Carbon/MastertableEmtrends.DG1CC.csv")

#_______________________________________________________________________________________________________________________________________


#___________________________Variable: Mean double bond analysis_________________________

DF_glm02 <-  Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Weight, Conc, SC1:SC4) %>% 
  dplyr::filter(Type != 'ELNLs' & Type != "ELPLs") %>% 
  tidyr::separate(SC1, into = c("SC1CC", "SC1DB"), remove = TRUE) %>% 
  tidyr::separate(SC2, into = c("SC2CC", "SC2DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("SC3CC", "SC3DB"), remove = TRUE) %>% 
  tidyr::separate(SC4, into = c("SC4CC", "SC4DB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(SC1DB.Wt = as.numeric(SC1DB)*Conc) %>% 
  dplyr::mutate(SC2DB.Wt = as.numeric(SC2DB)*Conc) %>%
  dplyr::mutate(SC3DB.Wt = as.numeric(SC3DB)*Conc) %>% 
  dplyr::mutate(SC4DB.Wt = as.numeric(SC4DB)*Conc) %>%
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Weight, SC1DB.Wt:SC4DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1DB.Wt', 'SC2DB.Wt', 'SC3DB.Wt', 'SC4DB.Wt'), names_to='Chain', values_to='ChainDB') %>% 
  dplyr::filter(ChainDB != "NA") %>% # to remove SC3 from DG
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(ChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Bond = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Weight, Mean.Bond) %>% 
  unique() %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  tidyr::unite(SubClass.Age, SubClass, Age, remove = FALSE) %>% 
  dplyr::select(SubClass.Age, Samples, Line, SubClass, Age, Time, Mean.Bond, Weight) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE)

by_SubClass.Age <- DF_glm02 %>% split(.$SubClass.Age)

# use full model to test significant interaction terms
M2 <- by_SubClass.Age %>% 
  map(~ lm(Mean.Bond ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", ifelse(p.value < 0.1, ".","")))))

# fit a reduced model to test significant terms
M2.R <- by_SubClass.Age %>% 
  map(~ lm(Mean.Bond ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", ifelse(p.value < 0.1, ".","")))))

M2.R.Emmeans <- by_SubClass.Age %>% 
  map(~ lm(Mean.Bond ~ as.numeric(Weight)+ LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~(LineTime)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id")  

# SAVE EMMEANS RESULTS for reduced models

rbind(M2.R.Emmeans) #%>% write_excel_csv("ScriptManuscript2/DoubleBond/EmmeansDBester.csv")


#_____________________________________

# Do FDR corrections for FULL MODEL to see if the interaction term i.e. as.numeric(Weight):LineTime is significant, if not then we move to fitting reduced model ijn next script

DF.p.values <- M2 %>% 
  dplyr::select(id, term, 'p.value') %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value) 

FDRcorrected.DB.fullmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = p.adjust(DF.p.values$`as.numeric(Weight):LineTime`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight):LineTime`))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate( LineTime = stars.pval(LineTime)) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = stars.pval(`as.numeric(Weight):LineTime`))

# save Anova results after FDR corrections full model
rbind(FDRcorrected.DB.fullmodel) #%>% write_excel_csv("ScriptManuscript2/DoubleBond/FDRcorrected.fullmodel.DB.csv")

# Save Anova results for full model
M2 <- by_SubClass.Age %>% 
  map(~ lm(Mean.Bond ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", ifelse(p.value < 0.1, ".",""))))) %>% 
  mutate(statistic = format(round(.$statistic, 1), nsmall = 1)) %>%
  dplyr::select(id,term,df,statistic) %>% 
  dplyr::mutate(F = statistic) %>% 
  dplyr::select(-statistic) %>% 
  pivot_wider(names_from = term, values_from = c(df, F)) %>% 
  dplyr::select(-2,-7,-11) %>%
  unite("Wt.", "df_as.numeric(Weight)" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("Str.", "df_LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("WtStr.", "df_as.numeric(Weight):LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>%
  dplyr::select(-"df_as.numeric(Weight)",-"df_LineTime",-"df_as.numeric(Weight):LineTime", -"df_Residuals") %>% 
  unite("Wt", "F_as.numeric(Weight)", "Wt.", sep = "(") %>% 
  dplyr::mutate(Wt = paste0(Wt,")")) %>%  
  unite("Str", "F_LineTime", "Str." , sep = "(") %>% 
  dplyr::mutate(Str = paste0(Str,")")) %>%
  unite("Wt.Str", "F_as.numeric(Weight):LineTime", "WtStr.", sep = "(") %>% 
  dplyr::mutate("Wt.Str" = paste0(Wt.Str,")")) #%>%
#write_excel_csv("ScriptManuscript2/DoubleBond/Anovaresults.fullmodel.DB.csv")


#_________________________________________

# Do FDR corrections for REDUCED MODEL on all the terms i.e. (Intercept), as.numeric(Weight), LineTime

DF.p.values <- M2.R %>% 
  dplyr::select(id, term, 'p.value') %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value) 

FDRcorrected.DB.reducedmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate( LineTime = stars.pval(LineTime)) 

# save Anova results after FDR corrections reduced model
rbind(FDRcorrected.DB.reducedmodel) %>% write_excel_csv("ScriptManuscript2/DoubleBond/FDRcorrected.reducedmodel.DB.csv")

# save Anova results for reduced model 
M2.R <- by_SubClass.Age %>% 
  map(~ lm(Mean.Bond ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*", ifelse(p.value < 0.1, ".",""))))) %>% 
  mutate(statistic = format(round(.$statistic, 1), nsmall = 1)) %>%
  dplyr::select(id,term,df,statistic) %>% 
  dplyr::mutate(F = statistic) %>% 
  dplyr::select(-statistic) %>% 
  pivot_wider(names_from = term, values_from = c(df, F)) %>% 
  dplyr::select(-2,-6,-9) %>%
  unite("Wt.", "df_as.numeric(Weight)" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  unite("Str.", "df_LineTime" , "df_Residuals", sep = ",", remove = FALSE) %>% 
  dplyr::select(-"df_as.numeric(Weight)",-"df_LineTime", -"df_Residuals") %>% 
  unite("Wt", "F_as.numeric(Weight)", "Wt.", sep = "(") %>% 
  dplyr::mutate(Wt = paste0(Wt,")")) %>%  
  unite("Str", "F_LineTime", "Str." , sep = "(") %>% 
  dplyr::mutate(Str = paste0(Str,")")) #%>% 
  #write_excel_csv("ScriptManuscript2/DoubleBond/Anovaresults.reducedmodel.DB.csv")


#____________________________________Create Master table for DB______________________________


rm(list=ls())

#_____________________
Emmeans <- readxl::read_excel("ScriptManuscript2/DoubleBond/EmmeansDBester.xlsx")

DF.p.Emmeans <- Emmeans 

A <- DF.p.Emmeans %>% 
  dplyr::select(-SE, -df) %>% 
  dplyr::mutate(emmean =  formatC(as.numeric(emmean), format = "e", digits = 2)) %>% 
  dplyr::mutate(lower.CL =  formatC(as.numeric(lower.CL), format = "e", digits = 2)) %>%
  dplyr::mutate(upper.CL =  formatC(as.numeric(upper.CL), format = "e", digits = 2)) %>%      
  dplyr::mutate(Unit = str_sub(emmean,-4,-1)) %>% 
  dplyr::mutate(emmean = str_sub(emmean, end=-5)) %>% 
  dplyr::mutate(lower.CL = str_sub(lower.CL, end=-5)) %>% 
  dplyr::mutate(upper.CL = str_sub(upper.CL, end=-5)) %>% 
  dplyr::mutate(left = "(", right = ")") %>% 
  tidyr::unite("A", "lower.CL", "upper.CL", sep = "-") %>% 
  tidyr::unite(B, "left", "A", "right", sep = "") %>%  
  tidyr::unite("C", "emmean", "B", sep = "")    %>% 
  tidyr::unite("P", "C", ".group", sep = " ")    %>%
  dplyr::filter(Unit !=  " NA") %>% 
  tidyr::pivot_wider(names_from = LineTime, values_from = P) #%>% 
  #write_excel_csv("ScriptManuscript2/DoubleBond/MastertableDB.csv")

#__________________________END_______________________________________

