#_________________________________________24/01/2025____________________________________

# Here we perform univariate analysis (Anova lm) for variable abundance each category, we have four categories: 
# Neutral lipids, phospholipids, ether neutral lipids and ether phospholipids.
# we 1st fit a full model and check if the interaction term is sig after FDR corrections, if not than we fit a 
# reduced model and we also generate marginal means for variables where strain was significant
# We use this analysis to generate Tables 1, S7-S10.

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

#__________

DF_glm01 <- Complete_DF %>% 
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
  tidyr::unite(Type.Age, Type, Age, remove = FALSE) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>% 
  dplyr::select(Type.Age, Samples, LineTime, Percentage, Weight) 

by_Type.Age <- DF_glm01 %>% split(.$Type.Age)

# fit a full model to test significant terms
M1 <- by_Type.Age %>% 
  map(~ lm(Percentage ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) 

# fit a reduced model to test significant terms
M1.R <- by_Type.Age %>% 
  map(~ lm(Percentage ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        car::Anova(type = 3) %>% broom::tidy()) %>% bind_rows(.id = "id") %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".",""))))) %>% 
  dplyr::mutate(type = "Type") 

#________________________________

# Do FDR corrections for variables to see if the interaction term i.e. as.numeric(Weight):LineTime is significant, if not then we move to fitting reduced model ijn next script

DF.p.values <- M1 %>% 
  dplyr::select(id, term, 'p.value') %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value) 

FDRcorrected.Category.fullmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = p.adjust(DF.p.values$`as.numeric(Weight):LineTime`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight):LineTime`))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate( LineTime = stars.pval(LineTime)) %>% 
  dplyr::mutate(`as.numeric(Weight):LineTime` = stars.pval(`as.numeric(Weight):LineTime`))

# NOTE: NONE OF THE WEIGHT BY LINE INTERACTION TERM IS SIGNIFICANT SO WE WILL USE A REDUCED MODEL

# save Anova results after FDR corrections full model

rbind(FDRcorrected.Category.fullmodel) #%>% write_excel_csv("ScriptManuscript2/Abundance/FDRcorrected.Category.fullmodel.csv")

# save Anova results for full model 

M1 <- by_Type.Age %>% 
  map(~ lm(Percentage ~ as.numeric(Weight)* LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
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
  dplyr::mutate("Wt.Str" = paste0(Wt.Str,")"))# %>%
  #write_excel_csv("ScriptManuscript2/Abundance/Anovaresults.CategoryAbundance.fullmodel.csv")

#______

# Do FDR corrections for reduced model on all the terms i.e. (Intercept), as.numeric(Weight), LineTime

DF.p.values <- M1.R %>% 
  dplyr::select(id, term, 'p.value', type) %>% 
  filter(!str_detect(term, "Residuals")) %>% 
  pivot_wider(names_from = term, values_from = p.value)

FDRcorrected.Category.reducedmodel <- DF.p.values %>% 
  dplyr::mutate(`(Intercept)` = p.adjust(DF.p.values$`(Intercept)`, method = "fdr", n = length(DF.p.values$`(Intercept)`))) %>% 
  dplyr::mutate(`as.numeric(Weight)` = p.adjust(DF.p.values$`as.numeric(Weight)`, method = "fdr", n = length(DF.p.values$`as.numeric(Weight)`))) %>% 
  dplyr::mutate(LineTime = p.adjust(DF.p.values$LineTime, method = "fdr", n = length(DF.p.values$LineTime))) %>% 
  dplyr::mutate(`(Intercept)` = stars.pval(`(Intercept)`)) %>% 
  dplyr::mutate(`as.numeric(Weight)` = stars.pval(`as.numeric(Weight)`)) %>% 
  dplyr::mutate(LineTime = stars.pval(LineTime)) 

# save Anova results after FDR corrections reduced model
rbind(FDRcorrected.Category.reducedmodel) #%>% write_excel_csv("ScriptManuscript2/Abundance/FDRcorrected.Category.reducedmodel.csv")

# save Anova results for reduced model
M1.R <- by_Type.Age %>% 
  map(~ lm(Percentage ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
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
  dplyr::mutate(Str = paste0(Str,")")) #%>%
  #write_excel_csv("ScriptManuscript2/Abundance/Anovaresults.Category.reducedmodel.csv")

#_____________________________Generate marginal means______________________________

by_Type.Age <- DF_glm01 %>% split(.$Type.Age)

M1.R <- by_Type.Age %>% 
  map(~ lm(Percentage ~ as.numeric(Weight) + LineTime, data = .x, contrasts = list(LineTime = contr.sum)) %>% 
        emmeans(~(LineTime)) %>% 
        cld(Letter = "abcdefgh") %>% as_tibble()) %>% 
  bind_rows(.id = "id")  #%>%
  #write_excel_csv("ScriptManuscript2/Abundance/EmmeansCategoryReducedmodel.csv")

#________________________________________________________________________________________________


#_________________________________Master table with emmeans____________________

# Here we generate Master table from Anova and emmeans analysis 

#rm(list=ls())

#Emmeans <- readxl::read_excel("ScriptManuscript2/Abundance/EmmeansCategory.xlsx")

#DF.p.Emmeans <- Emmeans  

# A <- DF.p.Emmeans %>% 
#   dplyr::select(-SE, -df) %>% 
#   dplyr::mutate(emmean =  formatC(as.numeric(emmean), format = "e", digits = 2)) %>% 
#   dplyr::mutate(lower.CL =  formatC(as.numeric(lower.CL), format = "e", digits = 2)) %>%
#   dplyr::mutate(upper.CL =  formatC(as.numeric(upper.CL), format = "e", digits = 2)) %>%      
#   dplyr::mutate(Unit = str_sub(emmean,-4,-1)) %>% 
#   dplyr::mutate(emmean = str_sub(emmean, end=-5)) %>% 
#   dplyr::mutate(lower.CL = str_sub(lower.CL, end=-5)) %>% 
#   dplyr::mutate(upper.CL = str_sub(upper.CL, end=-5)) %>% 
#   dplyr::mutate(left = "(", right = ")") %>% 
#   tidyr::unite("A", "lower.CL", "upper.CL", sep = "-") %>% 
#   tidyr::unite(B, "left", "A", "right", sep = "") %>%  
#   tidyr::unite("C", "emmean", "B", sep = "")    %>% 
#   tidyr::unite("P", "C", ".group", sep = " ")    %>%
#   dplyr::filter(Unit !=  " NA") %>% 
#   tidyr::pivot_wider(names_from = LineTime, values_from = P) #%>%
  #write_excel_csv("ScriptManuscript2/Abundance/MastertableCategory.csv")

#______________________________END___________________________________
