#_________________________________________24/01/2025____________________________________

# Here we calculate diversities, abundances, mean Carbon and mean double bond content
# in the seven strains to generate Figure 1, Table S3-S6 and we make correlation plots
# to generate Figure 6 and Figure S2

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("gtools")
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

#__________________________Calc conc of each lipid species in each sample____________________________________

DFConc1 <- DF.main %>% 
  ungroup %>% 
  dplyr::select(Samples, SubClass, Name, Conc) %>% 
  unique() %>% 
  pivot_wider(names_from = c(Samples), values_from = Conc) 

#__________________________Percent abundance in the seven strains_______________________________

# Calculate percentage conc and number of lipid species detected in each of the seven lines

# complete fxn fills in zeros for samples that have missing lipids

Complete_DF <- DF.main %>% ungroup() %>% 
  complete(SubClass, nesting(Samples, Age, Line, Time, Weight), fill = list(Conc = 0)) %>%
  dplyr::select(Samples, Age, Line, Time, SubClass, Conc, Weight)


DF01 <- Complete_DF %>% 
  ungroup %>% 
  dplyr::group_by(Samples,SubClass) %>% 
  dplyr::mutate(Total.Conc.by.SAsC = sum(Conc)) %>% # by samples, age and subclass 
  dplyr::ungroup() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Total.Conc.by.SA = sum(Conc)) %>% 
  dplyr::select(-Conc) %>%
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Age, Line, Time, SubClass, Total.Conc.by.SAsC, Total.Conc.by.SA) %>%
  dplyr::mutate(Percentage = (Total.Conc.by.SAsC/Total.Conc.by.SA)*100) %>% 
  unique() %>% 
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF01 %>% filter(Samples == "cbr_n04"))$Percentage)
  dplyr::select(-Total.Conc.by.SAsC, -Total.Conc.by.SA, -Samples) %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Percentage") %>% 
  dplyr::select(Line, Time, Age, SubClass, Percentage, se)

# Check to see if % adds to 100 for an individual sample, run the above codes till unique() the total % should come to 100. 
sum((DF01 %>% dplyr::filter(Samples == "s06_02"))$Percentage) #Note this works if you run codes till unique() for object DF01

# Check to see if percentage adds to 100 for a Line, Age, Time
DF01 %>% 
  dplyr::filter(Line == "ct" & Time == "New" & Age == "19 Day") %>% 
  dplyr::select(Percentage) %>% 
  sum()

# format and save as csv to produce Table S4
DF00 <- DF01 %>% 
  dplyr::mutate(Percentage = format(round(.$Percentage, 2), nsmall = 2)) %>% 
  dplyr::mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("PeSE", Percentage:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = PeSE) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>% 
#write.csv("ScriptManuscript2/MainTables/Lipid percent abundance each class in the seven strains.csv")


# Find the number of different lipid species present in each population at each age

unique(DF$SubClass)
#"TG"   "PC"   "PE"   "CL"   "PE p" "PS"   "LPC"  "PE e" "PI"   "PG"   "PC e" "TG e" "DG"   "DG e" "PS e"

Lipid.species <- DF %>% 
  ungroup() %>% 
  dplyr::filter(Line == "s06" & Age == "19 Day") %>% # for all the other lines
  dplyr::filter(SubClass == "TG e") %>% 
  dplyr::select(Name) %>%
  unique() %>% 
  dplyr::summarise(count = n()) 

# format and save as csv to produce Table S4

DF02 <- DF.main %>% 
  ungroup %>% 
  dplyr::select(Line, Time, Age, SubClass, Name) %>% 
  unique() %>% 
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Name") %>% 
  dplyr::select(Line, Time, Age, SubClass, N) %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = N) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>% 
#write.csv("ScriptManuscript2/MainTables/Lipid diversity each class in the seven strains.csv")


#___________________________Calculate total lipid titre in each strain__________________________________________

DF03 <- DF.main %>% 
  ungroup %>% 
  dplyr::group_by(Samples) %>% 
  mutate(LipidWT = (Conc*Weight*50)*1e-6) %>%  #Convert from ng to mg so LipidWt in mg/fly
  ungroup() %>% 
  dplyr::group_by(Samples) %>% 
  dplyr::mutate(Total.Lipid = sum(LipidWT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(Samples, Line, Age, Time, Weight, Total.Lipid) %>% 
  tidyr::unite(LineTime, Line, Time, remove = FALSE) %>%
  unique() %>% 
  summarySE(groupvars = c("Line", "Time", "Age"), measurevar = "Total.Lipid") %>% 
  dplyr::select(Line, Time, Age, Total.Lipid, se)


# format and save as csv for total lipid titre in mg/fly in Table S3

DF03A <- DF03 %>% 
  dplyr::mutate(Total.Lipid = format(round(.$Total.Lipid, 2), nsmall = 2)) %>% 
  dplyr::mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("TLSE", Total.Lipid:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = TLSE) %>% 
  dplyr::select("ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day") #%>% 
# write.csv("ScriptManuscript2/MainTables/Total lipid titre in mg per fly in the seven strains.csv")


#_____________________Calculate percent lipid abundance in each category for each strain_____________________________

# add another column to DF that has the following 4 lipid categories
PLs <- c("CL","PC", "PE", "PG", "PI", "PS","LPC")
NLs <- c("DG", "TG")
ELPLs <- c("PC e", "PE e", "PE p", "PS e")
ELNLs <- c("DG e", "TG e")

#__________

DF04 <- Complete_DF %>% 
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
  #run the codes till unique:total percentage should come to 100 for each individual sample to check: sum((DF04 %>% filter(Samples == "cbr_n04"))$Percentage)
  dplyr::select(-Total.Conc.by.Type, -Total.Conc.by.SA) %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Percentage") %>% 
  dplyr::select(Line, Time, Age, Type, Percentage, se)

# Check to see if percentage adds to 100 for a Line, Age, Time

DF04 %>% 
  dplyr::filter(Line == "ct" & Time == "New" & Age == "19 Day") %>% 
  dplyr::select(Percentage) %>% 
  sum()

# format and save as csv to produce Table S3
DF04A <- DF04 %>% 
  dplyr::mutate(Percentage = format(round(.$Percentage, 2), nsmall = 2)) %>% 
  dplyr::mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("PeSE", Percentage:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = PeSE) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>% 
#write.csv("ScriptManuscript2/MainTables/Lipid percent abundance each category in the seven strains.csv")


#______________________________Calculate lipid diversity in each category_________________________________________

# format and save as csv to produce Table S3 lipid diversity

DF04B <- DF.main %>% 
  ungroup %>% 
  dplyr::mutate(Type = ifelse(SubClass %in% PLs,"PLs", 
                              ifelse(SubClass %in% ELPLs,"ELPLs", ifelse(SubClass %in% ELNLs,"ELNLs", "NLs")))) %>%  
  dplyr::select(Line, Time, Age, Type, Name) %>% 
  unique() %>% 
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Name") %>% 
  dplyr::select(Line, Time, Age, Type, N) %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = N) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>% 
#write.csv("ScriptManuscript2/MainTables/Lipid diversity each category in the seven strains.csv")


#________________________Calculate PE/PC Ratio _____________________________________

# format and save as csv to produce Table S4 PE/PC ratio

DF05 <- DF.main %>%
  ungroup() %>% 
  dplyr::select(Samples, Line, SubClass, Age, Time, Conc) %>% 
  dplyr::filter(SubClass %in% c("PE", "PC")) %>% 
  dplyr::group_by(Samples, SubClass) %>%
  dplyr::mutate(S.Conc.C = sum(Conc)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(-Conc) %>% 
  unique() %>% 
  tidyr::pivot_wider(values_from = S.Conc.C, names_from = SubClass) %>% 
  dplyr::mutate(Ratio1 = PE/PC) %>%
  unique() %>% 
  dplyr::select(-PC, -PE) %>% 
  summarySE(groupvars = c("Line", "Time", "Age"), measurevar = "Ratio1") %>% 
  dplyr::select(Line, Time, Age, Ratio1, se) %>% 
  mutate(Ratio1 = format(round(.$Ratio1, 1), nsmall = 1)) %>% 
  mutate(se = format(round(.$se, 1), nsmall = 1)) %>% 
  tidyr::unite("PEPC", Ratio1:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = PEPC) %>% 
  dplyr::select("ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
# write.csv("ScriptManuscript2/MainTables/PEPC ratio in the seven strains.csv")



#______________Calculate average number of carbon and double bond in lipid category/class__________________

# First remove subclasses PSe and DGe from DF because these two classes have all lipid species that do not have acyl chain identified

DFnew <- DF.main %>% 
  dplyr::filter(SubClass != "PS e" & SubClass != "DG e")

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
  dplyr::filter(SubClass != "PS e" & SubClass != "DG e")#First remove subclasses PSe and DGe from DF because these two subclasses have all lipid species that do not have acyl chain identified


#___________________________

# Calculate the mean CC present in each lipid category in each of the lipid classes

CategoryCCester <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC1:SC4) %>% 
  dplyr::filter(Type != 'ELNLs' & Type != "ELPLs") %>% 
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
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, SC1CC.Wt:SC4CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1CC.Wt', 'SC2CC.Wt', 'SC3CC.Wt', 'SC4CC.Wt'), names_to='Chain', values_to='ChainCC') %>% 
  dplyr::filter(ChainCC != "NA") %>% # to remove SC3 from DG
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(ChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.Carbon = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.Carbon) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.Carbon", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, Type, Mean.Carbon, se) %>% 
  mutate(CC = format(round(.$Mean.Carbon, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.Carbon", CC:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.Carbon) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/Mean Carbons in ester lipid category in the seven strains.csv")


CategoryDBester <- Composition %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC1:SC4) %>% 
  dplyr::filter(Type != 'ELNLs' & Type != "ELPLs") %>% 
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
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, SC1DB.Wt:SC4DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('SC1DB.Wt', 'SC2DB.Wt', 'SC3DB.Wt', 'SC4DB.Wt'), names_to='Chain', values_to='ChainDB') %>% 
  dplyr::filter(ChainDB != "NA") %>% # to remove SC3 from DG
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(ChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.DB = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.DB) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.DB") %>% 
  dplyr::select(Line, Time, Age, Type, Mean.DB, se) %>% 
  mutate(DB = format(round(.$Mean.DB, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.DB", DB:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.DB) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/Mean DB in ester lipid category in the seven strains.csv")

#____________________

CategoryCCAlkyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "ALkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(AlkylCC.Wt = as.numeric(AlkylCC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, AlkylCC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylCC.Wt'), names_to='AlkylChain', values_to='AlkylChainCC') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.CarbonAlkyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.CarbonAlkyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, Type, Mean.CarbonAlkyl, se) %>% 
  mutate(CCAlkyl = format(round(.$Mean.CarbonAlkyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.CarbonAlkyl", CCAlkyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.CarbonAlkyl) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day")# %>%
#write.csv("ScriptManuscript2/MainTables/Mean Carbons in alkyl chains of lipid category in the seven strains.csv")

CategoryCCAcyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(Acyl1CC.Wt = as.numeric(Acyl1CC)*Conc) %>% 
  dplyr::mutate(Acyl2CC.Wt = as.numeric(Acyl2CC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Acyl1CC.Wt, Acyl2CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1CC.Wt', 'Acyl2CC.Wt'), names_to='AcylChain', values_to='AcylChainCC') %>% 
  dplyr::filter(AcylChainCC != "NA") %>% # to remove SC3 from DG
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AcylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.CarbonAcyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.CarbonAcyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, Type, Mean.CarbonAcyl, se) %>% 
  mutate(CCAcyl = format(round(.$Mean.CarbonAcyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.CarbonAcyl", CCAcyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.CarbonAcyl) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/Mean Carbons in acyl chains of lipid category in the seven strains.csv")


#___________________________

CategoryDBAlkyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(AlkylDB.Wt = as.numeric(AlkylDB)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, AlkylDB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylDB.Wt'), names_to='AlkylChain', values_to='AlkylChainDB') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.DBAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.DBAlkyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.DBAlkyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, Type, Mean.DBAlkyl, se) %>% 
  mutate(DBAlkyl = format(round(.$Mean.DBAlkyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.DBAlkyl", DBAlkyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.DBAlkyl) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/Mean DB in alkyl chains of lipid category in the seven strains.csv")

CategoryDBAcyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, Type) %>% 
  dplyr::mutate(Acyl1DB.Wt = as.numeric(Acyl1DB)*Conc) %>% 
  dplyr::mutate(Acyl2DB.Wt = as.numeric(Acyl2DB)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Type, Line, Time, Age, Acyl1DB.Wt, Acyl2DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1DB.Wt', 'Acyl2DB.Wt'), names_to='AcylChain', values_to='AcylChainDB') %>% 
  dplyr::filter(AcylChainDB != "NA") %>% # to remove SC3 from DG
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, Type) %>% 
  dplyr::mutate(Numerator = sum(AcylChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.DBAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, Type, Mean.DBAcyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "Type"), measurevar = "Mean.DBAcyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, Type, Mean.DBAcyl, se) %>% 
  mutate(DBAcyl = format(round(.$Mean.DBAcyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.DBAcyl", DBAcyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.DBAcyl) %>% 
  dplyr::select(Type, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/Mean DB in acyl chains of lipid category in the seven strains.csv")

#___________________________________________________________________________________

# Calculate the mean CC present in each lipid class

ClassCCester <- Composition %>% 
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
  dplyr::select(Samples, Line, Time, Age, SubClass, Mean.Carbon) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.Carbon") %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.Carbon, se) %>% 
  mutate(Mean.Carbon = format(round(.$Mean.Carbon, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.CC",Mean.Carbon:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.CC) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/CC ester lipid class in the seven strains.csv")

# 3.Variable: Mean double bond content in the ester lipid classes

ClassDBester <- Composition %>% 
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
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.Bond") %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.Bond, se) %>% 
  mutate(Mean.Bond = format(round(.$Mean.Bond, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.DB", Mean.Bond:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.DB) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") #%>%
#write.csv("ScriptManuscript2/MainTables/DB ester lipid class in the seven strains.csv")


# 4.Variable: Mean carbon content in acyl and alkyl chains of ether lipid classes 

ClassCCAlkyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(AlkylCC.Wt = as.numeric(AlkylCC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, AlkylCC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylCC.Wt'), names_to='AlkylChain', values_to='AlkylChainCC') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Mean.CarbonAlkyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.CarbonAlkyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.CarbonAlkyl, se) %>% 
  mutate(Mean.CarbonAlkyl = format(round(.$Mean.CarbonAlkyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("Mean.CCAlkyl", Mean.CarbonAlkyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = Mean.CCAlkyl) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") %>%
  dplyr::mutate(Chain = "Alkyl",.before=SubClass)

ClassCCAcyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(Acyl1CC.Wt = as.numeric(Acyl1CC)*Conc) %>% 
  dplyr::mutate(Acyl2CC.Wt = as.numeric(Acyl2CC)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Acyl1CC.Wt, Acyl2CC.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1CC.Wt', 'Acyl2CC.Wt'), names_to='AcylChain', values_to='AcylChainCC') %>% 
  dplyr::filter(AcylChainCC != "NA") %>% # to remove SC3 from DG
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AcylChainCC)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.CarbonAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Mean.CarbonAcyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.CarbonAcyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.CarbonAcyl, se) %>% 
  mutate(Mean.CarbonAcyl = format(round(.$Mean.CarbonAcyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("CCAcyl", Mean.CarbonAcyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = CCAcyl) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") %>%
  dplyr::mutate(Chain = "Acyl",.before=SubClass)

ClassEtherCC <- rbind(ClassCCAlkyl, ClassCCAcyl) #%>% 
# write.csv("ScriptManuscript2/MainTables/CC ether lipid class in the seven strains.csv")


# 4.Variable: Mean double bond content in acyl and alkyl chains of ether lipids 

ClassDBAlkyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Conc, SC1) %>% 
  tidyr::separate(SC1, into = c("AlkylCC", "AlkylDB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(AlkylDB.Wt = as.numeric(AlkylDB)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, AlkylDB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('AlkylDB.Wt'), names_to='AlkylChain', values_to='AlkylChainDB') %>% 
  dplyr::select(-AlkylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AlkylChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.DBAlkyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Mean.DBAlkyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.DBAlkyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.DBAlkyl, se) %>% 
  mutate(Mean.DBAlkyl = format(round(.$Mean.DBAlkyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("DBAlkyl", Mean.DBAlkyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = DBAlkyl) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") %>%
  dplyr::mutate(Chain = "Alkyl",.before=SubClass)

ClassDBAcyl <- Composition %>% 
  dplyr::filter(Type == 'ELNLs' | Type == "ELPLs") %>% 
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Conc, SC2:SC3) %>% 
  tidyr::separate(SC2, into = c("Acyl1CC", "Acyl1DB"), remove = TRUE) %>% 
  tidyr::separate(SC3, into = c("Acyl2CC", "Acyl2DB"), remove = TRUE) %>% 
  group_by(Samples, Name, SubClass) %>% 
  dplyr::mutate(Acyl1DB.Wt = as.numeric(Acyl1DB)*Conc) %>% 
  dplyr::mutate(Acyl2DB.Wt = as.numeric(Acyl2DB)*Conc) %>% 
  ungroup() %>%
  dplyr::select(Samples, Name, SubClass, Line, Time, Age, Acyl1DB.Wt, Acyl2DB.Wt, Conc) %>% 
  tidyr::pivot_longer(cols=c('Acyl1DB.Wt', 'Acyl2DB.Wt'), names_to='AcylChain', values_to='AcylChainDB') %>% 
  dplyr::filter(AcylChainDB != "NA") %>% # to remove SC3 from DG
  dplyr::select(-AcylChain) %>% 
  group_by(Samples, SubClass) %>% 
  dplyr::mutate(Numerator = sum(AcylChainDB)) %>% 
  dplyr::mutate(Denominator = sum(Conc)) %>% 
  dplyr::mutate(Mean.DBAcyl = Numerator/Denominator) %>% 
  dplyr::ungroup() %>%
  dplyr::select(Samples, Line, Time, Age, SubClass, Mean.DBAcyl) %>% 
  unique() %>%
  summarySE(groupvars = c("Line", "Time", "Age", "SubClass"), measurevar = "Mean.DBAcyl", na.rm = TRUE) %>% 
  dplyr::select(Line, Time, Age, SubClass, Mean.DBAcyl, se) %>% 
  mutate(Mean.DBAcyl = format(round(.$Mean.DBAcyl, 2), nsmall = 2)) %>% 
  mutate(se = format(round(.$se, 2), nsmall = 2)) %>% 
  tidyr::unite("DBAcyl", Mean.DBAcyl:se, sep = "±") %>% 
  tidyr::unite("Strain", Line:Time, sep = "") %>%
  pivot_wider(names_from = c(Strain, Age), values_from = DBAcyl) %>% 
  dplyr::select(SubClass, "ctNew_1 Day","ctNew_19 Day", "ctOld_1 Day", "ctOld_19 Day", "cbrNew_1 Day", "cbrNew_19 Day", 
                "cbrOld_1 Day", "cbrOld_19 Day", "sydNew_1 Day", "sydNew_19 Day", "sydOld_1 Day", "sydOld_19 Day", 
                "s06Older_1 Day", "s06Older_19 Day") %>%
  dplyr::mutate(Chain = "Acyl",.before=SubClass)

ClassEtherDB <- rbind(ClassDBAlkyl,ClassDBAcyl) #%>% 
#write.csv("ScriptManuscript2/MainTables/DB ether lipid class in the seven strains.csv")

#________________________________________________________________________________________________________________________


#__________________________Make correlation plots for CC DB for ELNLs & NLs______________________________________________

rm(list=ls())

# 1. Make a correlation plot for mean CC in alkyl vs acyl chains in ELNLs.

ELNLChainCC <- read_excel("ScriptManuscript2/MainTables/AlkylAcylENL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)") #%>% 
#tidyr::pivot_longer(cols=c(Akyl, Acyl), names_to='Chain', values_to='CC') 


# make correlation plot 

plotA_legend <- ELNLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain carbons ENL", y = "Alkyl chain carbons")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "bottom") +  guides(colour=guide_legend(nrow=1))+
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="Black"))+
  theme(strip.text = element_text(colour = 'white', size = 8))+
  stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)

plotA <- ELNLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain carbons ENL", y = "Alkyl chain carbons")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") +  
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="White"))+
  theme(strip.text = element_text(colour = 'Black', size = 8))+
  #stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)+
  stat_cor(aes(), color = "black", geom = "label", size = 2.5)

#ggsave(plot = plotA, width = 6, height = 4, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_Correlation_plot_ENLCC.jpg")                


#____________________

# 2. Make a correlation plot for mean CC in alkyl vs acyl chains in EPLs.

EPLChainCC <- read_excel("ScriptManuscript2/MainTables/AlkylAcylEPL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") 

plotB <- EPLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = Inf) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain carbons EPL", y = "Alkyl chain carbons")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(face="bold", size = 8))+
  #theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'Black', size = 8))+
  #stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)+
  stat_cor(aes(), color = "black", geom = "label", size = 2.5)


#ggsave(plot = plotB, width = 6, height = 4, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_Correlation_plot_EPLCC.jpg")                


#Arrange and save the correlation plots for CC

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(plotA_legend)

Figure6 <- grid.arrange(plotA, plotB, legend, ncol=2, nrow = 2, layout_matrix = rbind(c(1,2), c(3,3)),
                        widths = c(2.7, 2.7), heights = c(2.5, 0.2))

# Save Figure 6

#ggsave(plot = Figure1, width = 9, height = 6, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Figure6_Alkyl_vs_Acyl_CC_Correlation_plot_ELs.jpg")              

#_______________________________________________________________________________________________________________

#_______________________Correlation plot for mean DB in alkyl vs acyl chains in ELNLs___________________________

# 3. Make a correlation plot for mean DB in alkyl vs acyl chains in ELNLs.

ELNLChainDB <- read_excel("ScriptManuscript2/MainTables/AlkylAcylENL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)") 


# make correlation plot 

plotC <- ELNLChainDB %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain double bonds ENL", y = "Alkyl chain double bonds")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="White"))+
  theme(strip.text = element_text(colour = 'Black', size = 8))+
  #stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)
  stat_cor(aes(), color = "black", geom = "label", size = 2.5)


#ggsave(plot = plotC, width = 6, height = 4, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_Correlation_plot_ENLDB.jpg")                

#____________________

# 4. Make a correlation plot for mean DB in alkyl vs acyl chains in ELPLs.

ELPLChainDB <- read_excel("ScriptManuscript2/MainTables/AlkylAcylEPL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") 


# make correlation plot 

plotD <- ELPLChainDB %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = Inf) +
  #stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain double bonds EPL", y = "Alkyl chain double bonds")+
  theme_bw()+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #facet_grid(~ Age)+
  facet_grid(Age ~ .)+
  theme(strip.background =element_rect(fill="White"))+
  theme(strip.text = element_text(colour = 'Black', size = 8))+
  #stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)
  stat_cor(aes(), color = "black", geom = "label", size = 2.5)

#ggsave(plot = plotD, width = 6, height = 4, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_Correlation_plot_EPLDB.jpg")              

FigureS2 <- grid.arrange(plotC, plotD, legend, ncol=2, nrow = 2, layout_matrix = rbind(c(1,2), c(3,3)),
                         widths = c(2.7, 2.7), heights = c(2.5, 0.2))

# Save Figure S2

ggsave(plot = FigureS2, width =9, height = 6, units = "in", dpi = 300,filename = "Figures/FigureS2.jpg")              

#________________________________________END_________________________________________
