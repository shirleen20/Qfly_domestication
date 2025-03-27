#_________________________________________24/01/2025____________________________________

# Heatmaps of frequencies (%) of abundant acyl and alk(en)yl chains in two strains CTnew and 
# SDolder and classes TG, PE, TGe, PEe, PEp at Day 1 and Day 19
# We use this analysis to generate Figures 6, S2 and S3.

library("tidyverse")
library("ggplot2")
library("readxl")
library("Rmisc")
library("ggpubr")
library("gridExtra")
library("grid")

#____________________________Plot heat mapsfor Age D1____________________________________


# 1. CTnew D1 heatmap

rm(list=ls())

#__________

CTnewAbundance <- as.data.frame(read_csv("data/data_plot/CTnew_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = ct_New) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))

#CTnewAbundance %>% dplyr::select(AcylChain) %>% unique()

p1 <- CTnewAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CTnew") +
  labs(x= "Lipid classes CTnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

# Save a plot with legend for ether lipid

p1.Legend <- CTnewAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CTnew") +
  labs(x= "Lipid classes CTnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())


# 2. CTold D1 heatmap

CToldAbundance <- as.data.frame(read_csv("data/data_plot/CTold_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = ct_Old) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))


p2 <- CToldAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CTold") +
  labs(x= "Lipid classes CTold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())


# 3. SDnew D1 heatmap

SDnewAbundance <- as.data.frame(read_csv("data/data_plot/SDnew_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = syd_New) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))


p3 <- SDnewAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDnew") +
  labs(x= "Lipid classes SDnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())


# 4. SDold D1 heatmap

SDoldAbundance <- as.data.frame(read_csv("data/data_plot/SDold_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = syd_Old) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))


p4 <- SDoldAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  #theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDold") +
  labs(x= "Lipid classes SDold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


# 5. SDolder D1 heatmap

SDolderAbundance <- as.data.frame(read_csv("data/data_plot/SD0lder_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = s06_Older) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))

#SD0lderAbundance1 %>% dplyr::select(AcylChain) %>% unique()

p5 <- SDolderAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDolder") +
  labs(x= "Lipid subclasses SDolder", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


# 6. CNnew D1 heatmap

CNnewAbundance <- as.data.frame(read_csv("data/data_plot/CNnew_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = cbr_New) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))


p6 <- CNnewAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CNnew") +
  labs(x= "Lipid classes CNnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())



# 7. CNold D1 heatmap

CNoldAbundance <- as.data.frame(read_csv("data/data_plot/CNold_commonacylchains.csv")) %>% 
  dplyr::mutate(Age = paste0(Age, " Day")) %>% 
  tidyr::unite(AcylChain, Age, Class, remove = TRUE) %>% 
  dplyr::rename(PercentageAbundance = cbr_Old) %>% 
  dplyr::mutate(SideChain = gsub(":","\\.", SideChain) %>% as.numeric()) %>%  # Substitute ":" with "." and treat as numeric so that we can treat chains as number and arrange them in descending order, those smaller than 10 were appearing later when we treated them as character for e.g. 8:0 appeared after 10:0
  arrange(desc(SideChain)) %>%  # arrange in descending order
  dplyr::mutate(SideChain = format(SideChain, nsmall = 1) %>% gsub(" ","",.) %>% gsub("\\.",":",.) ) %>%
  dplyr::mutate(SideChain = as.factor(SideChain)) %>%
  dplyr::mutate(SideChain = factor(SideChain, levels = unique(SideChain) %>%
                                     gsub(":","\\.", .) %>%
                                     as.numeric() %>%
                                     sort(decreasing = TRUE) %>%
                                     format(nsmall =1) %>%
                                     gsub(" ","",.) %>%
                                     gsub("\\.",":",.)))


p7 <- CNoldAbundance %>%
  dplyr::filter(AcylChain %in% c("1 Day_TG","1 Day_TG e", "1 Day_TG e_ether", 
                                 "1 Day_PE", "1 Day_PE e", "1 Day_PE e_ether", 
                                 "1 Day_PE p", "1 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "1 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>% 
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CNold") +
  labs(x= "Lipid classes CNold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

#__________________

# First arrange and save D1 plots only for CTnew and SDolder to make Figure 6

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1.Legend)

Plot1 <- ggarrange(p1,p5, nrow = 2)

Plot1_annotated <- annotate_figure(Plot1, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

Figure1_arranged <- grid.arrange(Plot1_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                 widths = c(2.7, 2.7), heights = c(2.7, 0.2))

#ggsave(plot = Figure1_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "ScriptManuscript2/Commonchains/Common acyl chains CTnew & SDolder at D1.jpg")
ggsave(plot = Figure1_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "Figures/Figure 6.jpg")

#_______

# Second arrange and save D1 plots for all strains to make Figure S2

plot2.1 <- ggarrange(p1,p2, ncol = 2)

plot2.2 <- ggarrange(p6,p7, ncol = 2)

plot2.3 <- ggarrange(p3,p4, ncol = 2)

plot2.4 <- ggarrange(p5, legend, ncol = 2)

plot2 <- ggarrange(plot2.1, plot2.2, plot2.3, plot2.4, nrow = 4)


Plot2_annotated <- annotate_figure(plot2, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

#Accepted ms dimensions
#ggsave(plot = Plot2_annotated, width = 6.6, height = 8.8, units = "in", dpi = 300,filename = "Figures/Figure S2.jpg")
ggsave(plot = Plot2_annotated, width = 14, height = 12, units = "in", dpi = 300,filename = "Figures/Figure S2.jpg")


#_____________________________________________________________________________________________


#__________________________D19 HEATMAP ALL STRAINS____________________________________________

P1 <- CTnewAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CTnew") +
  labs(x= "Lipid classes CTnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

P2 <- CToldAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5)) +
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CTold") +
  labs(x= "Lipid classes CTold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

P3 <- SDnewAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDnew") +
  labs(x= "Lipid subclasses SDnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

P4 <- SDoldAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDold") +
  labs(x= "Lipid subclasses SDold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


P5 <- SDolderAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("SDolder") +
  labs(x= "Lipid subclasses SDolder", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

P6 <- CNnewAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CNnew") +
  labs(x= "Lipid subclasses CNnew", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

P7 <- CNoldAbundance %>%
  dplyr::filter(AcylChain %in% c("19 Day_TG","19 Day_TG e", "19 Day_TG e_ether", 
                                 "19 Day_PE", "19 Day_PE e", "19 Day_PE e_ether", 
                                 "19 Day_PE p", "19 Day_PE p_vinylether"))%>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG", "TG(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e", "TGe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_TG e_ether", "TGe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE", "PE(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e", "PEe(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE e_ether", "PEe(alk)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p", "PEp(acy)")) %>% 
  dplyr::mutate(AcylChain = replace(AcylChain, AcylChain == "19 Day_PE p_vinylether", "PEp(alk)")) %>% 
  dplyr::mutate(AcylChain = factor(as.factor(AcylChain), levels = c("TG(acy)","TGe(acy)", "TGe(alk)", 
                                                                    "PE(acy)", "PEe(acy)", "PEe(alk)", 
                                                                    "PEp(acy)", "PEp(alk)"))) %>%
  ggplot(aes(x = AcylChain, y = SideChain, fill = PercentageAbundance)) + # untransformed percentageabundance
  geom_tile() +
  #theme_bw()+
  scale_fill_gradient(low = "white",high = "darkblue",guide = "colorbar")+ 
  ggtitle("CNold") +
  labs(x= "Lipid subclasses CNold", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

#_______________

# First arrange and save D19 plots only for CTnew and SDolder 

Plot3 <- ggarrange(P1,P5, nrow = 2)

Plot3_annotated <- annotate_figure(Plot3, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

Figure3_arranged <- grid.arrange(Plot3_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                 widths = c(2.7, 2.7), heights = c(2.7, 0.2))


#ggsave(plot = Figure2_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "Figures/Figure S3.jpg")

#_______

# Second arrange and save D1 plots for all strains to make Figure S2

plot4.1 <- ggarrange(P1,P2, ncol = 2)

plot4.2 <- ggarrange(P6,P7, ncol = 2)

plot4.3 <- ggarrange(P3,P4, ncol = 2)

plot4.4 <- ggarrange(P5, legend, ncol = 2)

plot4 <- ggarrange(plot4.1, plot4.2, plot4.3, plot4.4, nrow = 4)


Plot4_annotated <- annotate_figure(plot4, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

#Accepted ms dimensions
#ggsave(plot = Plot4_annotated, width = 6.6, height = 8.8, units = "in", dpi = 300,filename = "Figures/Figure S3.jpg")

ggsave(plot = Plot4_annotated, width = 14, height = 12, units = "in", dpi = 300,filename = "Figures/Figure S3.jpg")


#______________________________________END_______________________________________


