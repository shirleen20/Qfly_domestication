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

#__________ SDolder D1 heatmap_________

SD0lderAbundance <- as.data.frame(read_csv("data/data_plot/SD0lder_commonacylchains.csv")) %>% 
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

p2 <- SD0lderAbundance %>%
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
  labs(x= "Lipid subclasses SDOlder", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())



#___
#Arrange and save D1 plots only for CTnew and SDolder to make Figure 6

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1.Legend)

Plot1 <- ggarrange(p1,p2, nrow = 2)

Plot1_annotated <- annotate_figure(Plot1, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

Figure1_arranged <- grid.arrange(Plot1_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                 widths = c(2.7, 2.7), heights = c(2.7, 0.2))

#ggsave(plot = Figure1_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "ScriptManuscript2/Commonchains/Common acyl chains CTnew & SDolder at D1.jpg")
ggsave(plot = Figure1_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "Figures/Figure 7.jpg")


#_____________________________________________________________________________________________

#__________________________D19 HEATMAP ALL STRAINS____________________________________________

pp3 <- CTnewAbundance %>%
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

p4 <- SD0lderAbundance %>%
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
  labs(x= "Lipid subclasses SDOlder", y = "Acyl/alk(en)yl chains", fill = "Percentage abundance") + 
  theme(axis.title.y = element_text(face = "bold"), axis.title.x = element_text(face="bold")) + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12)) + 
  theme(legend.position="none") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = (12))) + theme(axis.title = element_text(size = 8))+
  #theme(axis.text.x=element_blank()
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

Plot2 <- ggarrange(p3,p4, nrow = 2)

Plot2_annotated <- annotate_figure(Plot2, bottom = text_grob("Lipid classes", hjust = 0.5, face = "bold", size = 12),
                                   left = text_grob("Acyl/Alk(en)yl chains", rot = 90, face = "bold", size = 12))

Figure2_arranged <- grid.arrange(Plot2_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                 widths = c(2.7, 2.7), heights = c(2.7, 0.2))

#ggsave(plot = Figure2_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "Figures/Figure S3 Common acyl chains CTnew & SDolder at D19.jpg")
ggsave(plot = Figure2_arranged, width = 9.0, height = 8.5, units = "in", dpi = 300,filename = "Figures/Figure S3.jpg")

#______________________________________END_______________________________________

