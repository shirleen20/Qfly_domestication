#_________________________________________24/01/2025____________________________________

# Here we make correlation plots for carbon contents in alkyl and acyl chains of ether lipids
# We use this analysis to generate Figure 7

library("tidyverse")
library("ggrepel")
library("ggpubr")
library("readxl")
library("broom")
library("gtools")
library("multcomp")
library("multcompView")
library("Rmisc")
library("gridExtra")
library("grid")

rm(list=ls())

#_______________________________Correlation plots for CC at Day 1______________________________________________


# 1. Make a correlation plot for mean CC in alkyl vs acyl chains in ELNLs.

ELNLChainCC <- read_excel("data/data_plot/AlkylAcylENL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)") %>% 
  dplyr::filter(Age == "Day 1")

# make correlation plot 

plot_legend <- ELNLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  ggtitle("Ether neutral lipids") +
  theme_bw()+
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "bottom") +  guides(colour=guide_legend(nrow=1))+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #scale_y_continuous(limits = c(17.5, 20.5, by = 0.5))+
  scale_y_continuous(limits = c(17.5, 20.5, by = 0.5),labels = scales::number_format(accuracy = 0.1))+
  stat_cor(aes(), color = "black", geom = "label", size = 2.5)

#________________________DAY 1 PLOTS Acyl vs Alk(en)yl CC_____________________________


# Make plots for D1

plotENL <- ELNLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  ggtitle("Ether neutral lipids") +
  theme_bw()+
  scale_y_continuous(limits = c(18.1, 20.1, by = 0.5),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 20.1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#____________________

# 2. Make a correlation plot for mean CC in alkyl vs acyl chains in EPLs.

EPLChainCC <- read_excel("data/data_plot/AlkylAcylEPL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") %>% 
  dplyr::filter(Age == "Day 1")

plotEPL <- EPLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = Inf) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  ggtitle("Ether phospholipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.25))+
  scale_y_continuous(limits = c(14.1, 14.4, by = 0.2),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 14.4) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


#____________________________Correlation plots for CC for EPL classes____________________________

# make correlation plot for PEe

CC_PEe <- read_excel("data/data_plot/CC_EPLClasses.xlsx") %>% 
  dplyr::select(1,2,3,6,7,12,13) %>% 
  dplyr::rename(Alkyl = "PEe(alkyl)") %>% 
  dplyr::rename(Acyl = "PEe(acyl)") %>% 
  dplyr::rename(alkylSE = "PEe_se(alkyl)") %>% 
  dplyr::rename(acylSE = "PEe_se(acyl)") %>% 
  dplyr::filter(Age == "1 Day")

plotPEe <- CC_PEe %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-alkylSE, ymax=Alkyl+alkylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  ggtitle("PEe") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(limits = c(15.2, 15.8, by = 0.1),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") +  guides(colour=guide_legend(nrow=1)) + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 15.8) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# 3.make correlation plot for PEp

CC_PEp <- read_excel("data/data_plot/CC_EPLClasses.xlsx") %>% 
  dplyr::select(1,2,3,8,9,14,15) %>% 
  dplyr::rename(Alkyl = "PEp(alkyl)") %>% 
  dplyr::rename(Acyl = "PEp(acyl)") %>% 
  dplyr::rename(alkylSE = "PEp_se(alkyl)") %>% 
  dplyr::rename(acylSE = "PEp_se(acyl)")  %>% 
  dplyr::filter(Age == "1 Day")

plotPEp <- CC_PEp %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-alkylSE, ymax=Alkyl+alkylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  ggtitle("PEp") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(limits = c(13.15, 13.5, by = 0.1),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) +  
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 13.5) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


#Arrange and save the correlation plots for CC

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(plot_legend)

#___
#Arrange and save D1 plots 

Plot1 <- ggarrange(plotENL,plotEPL, plotPEe, plotPEp, ncol = 4)

Plot1_annotated <- annotate_figure(Plot1, top = text_grob("Day 1", hjust = 0.5, face = "bold", size = 14))
#                bottom = text_grob("Acyl chain(s)", hjust = 0.5, face = "bold", size = 12),
#                left = text_grob("Alk(en)yl chain", rot = 90, face = "bold", size = 12))

#_______________________________________________________________________________________________________________


#rm(list=ls())

#_______________________________Correlation plots for CC at Day 19______________________________________________


# 1. Make a correlation plot for mean CC in alkyl vs acyl chains in ELNLs.

ELNLChainCC <- read_excel("data/data_plot/AlkylAcylENL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)") %>% 
  dplyr::filter(Age == "Day 19")


#________________________DAY 1 PLOTS Acyl vs Alk(en)yl CC_____________________________


# Make plots for D19

plotENL <- ELNLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  #ggtitle("Ether neutral lipids") +
  theme_bw()+
  scale_y_continuous(limits = c(17.0, 20.0, by = 1),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 20.0) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#____________________

# 2. Make a correlation plot for mean CC in alkyl vs acyl chains in EPLs.

EPLChainCC <- read_excel("data/data_plot/AlkylAcylEPL_CC.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") %>% 
  dplyr::filter(Age == "Day 19")

plotEPL <- EPLChainCC %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = Inf) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  #ggtitle("Ether phospholipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.25))+
  scale_y_continuous(limits = c(13.8, 14.1, by = 0.1),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 14.1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


#____________________________Correlation plots for CC for EPL classes____________________________

# make correlation plot for PEe

CC_PEe <- read_excel("data/data_plot/CC_EPLClasses.xlsx") %>% 
  dplyr::select(1,2,3,6,7,12,13) %>% 
  dplyr::rename(Alkyl = "PEe(alkyl)") %>% 
  dplyr::rename(Acyl = "PEe(acyl)") %>% 
  dplyr::rename(alkylSE = "PEe_se(alkyl)") %>% 
  dplyr::rename(acylSE = "PEe_se(acyl)") %>% 
  dplyr::filter(Age == "19 Day")

plotPEe <- CC_PEe %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-alkylSE, ymax=Alkyl+alkylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  #ggtitle("PEe") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(limits = c(15.9, 16.6, by = 0.2),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") +  guides(colour=guide_legend(nrow=1)) + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 16.6) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# 3.make correlation plot for PEp

CC_PEp <- read_excel("data/data_plot/CC_EPLClasses.xlsx") %>% 
  dplyr::select(1,2,3,8,9,14,15) %>% 
  dplyr::rename(Alkyl = "PEp(alkyl)") %>% 
  dplyr::rename(Acyl = "PEp(acyl)") %>% 
  dplyr::rename(alkylSE = "PEp_se(alkyl)") %>% 
  dplyr::rename(acylSE = "PEp_se(acyl)")  %>% 
  dplyr::filter(Age == "19 Day")

plotPEp <- CC_PEp %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-alkylSE, ymax=Alkyl+alkylSE,colour = Strain), width = 0.001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain(s)", y = "Alk(en)yl chain")+
  #ggtitle("PEp") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  scale_y_continuous(limits = c(12.7, 13.0, by = 0.5),labels = scales::number_format(accuracy = 0.1))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + guides(colour=guide_legend(nrow=1)) +  
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 13.0) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#___
#Arrange and save D19 plots 

Plot2 <- ggarrange(plotENL,plotEPL, plotPEe, plotPEp, ncol = 4)

Plot2_annotated <- annotate_figure(Plot2, top = text_grob("Day 19", hjust = 0.5, face = "bold", size = 14))

#___
# Combine and save plot for D1 and D19 into one Figure

Figure <- ggarrange(Plot1_annotated, Plot2_annotated, nrow = 2)

Figure_annotated <- annotate_figure(Figure,
                                    #                                  top = text_grob("Day 1", hjust = 0.5, face = "bold", size = 14))
                                    bottom = text_grob("Acyl chain(s)", hjust = 0.5, face = "bold", size = 12),
                                    left = text_grob("Alk(en)yl chain", rot = 90, face = "bold", size = 12))


Figure_arranged <- grid.arrange(Figure_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                widths = c(2.7, 2.7), heights = c(2.5, 0.2))



#ggsave(plot = Figure_arranged, width = 9.0, height = 6.0, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_CC_Correlation_plot_ELs.jpg")              

#________________________________________END________________________________________________