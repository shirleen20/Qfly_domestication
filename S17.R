#_________________________________________24/01/2025____________________________________

# Here we make correlation plots for double bond contents in alkyl and acyl chains of ether lipids
# We use this analysis to generate Figure S2

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

#_______________________________Correlation plots for DB at Day 1______________________________________________

#Make a correlation plot for mean DB in alkyl vs acyl chains in ELNLs.

ELNLChainDB1 <- read_excel("data/data_plot/AlkylAcylENL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)")  %>% 
  dplyr::filter(Age == "Day 1")

plotELNLDB1 <- ELNLChainDB1 %>% 
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
  ggtitle("Ether neutral lipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  #stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", size = 3)
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 0.635)+
  scale_y_continuous(limits = c(0.00, 0.635, by = 0.20),labels = scales::number_format(accuracy = 0.001))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#____________________

# Make a correlation plot for mean DB in alkyl vs acyl chains in ELPLs.

ELPLChainDB1 <- read_excel("data/data_plot/AlkylAcylEPL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") %>% 
  dplyr::filter(Age == "Day 1") 


# make correlation plot 

plotELPLDB1 <- ELPLChainDB1 %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.000001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.000001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain double bonds EPL", y = "Alk(en)yl chain double bonds")+
  ggtitle("Ether phosoholipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 0.042, label.x = 1.726)+
  scale_y_continuous(limits = c(0.020, 0.042, by = 0.025),labels = scales::number_format(accuracy = 0.001))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


legend <- ELPLChainDB1 %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.000001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width =0.000001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5, max.overlaps = Inf) +
  labs(x= "Acyl chain double bonds EPL", y = "Alkyl chain double bonds")+
  ggtitle("Ether phospholipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "bottom") + guides(colour=guide_legend(nrow=1))+
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 0.042)+
  scale_y_continuous(limits = c(0.02, 0.042, by = 0.005),labels = scales::number_format(accuracy = 0.001))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())


#Arrange and save the correlation plots for DB

get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(legend)

#___
#Arrange and save D1 plots 

Plot1 <- ggarrange(plotELNLDB1,plotELPLDB1, ncol = 2)

Plot1_annotated <- annotate_figure(Plot1,top = text_grob("Day 1", hjust = 0.5, face = "bold", size = 14))

#______________________________________________________________________________________________________________

#_______________________________Correlation plots for CC at Day 19______________________________________________

# Make a correlation plot for mean DB in alkyl vs acyl chains in ELNLs.

ELNLChainDB19 <- read_excel("data/data_plot/AlkylAcylENL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "ENLs (akyl)") %>% 
  dplyr::rename(Acyl = "ENLs (acyl)")  %>% 
  dplyr::filter(Age == "Day 19")

plotELNLDB19 <- ELNLChainDB19 %>% 
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
  ggtitle("Ether neutral lipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 1.500)+
  scale_y_continuous(limits = c(0.00, 1.500, by = 0.25),labels = scales::number_format(accuracy = 0.001))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#____________________

# Make a correlation plot for mean DB in alkyl vs acyl chains in ELPLs.

ELPLChainDB19 <- read_excel("data/data_plot/AlkylAcylEPL_DB.xlsx") %>% 
  dplyr::rename(Alkyl = "EPLs (akyl)") %>% 
  dplyr::rename(Acyl = "EPLs (acyl)") %>% 
  dplyr::filter(Age == "Day 19") 


# make correlation plot 

plotELPLDB19 <- ELPLChainDB19 %>% 
  dplyr::mutate(Alkyl = as.numeric(Alkyl)) %>% 
  dplyr::mutate(Acyl = as.numeric(Acyl)) %>% 
  dplyr::mutate(Strain = factor(as.factor(Strain), levels = c("CTnew","CTold", "SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(Acyl, Alkyl, colour = Strain))+
  geom_point()+ 
  geom_errorbar(aes(xmin=Acyl-acylSE, xmax=Acyl+acylSE,colour = Strain), width = 0.000001, size = 0.2, alpha = 1) +
  geom_errorbar(aes(ymin=Alkyl-akylSE, ymax=Alkyl+akylSE,colour = Strain), width = 0.000001, size = 0.2, alpha = 1) +
  geom_text_repel(aes(label = paste0(Strain)), size = 3, min.segment.length = 0, seed = 42, box.padding = 0.5) +
  stat_smooth(method = "lm",col = "#C42126", se = FALSE,size = 0.5)+
  labs(x= "Acyl chain double bonds EPL", y = "Alk(en)yl chain double bonds")+
  ggtitle("Ether phosoholipids") +
  theme_bw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  theme(axis.title.y = element_text(face = "bold", size = 8), axis.title.x = element_text(face="bold", size = 8))+
  theme(axis.text.y = element_text(face = "bold", size = 8), axis.text.x = element_text(face="bold", size = 8))+
  theme(legend.position = "none") + 
  theme(legend.title=element_text(size= 8), legend.text = element_text(size = (8))) + 
  theme(axis.title = element_text(size = 10)) +
  stat_cor(aes(), color = "black", geom = "label", size = 2.5, label.y = 0.080)+
  scale_y_continuous(limits = c(0.060, 0.080, by = 0.005),labels = scales::number_format(accuracy = 0.001))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#___
#Arrange and save D19 plots 

Plot2 <- ggarrange(plotELNLDB19,plotELPLDB19, ncol = 2)

Plot2_annotated <- annotate_figure(Plot2, top = text_grob("Day 19", hjust = 0.5, face = "bold", size = 14))

#___
# Combine and save plot for D1 and D19 into one Figure

Figure <- ggarrange(Plot1_annotated, Plot2_annotated, nrow = 2)

Figure_annotated <- annotate_figure(Figure,
                                    #top = text_grob("Day 1", hjust = 0.5, face = "bold", size = 14))
                                    bottom = text_grob("Acyl chain(s)", hjust = 0.5, face = "bold", size = 12),
                                    left = text_grob("Alk(en)yl chain", rot = 90, face = "bold", size = 12))


Figure_arranged <- grid.arrange(Figure_annotated, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                                widths = c(2.7, 2.7), heights = c(2.5, 0.2))


#ggsave(plot = Figure_arranged, width = 9.0, height = 6.0, units = "in", dpi = 300,filename = "ScriptManuscript2/Figures/Alkyl_vs_Acyl_DB_Correlation_plot_ELs.jpg")              

#________________________________________END_________________________________________


