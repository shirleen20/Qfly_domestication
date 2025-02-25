#______________________________________24/01/2025_______________________________

# We make Figure 1 and Figure S1

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

rm(list=ls())

#______________ Create main data frame ________________________________________________________________

Variables <- read_excel("data/Data_F1.xlsx")

pd = position_dodge(0.4)

#_____________________________Figure S1__________________________

p <- Variables %>% 
  dplyr::select(Strains, Age, Wt, Wt.Se) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = Wt)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 2, position = pd) +
  geom_errorbar(aes(ymin = Wt-Wt.Se, ymax = Wt+Wt.Se, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  #theme(strip.background =element_rect(fill="white"))+
  #theme(strip.text = element_text(colour = 'Black', size =10))+
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)) +
  labs(x= "Strains", y = "Fly weight")+
  theme(legend.position = "bottom") +  guides(colour=guide_legend(nrow=1))+ 
  theme(legend.text=element_text(size=8)) +
  #theme(legend.title=element_text(size= 7), legend.text = element_text(size = (7))) +
  theme(legend.title = element_blank())+ 
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("Fly weight")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  #theme(axis.text.y = element_text(face = "bold", size = 7))+
  theme(axis.text.y = element_text(size = 10))+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(breaks = seq(1.5, 5.5, by=1), limits=c(1.5,5.5)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 

ggsave(plot = p, width = 8.0, height = 4.0, units = "in", dpi = 300,filename = "Figures/FigureS1.jpg")              

#_____________________________Figure 1__________________________

p1.legend <- Variables %>% 
  dplyr::select(Strains, Age, Pt, Pt.Se) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = Pt)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = Pt-Pt.Se, ymax = Pt+Pt.Se, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        #axis.text    = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x= "Strains", y = "Lipid weight (%)")+
  theme(legend.position = "bottom") +  guides(colour=guide_legend(nrow=1))+ 
  theme(legend.text=element_text(size=7)) +
  #theme(legend.title=element_text(size= 7), legend.text = element_text(size = (7))) +
  theme(legend.title = element_blank())+ 
  theme(axis.title.x =element_blank()) + theme(axis.text.x =element_blank()) +
  theme(axis.title.y = element_text(face = "bold", size = 7))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 

p1 <- Variables %>% 
  dplyr::select(Strains, Age, Pt, Pt.Se) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = Pt)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = Pt-Pt.Se, ymax = Pt+Pt.Se, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) + 
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  #theme(strip.background =element_rect(fill="white"))+
  #theme(strip.text = element_text(colour = 'Black', size =10))+
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "Lipid weight (%)")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("Total lipids (%)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  #theme(axis.text.y = element_text(face = "bold", size = 7))+
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.title.y = element_blank())+
  #scale_y_continuous(breaks = seq(1250, 2280, by=250), limits=c(1250, 2280)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 

p2 <- Variables %>% 
  dplyr::select(Strains, Age, TL, TL.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = TL)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = TL-TL.SE, ymax = TL+TL.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) + 
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  #theme(strip.background =element_rect(fill="white"))+
  #theme(strip.text = element_text(colour = 'Black', size =10))+
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "Total lipids")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("Total lipids (mg)")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  #theme(axis.text.y = element_text(face = "bold", size = 7))+
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.title.y = element_blank())+
  #scale_y_continuous(breaks = seq(1250, 2280, by=250), limits=c(1250, 2280)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 

#________________________________Plots for abundance_____________________________________

# We need these functions to plot p3

max_first1  <- 44.3   # Specify max of first y axis
max_second1 <- max(Variables$NLs.d) # Specify max of second y axis
min_first1  <- 0   # Specify min of first y axis
min_second1 <- min(Variables$NLs.d) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale1 = (max_second1 - min_second1)/(max_first1 - min_first1)
shift1 = min_first1 - min_second1

# Function to scale secondary axis
scale_function1 <- function(x, scale1, shift1){
  return ((x)*scale1 - shift1)
}

# Function to scale secondary variable values
inv_scale_function1 <- function(x, scale1, shift1){
  return ((x + shift1)/scale1)
}

p3 <- Variables %>% 
  dplyr::select(Strains, Age, NLs.d, NLs.ab, NLs.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, colour = Strains)) +
  geom_col(mapping = aes(y=NLs.ab), fill = NA, size = 0.09) +
  geom_errorbar(aes(ymin = NLs.ab-NLs.SE, ymax = NLs.ab+NLs.SE), width = 0.25, size = 0.5, alpha = 0.4) +
  geom_line(aes(y=inv_scale_function1(NLs.d, scale1, shift1), group=1, colour = "black"), size = 0.2) +
  geom_point(mapping = aes(y=inv_scale_function1(NLs.d, scale1, shift1), colour = Strains, shape  = Strains), size = 0.2, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(min_first1, max_first1), sec.axis = sec_axis(~scale_function1(., scale1, shift1), name="Diversity")) +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  #coord_cartesian(expand = FALSE) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "NL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("NL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15))

#______________________________________________
# We need these functions to plot p4

max_first2  <- 80   # Specify max of first y axis
max_second2 <- max(Variables$PLs.d) # Specify max of second y axis
min_first2  <- 0   # Specify min of first y axis
min_second2 <- min(Variables$PLs.d) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale2 = (max_second2 - min_second2)/(max_first2 - min_first2)
shift2 = min_first2 - min_second2

# Function to scale secondary axis
scale_function2 <- function(x, scale2, shift2){
  return ((x)*scale2 - shift2)
}

# Function to scale secondary variable values
inv_scale_function2 <- function(x, scale2, shift2){
  return ((x + shift2)/scale2)}

p4 <- Variables %>% 
  dplyr::select(Strains, Age, PLs.d, PLs.ab, PLs.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, colour = Strains)) +
  geom_col(mapping = aes(y=PLs.ab), fill = NA, size = 0.09) +
  geom_errorbar(aes(ymin = PLs.ab-PLs.SE, ymax = PLs.ab+PLs.SE), width = 0.25, size = 0.5, alpha = 0.4) +
  geom_line(aes(y=inv_scale_function2(PLs.d, scale2, shift2), group=1, colour = "black"), size = 0.2) +
  geom_point(mapping = aes(y=inv_scale_function2(PLs.d, scale2, shift2), colour = Strains, shape  = Strains), size = 0.2, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(min_first2, max_first2), sec.axis = sec_axis(~scale_function2(., scale2, shift2), name="Diversity")) +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "PL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("PL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15))

#______________________________________________
# We need these functions to plot p5

max_first3  <- 1.7   # Specify max of first y axis
max_second3 <- max(Variables$ELNLs.d) # Specify max of second y axis
min_first3  <- 0   # Specify min of first y axis
min_second3 <- min(Variables$ELNLs.d) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale3 = (max_second3 - min_second3)/(max_first3 - min_first3)
shift3 = min_first3 - min_second3

# Function to scale secondary axis
scale_function3 <- function(x, scale3, shift3){
  return ((x)*scale3 - shift3)
}

# Function to scale secondary variable values
inv_scale_function3 <- function(x, scale3, shift3){
  return ((x + shift3)/scale3)}

p5 <- Variables %>% 
  dplyr::select(Strains, Age, ELNLs.d, ELNLs.ab, ELNLs.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, colour = Strains)) +
  geom_col(mapping = aes(y=ELNLs.ab), fill = NA, size = 0.09) +
  geom_errorbar(aes(ymin = ELNLs.ab-ELNLs.SE, ymax = ELNLs.ab+ELNLs.SE), width = 0.25, size = 0.5, alpha = 0.4) +
  geom_line(aes(y=inv_scale_function3(ELNLs.d, scale3, shift3), group=1, colour = "black"), size = 0.2) +
  geom_point(mapping = aes(y=inv_scale_function3(ELNLs.d, scale3, shift3), colour = Strains, shape  = Strains), size = 0.2, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(min_first3, max_first3), sec.axis = sec_axis(~scale_function3(., scale3, shift3), name="Diversity", labels = scales::number_format(accuracy = 0.1))) +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELNL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELNL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(face = "bold", size = 7))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15))


# We need these functions to plot p6

max_first4  <- 6.0   # Specify max of first y axis
max_second4 <- max(Variables$ELPLs.d) # Specify max of second y axis
min_first4  <- 0   # Specify min of first y axis
min_second4 <- min(Variables$ELPLs.d) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale4 = (max_second4 - min_second4)/(max_first4 - min_first4)
shift4 = min_first4 - min_second4

# Function to scale secondary axis
scale_function4 <- function(x, scale4, shift4){
  return ((x)*scale4 - shift4)
}

# Function to scale secondary variable values
inv_scale_function4 <- function(x, scale4, shift4){
  return ((x + shift4)/scale4)}

p6 <- Variables %>% 
  dplyr::select(Strains, Age, ELPLs.d, ELPLs.ab, ELPLs.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, colour = Strains)) +
  geom_col(mapping = aes(y=ELPLs.ab), fill = NA, size = 0.09) +
  geom_errorbar(aes(ymin = ELPLs.ab-ELPLs.SE, ymax = ELPLs.ab+ELPLs.SE), width = 0.25, size = 0.5, alpha = 0.4) +
  geom_line(aes(y=inv_scale_function4(ELPLs.d, scale4, shift4), group=1, colour = "black"), size = 0.2) +
  geom_point(mapping = aes(y=inv_scale_function4(ELPLs.d, scale4, shift4), colour = Strains, shape  = Strains), size = 0.2, position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(min_first4, max_first4), sec.axis = sec_axis(~scale_function4(., scale4, shift4), name="Diversity")) +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +     
  labs(x= "Strains", y = "ELPL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELPL Ab.")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15))


#_____________________________________________Plots for carbon contents_________________________________________________________

p7 <- Variables %>% 
  dplyr::select(Strains, Age, NLs.CC, NLs.CC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = NLs.CC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = NLs.CC-NLs.CC.SE, ymax = NLs.CC+NLs.CC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +   
  labs(x= "Strains", y = "NL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("NL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(15.9, 16.7, by=0.2), limits=c(15.9, 16.7)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p8 <- Variables %>% 
  dplyr::select(Strains, Age, PLs.CC, PLs.CC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = PLs.CC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = PLs.CC-PLs.CC.SE, ymax = PLs.CC+PLs.CC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "PL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("PL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(17.1, 17.5, by=0.1), limits=c(17.1, 17.5)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p9 <- Variables %>% 
  dplyr::select(Strains, Age, ELNLs.AlkylCC , ELNLs.AlkylCC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELNLs.AlkylCC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELNLs.AlkylCC-ELNLs.AlkylCC.SE, ymax = ELNLs.AlkylCC+ELNLs.AlkylCC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(strip.text = element_blank()) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELNL alk(en)yl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELNL alk(en)yl chain CC ")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(face = "bold", size = 7))+
  scale_y_continuous(breaks = seq(17.1, 20.5, by=0.8), limits=c(17.0,20.5)) +
  theme(axis.text.y = element_text(size = 8))+
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p10 <- Variables %>% 
  dplyr::select(Strains, Age, ELNLs.AcylCC , ELNLs.AcylCC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELNLs.AcylCC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELNLs.AcylCC-ELNLs.AcylCC.SE, ymax = ELNLs.AcylCC+ELNLs.AcylCC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(strip.text = element_blank()) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELNL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELNL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(16, 18, by=0.5), limits=c(16,18)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p11 <- Variables %>% 
  dplyr::select(Strains, Age, ELPLs.AlkylCC , ELPLs.AlkylCC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELPLs.AlkylCC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELPLs.AlkylCC-ELPLs.AlkylCC.SE, ymax = ELPLs.AlkylCC+ELPLs.AlkylCC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(strip.text = element_blank()) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELPL alk(en)yl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELPL alk(en)yl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(13.9, 14.3, by=0.1), limits=c(13.9,14.3)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p12 <- Variables %>% 
  dplyr::select(Strains, Age, ELPLs.AcylCC , ELPLs.AcylCC.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELPLs.AcylCC)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELPLs.AcylCC-ELPLs.AcylCC.SE, ymax = ELPLs.AcylCC+ELPLs.AcylCC.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(strip.text = element_blank()) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELPL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELPL acyl chain CC")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(21.8, 22.2, by=0.1), limits=c(21.8,22.2)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


#______________________________________Plots for double bond contents___________________________________

p13 <- Variables %>% 
  dplyr::select(Strains, Age, NLs.DB, NLs.DB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = NLs.DB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = NLs.DB-NLs.DB.SE, ymax = NLs.DB+NLs.DB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "NL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("NL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(0.4, 1.2, by=0.2), limits=c(0.4, 1.2)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p14 <- Variables %>% 
  dplyr::select(Strains, Age, PLs.DB, PLs.DB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = PLs.DB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = PLs.DB-PLs.DB.SE, ymax = PLs.DB+PLs.DB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +      
  labs(x= "Strains", y = "PL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("PL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(1.1, 1.5, by=0.1), limits=c(1.1,1.5)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p15 <- Variables %>% 
  dplyr::select(Strains, Age, ELNLs.AlkylDB , ELNLs.AlkylDB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELNLs.AlkylDB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELNLs.AlkylDB-ELNLs.AlkylDB.SE, ymax = ELNLs.AlkylDB+ELNLs.AlkylDB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELNL alkyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELNL alk(en)yl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(0, 1.3, by=0.3), limits=c(0,1.3)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p16 <- Variables %>% 
  dplyr::select(Strains, Age, ELNLs.AcylDB , ELNLs.AcylDB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELNLs.AcylDB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELNLs.AcylDB-ELNLs.AcylDB.SE, ymax = ELNLs.AcylDB+ELNLs.AcylDB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELNL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELNL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(0.4, 1.2, by=0.2), limits=c(0.4,1.2)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p17 <- Variables %>% 
  dplyr::select(Strains, Age, ELPLs.AlkylDB , ELPLs.AlkylDB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELPLs.AlkylDB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELPLs.AlkylDB-ELPLs.AlkylDB.SE, ymax = ELPLs.AlkylDB+ELPLs.AlkylDB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +
  labs(x= "Strains", y = "ELPL alkyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELPL alk(en)yl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(0.00, 0.08, by=0.02), limits=c(0.00, 0.08)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


p18 <- Variables %>% 
  dplyr::select(Strains, Age, ELPLs.AcylDB , ELPLs.AcylDB.SE) %>% 
  dplyr::mutate(Strains = factor(as.factor(Strains), levels = c("CTnew","CTold","SDnew","SDold","SDolder","CNnew","CNold"))) %>%
  ggplot(aes(x = Strains, y = ELPLs.AcylDB)) +
  geom_point(mapping = aes(colour = Strains, shape  = Strains), size = 1, position = pd) +
  geom_errorbar(aes(ymin = ELPLs.AcylDB-ELPLs.AcylDB.SE, ymax = ELPLs.AcylDB+ELPLs.AcylDB.SE, colour = Strains), width = 0.001, size = 0.4, alpha = 0.5) +
  theme_bw() +
  facet_wrap(~Age) +
  theme(axis.ticks.length=unit(0.015, "cm")) +
  theme(axis.ticks=element_line(size= 0.05)) +
  theme(strip.text = element_blank()) +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 7, hjust = 0.5)) +     
  labs(x= "Strains", y = "ELPL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  ggtitle("ELPL acyl chain DB")+
  theme(legend.position = "none")+ theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  theme(axis.title.x =element_blank())+ theme(axis.text.x =element_blank()) +
  #theme(axis.title.y = element_text(face = "bold", size = 8))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  #panel.grid.minor = element_blank(),
  #panel.border = element_rect(colour = "black", fill = NA, size=0.2))+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 7))+
  scale_y_continuous(breaks = seq(1.0, 1.9, by=0.2), limits=c(1.0, 1.9)) +
  scale_colour_manual(values = c("CTnew" = "red", "CTold" = "red","CNnew" = "green", "CNold" = "green",
                                 "SDnew" = "blue","SDold" = "blue", "SDolder" = "blue")) +
  scale_shape_manual(values = c("CTnew" = 16, "CTold" = 17,"CNnew" = 16, "CNold" = 17,
                                "SDnew" = 16,"SDold" = 17, "SDolder" = 15)) 


#_______________________________Arrange all plots_____________________________________________

F1 <- ggarrange(p2,p1,p3, nrow = 1, ncol=3)
F2 <- ggarrange(p4,p5,p6, nrow = 1, ncol=3)
F3 <- ggarrange(p7,p8,p9, nrow = 1, ncol=3)
F4 <- ggarrange(p10,p11,p12, nrow = 1, ncol=3)
F5 <- ggarrange(p13,p14,p15, nrow = 1, ncol=3)
F6 <- ggarrange(p16,p17,p18, nrow = 1, ncol=3)

# so first extract a legend that is laid out horizontally to add to combine plots

legend <- get_legend(p1.legend)
Figure <- ggarrange(F1,F2,F3,F4,F5,F6, ncol =1, nrow = 6)
Figure1 <- grid.arrange(Figure, legend, nrow = 2, layout_matrix = rbind(c(1,1), c(3,3)),
                        widths = c(2.7, 2.7), heights = c(2.5, 0.2))

ggsave(plot = Figure1, width = 170, height = 170, units = "mm", dpi = 300,filename = "Figures/Figure1.jpg")              

#________________________________END____________________________________________________________

