library(Rcpp)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(bio3d)
library(reshape2)
library(seqinr)
library(cowplot)

setwd("C:/Users/camil/OneDrive - Université Laval/MSc - Session 1 - Aut 2021/Analyse structure ERG11")
rm(list = ls())
gc()

#goal: get distances between amino acids in the protein and the substrate or the drugs

#fonction that takes a distance matrix and turns
#it into a two column table
flattenCorrMatrix2 <- function(cormat) {
  out <- melt(as.matrix(cormat), varnames = c("row", "col"))
  out$row <- as.numeric(as.character(out$row))
  out$col <- as.numeric(as.character(out$col))
  return(out)
}
#function that measures the euclidian distance 
#among points in a 3D space
measurePairwiseDist <- function(PDBatoms) {
  mat <- PDBatoms %>% select(x,y,z)
  PDBatoms$posnum <- seq(1:nrow(PDBatoms))
  t <- as.matrix(dist(mat), method="euclidian")
  tflat <- flattenCorrMatrix2(t)
  tflat <- left_join(tflat,PDBatoms, by=c("row"="posnum"))
  tflat <- left_join(tflat,PDBatoms, by=c("col"="posnum"))
  return(tflat) 
}


#5eqb
#Crystal structure of lanosterol 14-alpha demethylase with intact transmembrane domain bound to itraconazole
#4LXJ
#Saccharomyces cerevisiae lanosterol 14-alpha demethylase with lanosterol bound
#4WMZ
#S. cerevisiae CYP51 complexed with fluconazole in the active site


#for each dataset, keep the molecules, substrate, heme and inhibitor
#calculate for each amino acid the shortest distance to these molecules

#### LXJ
lxj <- read.pdb("4lxj") #two bound molecules; in resid: HEM and LAN (substrate)
#eliminate HOH and OXY resisdues and measure pairwise distance
dis_lxj <- lxj$atom %>% filter(resid !="HOH", resid !="OXY") %>% measurePairwiseDist()
#calculate minimal distance between the molecules and amino acids 
dis_lxj %<>% filter(type.x=="HETATM", type.y=="ATOM")
#for each type of molecule, measure the minimum distance for each residue 
min_dis_lxi <- dis_lxj %>% group_by(resid.x,resno.y,resid.y) %>% top_n(-1, value)
min_dis_lxi %<>% select(resno.y, resid.y,resid.x, value)                        

###eqp 
eqb <- read.pdb("5eqb") #two bound molecules; HEM and 1YN (itraconazole)
dis_eqb <- eqb$atom %>% filter(resid !="HOH", resid !="HEM") %>% measurePairwiseDist()
dis_eqb %<>% filter(type.x=="HETATM", type.y=="ATOM")
min_dis_eqb <- dis_eqb %>% group_by(resid.x,resno.y,resid.y) %>% top_n(-1, value)
min_dis_eqb %<>% select(resno.y, resid.y,resid.x, value)

##wmz
wmz <- read.pdb("4wmz") #two bound molecules; HEM and TPF (fluconazole) 
dis_wmz <- wmz$atom %>% filter(resid !="HOH", resid !="HEM") %>% measurePairwiseDist()
dis_wmz %<>% filter(type.x=="HETATM", type.y=="ATOM")
min_dis_wmz <- dis_wmz %>% group_by(resid.x,resno.y,resid.y) %>% top_n(-1, value)
min_dis_wmz %<>% select(resno.y, resid.y,resid.x, value)

#combine the three datasets and export to file

min_dis <- rbind(min_dis_eqb, min_dis_lxi, min_dis_wmz)
names(min_dis) <- c("res_num", "res", "molecule", "min_dist")

#a turn three letter codes to amino acid, a is a function in seqinr
min_dis$res <- a(str_to_title(min_dis$res))

write.table(min_dis, file="summary_structure_dist_substrate.txt", quote=FALSE, col.names = FALSE)
write.table(min_dis, file="summary_structure_dist_substrate2.txt", sep=",", quote=FALSE, col.names = FALSE)

library(openxlsx)
write.xlsx(min_dis,file="summary_structure_dist_substrate3.xlsx")
########################
#see distance on heatmap
########################

#here, keep HEM and LAN
dis_lxj <- lxj$atom %>% filter(resid !="HOH", resid !="OXY", resid !="HEM") %>% measurePairwiseDist()
#calculate minimal distance between the molecules and other molecules
dis_lxj %<>% filter(resid.x=="LAN", type.y=="ATOM")
#for each type of molecule, measure the minimum distance for each residue 
min_dis_lxi <- dis_lxj %>% group_by(resno.y,resid.y,eleno.x) %>% top_n(-1, value)
min_dis_lxi %<>% select(resno.y, resid.y,eleno.x, value)          


###eqp 
eqb <- read.pdb("5eqb") #two bound molecules; HEM and 1YN (itraconazole)
dis_eqb <- eqb$atom %>% filter(resid !="HOH", resid !="HEM") %>% measurePairwiseDist()
dis_eqb %<>% filter(resid.x=="1YN", type.y=="ATOM")
min_dis_eqb <- dis_eqb %>% group_by(resno.y,resid.y,eleno.x) %>% top_n(-1, value)
min_dis_eqb %<>% select(resno.y, resid.y,eleno.x, value)   


wmz <- read.pdb("4wmz") #two bound molecules; HEM and TPF (fluconazole) 
dis_wmz <- wmz$atom %>% filter(resid !="HOH", resid !="HEM") %>% measurePairwiseDist()
dis_wmz %<>% filter(resid.x=="TPF", type.y=="ATOM")
min_dis_wmz <- dis_wmz %>% group_by(resno.y,resid.y,eleno.x) %>% top_n(-1, value)
min_dis_wmz %<>% select(resno.y, resid.y,eleno.x, value)

h1 <- length(unique(min_dis_lxi$eleno.x))
l1 <- length(unique(min_dis_lxi$resno.y))

h2 <- length(unique(min_dis_eqb$eleno.x))
l2 <- length(unique(min_dis_eqb$resno.y))

h3 <- length(unique(min_dis_wmz$eleno.x))
l3 <- length(unique(min_dis_wmz$resno.y))


distance_lanosterol <- ggplot(min_dis_lxi, aes(x=resno.y, y=eleno.x))+
  geom_tile(data=min_dis_lxi, aes(fill=log2(value)),color='black')+
  scale_fill_gradient(low="red", high="black")+
  labs(x="Position", y="Position")+
  theme_cowplot()+
  theme(text = element_text(size=10))+
  coord_equal(ratio=3)+
  scale_x_continuous(n.breaks = 30,expand=c(0,0))+
  scale_y_continuous(n.breaks = 30,expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(title="Distance to lanosterol (log2)"))+
  theme(legend.position='top', 
        axis.text.y = element_text(size=4))
distance_lanosterol


distance_fluconazole <- ggplot(min_dis_wmz, aes(x=resno.y, y=eleno.x))+
  geom_tile(data=min_dis_wmz, aes(fill=log2(value)),color='black')+
  scale_fill_gradient(low="yellow", high="black")+
  labs(x="Position", y="Position")+
  theme_cowplot()+
  theme(text = element_text(size=10))+
  coord_equal(ratio=3*h1*l2/(l1*h2))+
  scale_x_continuous(n.breaks = 30,expand=c(0,0))+
  scale_y_continuous(n.breaks = 30,expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(title="Distance to Fluconazole (log2)"))+
  theme(legend.position='top', 
        axis.text.y = element_text(size=4))
distance_fluconazole

distance_itraconazole <- ggplot(min_dis_eqb, aes(x=resno.y, y=eleno.x))+
  geom_tile(data=min_dis_eqb, aes(fill=log2(value)),color='black')+
  scale_fill_gradient(low="blue", high="black")+
  labs(x="Position", y="Position")+
  theme_cowplot()+
  theme(text = element_text(size=10))+
  coord_equal(ratio=h1*l3/(l1*h3))+
  scale_x_continuous(n.breaks = 30,expand=c(0,0))+
  scale_y_continuous(n.breaks = 30,expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(title="Distance to Itraconazole (log2)"))+
  theme(legend.position='top', 
        axis.text.y = element_text(size=4))
distance_itraconazole

plot_grid(distance_lanosterol, distance_itraconazole, distance_fluconazole, 
          labels = c('A', 'B', "C"), nrow=3, label_size = 12,align = "v")



min_dis_lxi$drug <- rep("LAN", nrow(min_dis_lxi))
min_dis_eqb$drug <- rep("IT", nrow(min_dis_eqb))
min_dis_wmz$drug <- rep("FZ", nrow(min_dis_wmz))

all <- rbind(min_dis_lxi,min_dis_eqb,min_dis_wmz)

all_g <- all %>% group_by(resno.y, drug) %>% top_n(-1, value)


ggplot(all_g, aes(x=resno.y, y=as.factor(drug)))+
  geom_tile(data=all_g, aes(fill=log2(value)),color='black')+
  scale_fill_gradient(low="Red", high="black")+
  labs(x="Position", y="Molecule")+
  theme_classic()+
  coord_equal(ratio=20)+
  scale_x_continuous(n.breaks = 30,expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15),
        axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1, size=15))+
  guides(fill=guide_legend(title="Distance to molecule(log2)"))+
  theme(legend.position='top', 
        axis.text.y = element_text(size=4))

