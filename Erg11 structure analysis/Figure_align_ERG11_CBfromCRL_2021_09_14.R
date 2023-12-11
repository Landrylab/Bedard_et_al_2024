library(tidyverse)
library(magrittr)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(data.table)
library(seqinr)
library(ape)
library("Biostrings")
library(httr)


.libPaths()
library(clock)
tzdb::tzdb_names()


setwd("C:/Users/camil/OneDrive - Université Laval/MSc - Session 1 - Aut 2021/Analyse structure ERG11")

##### analyze alignment ####
#import
df <- read.alignment("ERG11_msa.fasta", format="fasta")

#calculate distance 
d <- as.matrix(dist.alignment(df, "identity"))
# the squared root of the pairwise distances. 
#For example, if identity between 2 sequences is 80 the squared root of (1.0 - 0.8) i.e. 0.4472136.

id_from_scer <- as.data.frame(d[1:ncol(d),"saccharomyces_cerevisiae", drop=FALSE])
id_from_scer <- tibble::rownames_to_column(id_from_scer, "VALUE")
names(id_from_scer) <- c("species", "distance")
id_from_scer %<>% mutate(distance = distance^2)

pl_id_cer <- id_from_scer %>% ggplot(aes(x=1-distance))+
  geom_histogram(fill="#00BFFF")+
  labs(
    x = "AA similarity",
    y = "Number of orthologs"
  )+
  theme_cowplot()+
  theme(legend.position = c(0.8, 0.8),
        text = element_text(size=20))
pl_id_cer

#work with sequences
df <- data.frame(t(as.matrix(df)))
df %<>% dplyr::rename(erg11=saccharomyces_cerevisiae)

#list of species used in the database MARDY
mardy <- c("Aspergillus_fumigatus",
"Pyrenopeziza_brassicae",
"Candida_albicans",
"Candida_auris",
"Cryptococcus_neoformans",
"Candida_tropicalis",
"Molilinia_fructicola",
"erg11",
"Candida_parapsilosis",
"Candida_glabrata",
"Candida_krusei",
"Aspergillus_flavus",
"Mycosphaerella_graminicola",
"Mycosphaerella_fijiensis")

mardy <- tolower(mardy)

mardy_species <- intersect(mardy, names(df))

#change order
df %<>% select(erg11, 
               candida_glabrata,
               candida_tropicalis,
               candida_auris,
               candida_albicans,
               candida_parapsilosis,
               cryptococcus_neoformans,
               aspergillus_fumigatus,
               everything())

df %<>% mutate_all(as.character)

#define positions in the alignment and put as first column
df$pos <- seq(1:nrow(df))
df %<>% select(pos, everything())


#create a function to have pos number of each seq rel to alignment
#this will allow to know which residue of each ortholog
#corresponds to each in the alignment

pos_seq <- function(df_used, name_species_input){
  temp_df <- df_used %>% select(pos, name_species_input)
  names(temp_df) <- c("pos", "seq")
  temp_df %<>% filter(seq != "-")  
  temp_df$pos_seq <- seq(1:nrow(temp_df))
  temp_df$pos_seq <- paste(temp_df$pos_seq, temp_df$seq, sep="_")
  names(temp_df) <- c("pos", name_species_input, paste(name_species_input, "_pos", sep=""))
  return(temp_df)
}

#process all sequences
for (spe in names(df)[2:length(names(df))]) {
  temp <- pos_seq(df, spe)
  df <- merge(df, temp, by = intersect(names(df), names(temp)), all=TRUE)
}

#use columns with positions and then split to have position and aa
dfrec <- df %>% select(pos, erg11, ends_with("_pos"))

#bring to long format
dfl <- dfrec %>% pivot_longer(cols=3:ncol(dfrec), names_to="species", values_to="aa") %>%
              mutate(aa=toupper(aa))
names(dfl)[1] <- "aln"

dfl %<>% separate(aa, c("pos", "res"), sep="_") 

dfl %<>% mutate(species = str_replace(species, "_pos",""))

#error in naming columns
dfl$erg11 = str_to_upper(dfl$erg11)

#will need a column with the position of the Erg11 residue pos when I need to add information
#on mutation effects for instance
new_cord <- dfl %>% filter(erg11 !="-", species=="erg11") %>% select(aln, pos)
names(new_cord) <- c("aln", "pos_erg11_aln")

dfl <- left_join(dfl, new_cord, by=c("aln"="aln"))

dfl %<>% select(pos_erg11_aln, everything())

##### Read Mardy DATA####

mardy <- read_csv("DB_by_gene.csv", skip=1)

genes <- unique(mardy$`Gene name`)

erg11ortho <- c(genes[1], genes[2], genes[3], genes[7], genes[8], genes[9], genes[13])

#gene names to look at

mardy %<>% filter(`Gene name` %in% erg11ortho) %>%
                filter(is.na(`Tandem repeat name`)) %>%
                mutate(species = str_replace(Organism, " ","_")) %>%
                mutate(species = tolower(species)) 

#extract positions of mutations
#skip lines where more than one mutation
#split 

code_mut <- mardy %>% filter(str_length(`AA mutation`)<6) %>%
                mutate(wt_res = str_sub(`AA mutation`, 1, 1),
                       mut_res = str_sub(`AA mutation`, -1, -1),
                       pos_mut = str_remove_all(`AA mutation`, "[A-Za-z]+")) %>%
                       select(species, Drug, wt_res, mut_res, pos_mut)

dfl <- left_join(dfl, code_mut, by=c("species"="species", "pos"="pos_mut"))
dfl <- as_tibble(dfl)

library(rJava)
library(xlsx)
# Write the first data set in a new workbook
write.xlsx(dfl, file = "mardyalign.xlsx", append = FALSE)

####import structure data, distance from molecules####
PDB <- read_delim("summary_structure_dist_substrate.txt", delim=" ", col_names=FALSE)
names(PDB) <- c("id", "pos_erg11", "res", "molecule", "min_dist")

mdid <- PDB %>% group_by(pos_erg11) %>% summarise(mean_dist = mean(min_dist)) 

  

#test if coordinates match
test <- dfl %>% filter(erg11 != "-", species=="erg11")

test <- left_join(test, PDB, by=c("pos"="pos_erg11"))

#need PDB to be wide before continuing
PDBw <- PDB %>% select(pos_erg11, res, molecule, min_dist) %>% pivot_wider(names_from=molecule, values_from = min_dist)

names(PDBw) <- c("pos_erg11", "res", "d_itraco", "d_heme", "d_lano", "d_fluco")
PDBw %<>% select(-res)
PDBw$pos_erg11 <- as.character(PDBw$pos_erg11)

dfl <- left_join(dfl, PDBw, by=c("pos_erg11_aln"="pos_erg11"))

#plot distance from molecules for resistance mutations
pos_erg11_w_mutres <- dfl %>% filter(!is.na(mut_res)) %>% select(aln) 

for_plot_dist <- dfl %>% 
      select(aln, d_itraco, d_heme, d_lano, d_fluco) %>% 
      distinct(aln, d_itraco, d_heme, d_lano, d_fluco) %>%
      pivot_longer(cols=2:5, names_to="dist_type", values_to="dist") %>%
      mutate(resis_mut = ifelse(aln %in% pos_erg11_w_mutres$aln, "resistance","no_resistance")) %>%
      filter(!is.na(dist))

test <- (table(for_plot_dist$aln, for_plot_dist$resis_mut))
intersect(for_plot_dist$aln[for_plot_dist$resis_mut=="no_resistance"],for_plot_dist$aln[for_plot_dist$resis_mut=="resistance"])

plt_dis_res <- ggplot(for_plot_dist, aes(x=resis_mut, y=dist, fill=dist_type)) +
              geom_boxplot(outlier.shape = NA)+
              scale_fill_manual(values=c("#0000CD", "#8B4513","#00BFFF","#F4A460"))+
              geom_jitter(alpha=0.1)+
              labs(x="Site in ScErg11", y="Distance from molecule",
                  fill="Molecule")+
              theme_cowplot()+
              theme(legend.position = c(0.8, 0.8),
                    text = element_text(size=20))
plt_dis_res

####GEMME virtual ####

G <- read_delim("saccharomyces_normPred_evolCombi.txt", delim=" ", skip=1, col_names = FALSE)
names(G) <- seq(1:(ncol(G)))-1
Gl <- G %>% pivot_longer(cols=2:ncol(G), names_to="pos", values_to="Geffect")
names(Gl) <- c("aa", "pos", "Geffect") 


ggplot(Gl, aes(x=as.numeric(as.character(pos)), y=aa, fill=Geffect))+
  geom_tile(color='black')+
  scale_fill_gradient(low='red', high='blue')+
  labs(x="Position", y="Amino Acid")+
  theme_bw()

#mark positions with resistance mutations
pos_mutres <- dfl %>% filter(!is.na(mut_res)) %>% 
              select(pos_erg11_aln) %>%
              distinct(pos_erg11_aln)

library(readxl)

posDMS<- read_excel("C:/Users/CABED117/OneDrive - Université Laval/MSc - Session 1 - Aut 2021/Analyse structure ERG11/positions_for_ScERG11_DMS_CB_2021-09-20.xlsx", 
                                                      col_types = c("numeric", "text", "numeric"))

pos_mutres$y <- rep(0.5, nrow(pos_mutres))              


Gl <- left_join(Gl, pos_mutres, by=c("pos"= "pos_erg11_aln"))              

plt_heatm <- ggplot(Gl, aes(x=as.numeric(as.character(pos)), y=aa))+
  geom_raster(data=Gl, aes(fill = Geffect), color='black')+
  scale_fill_gradient(low="yellow", high="blue")+
  labs(x="Position in ScErg11p", y="Amino Acid")+
  theme_cowplot()+  scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400,450,500))+
  geom_point(data = posDMS, 
             aes(x = pos, 
                 y = 20.9),
             shape=6, color="red")+
    theme(text = element_text(size=20))+
   coord_equal(ratio=10)

plt_heatm

#plot with distance from drugs
to_plot <- dfl %>% filter(!is.na(pos_erg11_aln), species=="erg11") %>% 
  select(pos_erg11_aln, d_itraco, d_heme, d_lano, d_fluco) %>%
  distinct(pos_erg11_aln, d_itraco, d_heme, d_lano, d_fluco) %>%
  pivot_longer(cols=2:5, names_to="dist_type", values_to="dist") 

pos_mutres$y <- rep(2, nrow(pos_mutres)) 

to_plot <- left_join(to_plot, pos_mutres, by=c("pos_erg11_aln"= "pos_erg11_aln"))  


library(readxl)

pos_unique_heme <- read_excel("C:pos_unique_heme.xlsx")



line_plot <- ggplot(to_plot, aes(x=as.numeric(as.character(pos_erg11_aln)), y=dist, col=as.factor(dist_type)))+
       geom_line(lwd=0.5)+ geom_hline(yintercept = 12) +
       scale_color_manual(values=c("#6C44A6","#CF6448","#417471", "#E9D40D"), breaks = c("d_fluco","d_itraco","d_lano","d_heme"), labels = c("Fluconazole","Itraconazole","Lanosterol", "Heme") )+
       theme_cowplot()+  theme(text=element_text(family = "Times New Roman"))+
       labs(x="Position in ScErg11", y="Distance from molecule (Å)", col="Molecules")+ scale_y_continuous(breaks = c(0,5,10,15,20,40,60)) +
       scale_x_continuous(breaks = c(0,50,100,150,200,250,300,350,400,450,500,530)) +
       geom_point(data = to_plot, 
           aes(x = as.numeric(as.character(pos_erg11_aln)), 
               y = y),
           shape=4, color="black") 


line_plot 

#add data on sensitivity to mutations
ave_effect <- Gl %>% select(pos, Geffect) %>%
              group_by(pos) %>% summarise(ave_g=mean(Geffect, na.rm=T))

##combine data with Dfl to look at the expected effects of resistance mutations 
ggplot(ave_effect, aes(x=as.numeric(as.character(pos)), y=ave_g))+
  geom_line()+
  theme_cowplot()+
  labs(x="Position", y="Sensitivity to mutations")

ave_effect %<>% mutate(res = ifelse(pos %in% pos_mutres$pos_erg11_aln, 1,0))

eff_mut_side <- ggplot(ave_effect, aes(x=as.factor(res), y=ave_g))+
  geom_boxplot(outlier.shape = NA, fill="#00BFFF")+ scale_x_discrete(labels = c("no", "yes"))+
  geom_jitter(alpha=0.1)+
  theme_cowplot()+
  labs(x="Resistance site", y="Robusness to mutations")+
  theme(text = element_text(size=20))

eff_mut_side

#combine dfl with Gl to look at whether the effect of the substitutions associated 
#with resistance are expected to be strong or not. 
Gl %<>% select(-y) %>% mutate(aa = toupper(aa))

#positions in Gl correspond to positions of the Scerg11  
#select positions where the WT residue of the orthologs are the same as erg11
#remove erg11 from line

species_w_res_mu <- dfl %>% filter(!is.na(mut_res)) %>% select(species) %>% distinct(species)

res_evo <- dfl %>% filter(!is.na(pos_erg11_aln), species !="erg11" ) %>%
                   filter(erg11 == res, species %in% species_w_res_mu$species)

res_evo <- left_join(res_evo, Gl, by=c("pos_erg11_aln"="pos", "mut_res"="aa"))
                   
##check super conserved sites

sumevo <-dfl %>% filter(!is.na(pos_erg11_aln)) %>% group_by(pos_erg11_aln, res) %>% tally() %>% mutate(res = ifelse(is.na(res), "gap", res))

#sumevo %>% pivot_wider(names_from=res, values_from = n) %>% rename(gap=`NA`)

ggplot(sumevo , aes(x=as.numeric(as.character(pos_erg11_aln)), y=res))+
  geom_raster(data=sumevo, aes(fill = log2(n)), color='black')+
  scale_fill_gradient(low="#FFD700", high="#00BFFF")+
  labs(x="Position", y="Amino Acid")+
  theme_cowplot()+
  theme(text = element_text(size=20))+
  coord_equal(ratio=10)+
  scale_x_discrete(expand=c(0,0))

ggplot(sumevo , aes(x=as.numeric(as.character(pos_erg11_aln)), y=res))+
  geom_col(data=sumevo, aes(fill = n))+
  labs(x="Position", y="Amino Acid")+
  theme_cowplot()+
  theme(text = element_text(size=20))+
  coord_equal(ratio=10)+
  scale_x_discrete(expand=c(0,0))


