
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(patchwork)


setwd("~/../Downloads/nov_fams_data/combined")


######
# dnds distribution (fig. 1B)
######


data = read.table("1B.tab",sep = '\t',header = T)
ggplot(data)+
  geom_histogram(aes(x = dnds),bins = 20,alpha = 0.2,color = "#56B4E9",fill = "#56B4E9")+
  xlab("dN / dS")+
  ylab("FESNov protein families")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))


#######
# average identity plot (fig. 1C)
#########

data = read.table("1C.tab",sep = '\t',header = T)
head(data)
ggplot(data)+
  geom_histogram(aes(x = av_id),alpha = 0.2,bins =  20,color = "#56B4E9",fill = "#56B4E9")+
  xlab("Average aminoacid identity")+
  ylab("FESNov protein families")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))+
  xlim(0,100)


####
# size distribution novel fams vs eggnog and small peptides (fig. 1D)
#####

data = read.table("1D.tab",header = T)

# plot lengh distribution
ggplot(data)+
#  coord_cartesian(xlim=c(0,1000)) +
  geom_histogram(aes(color = dataset,fill = dataset,x = gene_len),alpha = 0.3,position = "identity")+
  xlim(c(0,1000))+
  labs(fill="",color = "")+
  ylab(label = "Protein families")+
  xlab("Gene length")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))+
  scale_color_manual(values=c("#E69F00", "#56B4E9",'Salmon'),labels = c("Eggnog","Novel", "Small"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9",'Salmon'),labels = c("Eggnog","Novel",  "Small"))


####
# Number of species of FESNov fams vs eggnog (fig. 1E)
#####

data = read.table("1E.tab",header = T)

ggplot(data)+
  #coord_cartesian(xlim=c(0,50)) +
  geom_histogram(aes(color = dataset,fill = dataset,x = n_species),alpha = 0.3,position = "identity",bins = 50)+
  xlim(c(0,50))+
  labs(fill="",color = "")+
  ylab(label = "")+
  xlab("Number species")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))+
  scale_color_manual(values=c("#E69F00","#56B4E9"),labels = c("Eggnog","Novel"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9"),labels = c("Eggnog","Novel"))


####
# habitat distribution (fig. 3B)
#####

data_habs = read.table("3B_histogram.tab",header = T,sep = '\t')
head(data_habs)
data_mot = read.table("3B_mot.tab",header = T,sep = '\t')
head(data_mot)

ggplot()+
  geom_histogram(data = data_habs,aes(x = num_biomes,fill = rare),alpha = 0.6,position = "identity",bins = 70)+
  xlim(c(1,70))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac","grey"),name= "Rareness",labels = c("1-10 samples", "11-50 samples", "> 50 samples"))+
  geom_line(data = data_mot,aes(x = num_biomes, y = prop*100000,color = origin),size = 1)+
  scale_y_continuous(sec.axis = sec_axis(~.*1/100000, name = "Proportion in mobile"))+
  xlab("Habitats")+
  ylab("Unknown protein families")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))+
  scale_color_manual(labels = c("Plasmid","Viral") ,values = c("#CC6666", "#9999CC"),name= "Mobile")


####
# lineage specificity & motility (fig. 3C), also for Extended Data Figure S6
####

data_lca = read.table("3C_histogram.tab",header = T,sep = '\t')
head(data_lca)
data_mot = read.table("3C_mot.tab",header = T,sep = '\t')
head(data_mot)

data_lca$lin_sp = factor(data_lca$lin_sp,levels = c('r',"d","p","c","o","f","g"))
ggplot() +
  geom_bar(data = data_lca,aes(x = lin_sp, y = n),stat = 'identity',fill = 'grey80')+
  geom_line(data = data_mot,aes(x = lin_sp, y = prop*1000000,color = origin,group = origin),size = 2)+
  scale_y_continuous(sec.axis = sec_axis(~.*0.0000001, name = "Proportion in mobile"))+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20)) +
  scale_color_manual(labels = c("Plasmid","Viral") ,values = c("#CC6666", "#9999CC"),name= "Mobile")+
  ylab("Unknown protein families")+
  xlab("Lineage specificity")+
  scale_x_discrete(name ="Lineage specificity", 
                   labels=c("Multi\ndomain","Multi\nphyla","Multi\nclass","Multi\norder","Multi\nfamily","Multi\ngenus","Multi\nspecies"))



#####
# dnds synapomorphic vs non synapomorphic (Fig. 4B)
#####

data = read.table("4B.tab",header = T)
head(data)

ggplot(data,aes(x = syn,y = dnds)) +
  geom_boxplot(aes(x = syn,y = dnds))+
  geom_signif(comparisons = list(c("TRUE","FALSE")), 
              map_signif_level=TRUE)+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=20))+
  xlab("")+
  ylab("dN/dS")+
  scale_x_discrete(labels=c("TRUE" = "Synapomorfic", "FALSE" = "Non-synapomorfic"))+
  coord_flip()


#####
# maximum conservation score synapomorphic vs non synapomorphic  (Fig. 4C)
#####

data = read.table("4C.tab",header = T)
head(data)

ggplot(data,aes(x = syn, y = max_cons))+
  geom_boxplot()+
  geom_signif(
    comparisons = list(c("TRUE","FALSE")),
    map_signif_level = TRUE)+
  ylab("Maximum genomic context conservation score")+
  xlab("")+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=20))+
  scale_x_discrete(labels=c("TRUE" = "Synapomorfic", "FALSE" = "Non-synapomorfic"))
  

####
# CRC DA gene families (Fig 5A)
####

# load abs
d_plot = read.table("5.tab",header = T)
head(d_plot)


# plot
p1 = ggplot(d_plot %>% 
              arrange(mean_diff),aes(y = reorder(cluster,-mean_diff), x =  diff_crc_ctr))+
  geom_point(aes(shape = Study, color = direction))+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    size = 2
  )+
  stat_summary(position=position_dodge(0.95),geom="errorbar")+
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=15),
        axis.title=element_text(size=15),
  )  + 
  geom_vline(xintercept = 0)+
  xlab("Mean abundance difference")+
  ylab("Novel gene family")+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette="Dark2")



# p-values
d_prev = d_plot %>% 
  dplyr::select(cluster,pval,direction,mean_diff) %>% 
  unique()

p2 = ggplot(d_prev)+
  geom_bar(aes(y = pval, x = reorder(cluster,-mean_diff),fill = direction),stat = "identity")+
  coord_flip()+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=15),
        axis.title=element_text(size=15))+
  ylab("q-value")+
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(0, 0.01),breaks = c(0,0.01))+
  theme(legend.position = "none")


# prevalence
d_prev = d_plot %>% 
  dplyr::select(cluster,n_zeros,direction,mean_diff) %>% 
  unique()

p3 = ggplot(d_prev)+
  geom_bar(aes(y = (575-n_zeros) / 575, x = reorder(cluster,-mean_diff),fill = direction),stat = "identity")+
  coord_flip()+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=15),
        axis.title=element_text(size=15))+
  ylab("Prevalence")+
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(0, 1),breaks = c(0,0.5,1))+
  theme(legend.position = "none")


p1 + p2 + p3 +  plot_layout(guides = 'collect')+ plot_layout(widths = c(5,1,1))

#####
# Extended Data Figure S5
######

d_plot = read.table("S5.tab",header = T)

ggplot(d_plot) +
  geom_point(aes(x = num_habs, y = prop_mob,color = "red"))+
  ylim(c(0,0.4))+
  theme_classic()+
  annotate("text", x = 20, y = 0.3, label = paste("R = ",as.character(round(cor(d_plot$num_habs,d_plot$prop_mob,method = "spearman"),2))),size = 10)+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20),
    legend.position = "none")+
  xlab("Number of habitats")+
  ylab("Proportion of unknown protein families in\n plasmid / viral contigs")


######
# Extended Figure S6
######

# S6a

data_lca = read.table("S6A_histogram.tab",header = T,sep = '\t')
head(data_lca)
data_mot = read.table("S6A_mot.tab",header = T,sep = '\t')
head(data_mot)

data_lca$lin_sp= factor(data_lca$lin_sp,levels = c('r',"d","p","c","o","f","g"))
ggplot() +
  geom_bar(data = data_lca,aes(x = lin_sp, y = n),stat = 'identity',fill = 'grey80')+
  geom_line(data = data_mot,aes(x = lin_sp, y = prop*1000000,color = origin,group = origin),size = 2)+
  scale_y_continuous(sec.axis = sec_axis(~.*0.0000001, name = "Proportion in mobile"))+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20)) +
  scale_color_manual(labels = c("Plasmid","Viral") ,values = c("#CC6666", "#9999CC"),name= "Mobile")+
  ylab("Unknown protein families")+
  xlab("Lineage specificity")+
  scale_x_discrete(name ="Lineage specificity", 
                   labels=c("Multi\ndomain","Multi\nphyla","Multi\nclass","Multi\norder","Multi\nfamily","Multi\ngenus","Multi\nspecies"))


# S6b
data_lca = read.table("S6B_histogram.tab",header = T,sep = '\t')
head(data_lca)
data_mot = read.table("S6B_mot.tab",header = T,sep = '\t')
head(data_mot)

data_lca$lin_sp= factor(data_lca$lin_sp,levels = c('r',"d","p","c","o","f","g"))
ggplot() +
  geom_bar(data = data_lca,aes(x = lin_sp, y = n),stat = 'identity',fill = 'grey80')+
  geom_line(data = data_mot,aes(x = lin_sp, y = prop*1000000,color = origin,group = origin),size = 2)+
  scale_y_continuous(sec.axis = sec_axis(~.*0.0000001, name = "Proportion in mobile"))+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20)) +
  scale_color_manual(labels = c("Plasmid","Viral") ,values = c("#CC6666", "#9999CC"),name= "Mobile")+
  ylab("Unknown protein families")+
  xlab("Lineage specificity")+
  scale_x_discrete(name ="Lineage specificity", 
                   labels=c("Multi\ndomain","Multi\nphyla","Multi\nclass","Multi\norder","Multi\nfamily","Multi\ngenus","Multi\nspecies"))



######
# Extended Data Figure S7
######

# S7a
tsne = read.table("S7A.tab",sep = '\t',header = T)
head(tsne)

cols_c = c("darkorange3","cornsilk4","blue3","plum3","lawngreen","gold1","azure4","mediumpurple","lightblue1","cadetblue1","orchid4","cornflowerblue","seashell2","yellowgreen","lightskyblue4","pink3","olivedrab","seashell4","ivory","paleturquoise1","lightskyblue3","cyan1","pink1","lightpink","mistyrose1","orangered4","palevioletred1","goldenrod3","thistle","dodgerblue","mintcream","darksalmon","honeydew1","tan1","maroon1","springgreen2","mediumorchid2","goldenrod1","white","darkolivegreen2","mediumslateblue","thistle1","mediumpurple4","peachpuff3","thistle4","turquoise","moccasin","orangered1","rosybrown")
p1 = ggplot(tsne)+
  geom_point(aes(x = tsne1,y = tsne2,col = country),size = 1)+
  theme_classic()+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  theme(
    text=element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=15))+
  scale_color_manual(values = cols_c)

# S7b
tsne = read.table("S7B.tab",sep = '\t',header = T)
head(tsne)

p2 = ggplot(tsne)+
  geom_point(aes(x = tsne1,y = tsne2,col = country),size = 1)+
  theme_classic()+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  theme(
    text=element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=15))+
  scale_color_manual(values = cols_c)


# FS7c
tsne = read.table("S7C.tab",sep = '\t',header = T)
head(tsne)
col_c = c("olivedrab4","hotpink3","gold1","mediumpurple3","deepskyblue","antiquewhite4",
          "blue", "hotpink1","deeppink3","orange4","chartreuse4","snow3", "goldenrod")

p3 = ggplot(tsne)+
  geom_point(aes(x = tsne1,y = tsne2,col = habitat),size = 2)+
  theme_classic()+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  theme(
    text=element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=15))+
  scale_color_manual(values = col_c)  


# S7D
tsne = read.table("S7D.tab",sep = '\t',header = T)

p4 = ggplot(tsne)+
  geom_point(aes(x = tsne1,y = tsne2,col = habitat),size = 2)+
  theme_classic()+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  theme(
    text=element_text(size=15),
    axis.title=element_text(size=15),
    axis.text=element_text(size=15))+
  scale_color_manual(values = col_c)

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")


######
# Extended Data Figure S9
######


# S9A

data = read.table("S9A.tab",header = T,sep = "\t")


data$dataset = factor(data$dataset,levels = c("KOs", "FESNov\n fams","KOs + \nFESNov fams"))
p1 = ggplot(data,aes(x = dataset,y = auc))+
  geom_boxplot()+
  theme_classic()+
  geom_signif(
    comparisons = list(c("KOs", "KOs + \nFESNov fams")),
    map_signif_level = TRUE
  ) +
  theme(text=element_text(size=15),
        axis.title=element_text(size=15),
  )+
  ylab ("AUC")+
  xlab("Dataset")


# S9B

d = read.table("S9B.tab",header = T,sep = "\t")
head(data)

d$dataset = factor(d$dataset,levels = c("KOs", "FESNov\n  fams","KOs + FESNov\n  fams"))
p2 = ggplot(d,aes(x = dataset,y = auc))+
  geom_boxplot()+
  theme_classic()+
  theme(text=element_text(size=15),
        axis.title=element_text(size=15),
  )+
  ylab ("AUC")+
  xlab("Dataset")+
  scale_fill_brewer(palette="Spectral")+
  geom_signif(
    comparisons = list(c("KOs", "FESNov\n  fams")),
    map_signif_level = TRUE
  ) 

?wilcox.test
p1 + p2
