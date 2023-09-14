
library(ggplot2)
library(dplyr)
library(tidyr)

#####
# load data
#####

# filtered fams
f_fams = read.table("ffams.txt")


# dnds for each fam
dnds_data = read.table("busted_final.tab")
names(dnds_data) = c("fan","dnds")
dnds_data = dnds_data %>% 
  filter(fan %in% f_fams$V1)
dnds_data$dnds = as.numeric(dnds_data$dnds)

# number members per fam
size_dist_data = read.table("microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.rep_clu_seqs.size_dist.tab")
names(size_dist_data) = c("cluster","value")
size_dist_data$variable = "Novel"
size_dist_data = size_dist_data %>% 
  filter(cluster %in% f_fams$V1)

# plasflow predictions per fam
plasflow = read.table("microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.plasflow_0.95.tab")
head(plasflow)
plasflow = plasflow %>% 
  mutate(plasmid = case_when(V4 > 0 ~ TRUE,
                            V4 == 0  ~ FALSE)) %>%
  select(V1, plasmid)

# seeker predictions per fam 
seeker = read.table("microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.seeker_0.95.tab")
seeker = seeker %>% 
  mutate(viral = case_when(V4 >0 ~ TRUE,
                             V4==0 ~ FALSE)) %>% 
  select(V1, viral)

# lcas per fam 
lcas = read.table("lca_per_fam.rank.tab")
names(lcas) = c("V1","lca")

# number habitats per fam  
num_biomes = read.table("number_biomes_per_fam.eval_1e-3_cov50.tab") #
names(num_biomes) = c("V1","num_biomes")
head(num_biomes)
nrow(num_biomes)

# number of detections in per fam
number_samples = read.table("number_samples_per_fam.eval_1e-3_cov_50.tab")
names(number_samples) = c("V1","num_samples")

# number habitats per fam (strict filters)
number_biomes_90 = read.table("number_biomes_per_fam.id_90_cov_50.tab")
names(number_biomes_90) = c("V1","num_biomes_90")

# number of detections per fam (strong filters)
number_samples_90 = read.table("number_samples_per_fam.id_90_cov_50.tab")
names(number_samples_90) = c("V1","num_samples_90")

# alistat stats for joining with other data
data_alistat  = read.table("../alg_stats/alistat_stats.tab",header = T)
n = names(data_alistat)
n[1] = "V1"
names(data_alistat) = n
data_alistat = data_alistat %>% 
  filter(V1 %in% f_fams$V1)


# join  data
data = filt_fams %>% 
  left_join(plasflow,by = "V1") %>% 
  left_join(seeker,by = "V1") %>% 
  left_join(lcas,by = "V1") %>% 
  left_join(num_biomes,by = "V1") %>% 
  mutate(multibiome = case_when(num_biomes>1 ~ TRUE,
                           num_biomes<=1 ~ FALSE)) %>% 
  mutate(mobile_vir = case_when(plasmid == T | viral == T ~ T,
                                plasmid == F & viral == F ~ F)) %>% 
  mutate(synapo = case_when(V1 %in% synapo_fams$V1 ~ T,
                            ! V1 %in% synapo_fams$V1 ~ F)) %>% 
  left_join(synapo_fams,by = "V1") %>% 
  left_join(data_alistat, by = "V1") %>% 
  left_join(number_samples,by = "V1") %>% 
  mutate(rare = case_when(num_samples < 10 ~'< 10 samples',
                   num_samples>= 10 & num_samples <=50 ~ '11-50 samples',
                   num_samples >50 ~ ">50 samples")) %>% 
  left_join(number_biomes_90,by = "V1") %>% 
  left_join(number_samples_90,by = "V1") %>% 
  mutate(rare_90 = case_when(num_samples_90 < 10 ~'< 10 samples',
                          num_samples_90 >= 10 & num_samples_90 <=50 ~ '11-50 samples',
                          num_samples_90 >50 ~ ">50 samples"))





######
# dnds distribution (fig. 1B)
######

ggplot(dnds_data)+
  geom_histogram(aes(x = dnds),bins = 20,alpha = 0.2,color = "#56B4E9",fill = "#56B4E9")+
  xlab("dN / dS")+
  ylab("Unknown protein families")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))

#######
# average identity plot (fig. 1C)
#########

ggplot(stats)+
  geom_histogram(aes(x = av_id),alpha = 0.2,bins =  20,color = "#56B4E9",fill = "#56B4E9")+
  xlab("Average identity")+
  ylab("") + 
  theme_classic()+
  theme(
    text=element_text(size=40),
    axis.title=element_text(size=30),
    axis.text=element_text(size=20))+
  xlim(0,100)

####
# size distribution novel fams vs eggnog and small peptides (fig. 1D-E)
#####

# load short pep and eggnog lens
short_pep = read.table("sequence_len_small_peptide.txt",sep = '\t',header = T)
short_pep$cluster = "small"
short_pep$variable = "Small"
short_pep$value = short_pep$Length.of.peptide
short_pep$Length.of.peptide = NULL

eggnog = read.table("per_fam_len.all.tab")
names(eggnog) = c("cluster","value")
eggnog$variable = "Eggnog"

# combine
data_all = rbind(size_dist_data,short_pep,eggnog) %>% 
  filter(!variable %in%  "conensous_len")

# plot
ggthemr("dust")
ggplot(data_all[order(data_all$variable, decreasing = F),])+
  geom_histogram(aes(color = variable,fill = variable,x = value),alpha = 0.3,position = "identity")+
  xlim(c(0,1000))+
  labs(fill="",color = "")+
  ylab(label = "Unknown protein families")+
  xlab("Gene length")+
  theme_classic()+
  theme(
    text=element_text(size=20),
    axis.title=element_text(size=20),
    axis.text=element_text(size=20))+
  scale_color_manual(values=c("#E69F00", "#56B4E9",'Salmon'),labels = c("Eggnog","Novel", "Small"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9",'Salmon'),labels = c("Eggnog","Novel",  "Small"))


# sub plot 
data_num_sp = read.table("microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.num_sp.tsv")
names(data_num_sp) = c("cluster","value")
data_num_sp$variable = "Novel"

data_num_sp = data_num_sp %>% 
  filter(cluster %in% f_fams$V1)

# eggnog
eggnog = read.table("eggnog_num_tax.all.tab")
names(eggnog) = c("cluster","value")
eggnog$variable = "Eggnog"

# combine & plot
data_all = rbind(data_num_sp,eggnog) %>% 
  filter(!variable %in%  "conensous_len")
ggplot(data_all[order(data_all$variable, decreasing = F),])+
  geom_histogram(aes(color = variable,fill = variable,x = value),alpha = 0.3,position = "identity",bins = 50)+
  xlim(c(0,50))+
  labs(fill="",color = "")+
  ylab(label = "")+
  xlab("Number species")+
  theme_classic()+
  theme(
    text=element_text(size=40),
    axis.title=element_text(size=30),
    axis.text=element_text(size=20))+
  scale_color_manual(values=c("#E69F00","#56B4E9"),labels = c("Eggnog","Novel"))+
  scale_fill_manual(values=c("#E69F00","#56B4E9"),labels = c("Eggnog","Novel"))

####
# habitat distribution (fig. 3B)
#####

data_multi_p = data %>% 
  group_by(plasmid,num_biomes) %>% 
  summarise(n = n()) %>% 
  group_by(num_biomes) %>% 
  mutate(prop_p = proportions(n)) %>% # calculate proportions num_biomes  
  filter(! is.na(num_biomes)) %>% 
  filter(plasmid == TRUE) %>% 
  filter(num_biomes<80)


data_multi_v = data %>% 
  group_by(viral,num_biomes) %>% 
  summarise(n = n()) %>% 
  group_by(num_biomes) %>% 
  mutate(prop_v = proportions(n)) %>% # calculate proportions num_biomes  
  filter(! is.na(num_biomes)) %>% 
  filter(viral == TRUE) %>% 
  filter(num_biomes<80)

data_multi_mov =  merge(data_multi_p,data_multi_v,by = "num_biomes") %>% 
  gather(key = "origin",value = "prop",prop_v,prop_p)


ggplot()+
  geom_histogram(data = data,aes(x = num_biomes,fill = rare),alpha = 0.6,position = "identity",bins = 70)+
  xlim(c(1,70))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac","grey"),name= "Rareness",labels = c("1-10 samples", "11-50 samples", "> 50 samples"))+
  geom_line(data = data_multi_mov,aes(x = num_biomes, y = prop*100000,color = origin),size = 1)+
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
# lineage specificity (fig. 3C)
####


# number of fams spaning multiple phylums / in the same genus
(nrow(data[data$lca == 's',]) + nrow(data[data$lca == 'g',]))  / nrow(data)
(nrow(data[data$lca == 'd',]) + nrow(data[data$lca == 'r',]))  / nrow(data)
(nrow(data[data$lca == 'r',]))  / nrow(data)

# mobile stats for root and genus
nrow(data[data$plasmid==T,]) / nrow(data)
nrow(data[data$lca=='r' & data$plasmid==T,]) / nrow(data[data$lca=='r',])
nrow(data[(data$lca=='r' | data$lca == 'd') & data$plasmid==T,]) / nrow(data[data$lca=='r' | data$lca == 'd' ,])
nrow(data[data$lca=='g' & data$plasmid==T,]) / nrow(data[data$lca=='g',])

# viral stats for root and genus 
nrow(data[data$viral==T,]) / nrow(data)
nrow(data[data$lca=='r' & data$viral==T,]) / nrow(data[data$lca=='r',])
nrow(data[(data$lca=='r' | data$lca == 'd') & data$viral==T,]) / nrow(data[data$lca=='r' | data$lca == 'd' ,])
nrow(data[data$lca=='g' & data$viral==T,]) / nrow(data[data$lca=='g',])

# identity in different lcas
mean(data[data$lca %in% c('s','g'),]$av_id)
mean(data[data$lca %in% c('d','r'),]$av_id)
(nrow(data[data$lca == 'r',]))  / nrow(data)

# stats for multibiome
data_st = data %>% 
  filter(!is.na(data$num_biomes))
nrow(data_st[data_st$num_biomes == 40 & data_st$viral ==T,]) / nrow(data_st[data_st$num_biomes==40,])


# plot
data_multi_p = data %>% 
  group_by(lca,plasmid) %>% 
  summarise(n = n()) %>% 
  group_by(lca) %>% 
  mutate(prop_p = proportions(n)) %>% 
  filter(!lca %in% 's') %>% 
  group_by(lca) %>% 
  mutate(n1 = sum(n)) %>% 
  filter(plasmid == TRUE)

data_multi_v = data %>% 
  group_by(lca,viral) %>% 
  summarise(n = n()) %>% 
  group_by(lca) %>% 
  mutate(prop_v = proportions(n)) %>% 
  filter(!lca %in% 's') %>% 
  group_by(lca) %>% 
  mutate(n1 = sum(n)) %>% 
  filter(viral == TRUE)

data_multi_mov =  merge(data_multi_p,data_multi_v,by = "lca") %>% 
  gather(key = "origin",value = "prop",prop_v,prop_p)

data_multi_p$lca = factor(data_multi_p$lca,levels = c('r',"d","p","c","o","f","g"))
ggplot() +
  geom_bar(data = data_multi_p,aes(x = lca, y = n1),stat = 'identity',fill = 'grey80')+
  geom_line(data = data_multi_mov,aes(x = lca, y = prop*1000000,color = origin,group = origin),size = 2)+
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
# dnds synapomorphic vs non synapomorphic (Fig. 3C)
#####

# mobile stats for synapos / non synapos
nrow(data[data$plasmid==T & data$synapo == T & data$synapo_lev %in% c('p','c','o'),]) / nrow(data[data$synapo == T & data$synapo_lev %in% c('p','c','o'),])
nrow(data[data$plasmid==T,]) / nrow(data)

syn_fams = read.table("synapo_fams.tab",sep = '\t')
names(syn_fams) = c("lev","lev_name","fam")
syn_fams = syn_fams %>% 
  mutate(syn = case_when(lev %in% c('p','c','o')~TRUE,
                         !lev %in% c('p','c','o')~FALSE)) %>% 
  filter(syn == T)

data_c = dnds_data %>% 
  mutate(syn = case_when(fam %in% syn_fams$fam ~ TRUE,
                         !fam %in% syn_fams$fam ~ FALSE)) %>% 
    filter( fam %in% filt_fams$V1 )

ggplot(all_d,aes(x = syn,y = dnds)) +
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
  scale_x_discrete(labels=c("TRUE" = "Synapomorfic", "FALSE" = "Non-synapomorfic"))

#####
# maximum conservation score synapomorphic vs non synapomorphic  (Fig. 3D)
#####

data = read.table("max_score_prev_pos_per_fam.tab")
data = data %>% 
  filter(fam %in% f_fams$V1) %>% 
  mutate(syn = case_when(V1 %in% syn$V1 ~"synapomorphic",
                         !V1 %in% syn$V1 ~ "non-synapomorphic"))

ggplot(data,aes(x = syn, y = max_cons))+
  geom_boxplot(aes(x = syn, y = V2))+
  geom_signif(
    comparisons = list(c("synapomorphic","non-synapomorphic")),
    map_signif_level = TRUE)+
  ylab("Maximum genomic context conservation score")+
  xlab("")+
  coord_flip()+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(size=20),
        axis.title=element_text(size=20),
        axis.text=element_text(size=20))


####
# CRC DA gene families (Fig 5A)
####

kos = read.csv("kos_DE.abs.tab",header = T)
nfams = read.csv("nfams_DE.abs.txt", header = T)
samples = read.table("metadata.tsv",header = T)
samples = samples[,c(1,6,11)] #  CRC - CTR


# load abs
kos = kos %>% 
  gather(key = "Sample_ID", value = "abs", 2:ncol(kos))
nfams = nfams %>% 
  gather(key = "Sample_ID", value = "abs", 2:ncol(nfams))

# load pvals / fcs
pvals = read.table("nov_fams_kos.CRC.pvalues_block.f_fams.tab",header = T)
fc = read.table("nov_fams_kos.CRC.fold_change.f_fams.tab",header = T)
names(fc) = c("clusters","fc")
dataP = pvals %>% 
  left_join(fc,by = "clusters")
names(dataP) = c("cluster","pval","fc")

# join
data = rbind(kos,nfams)
data = merge(data,samples,by = "Sample_ID") %>% 
  mutate(origin = case_when(str_detect(cluster,"ko")~"KO",
                            !str_detect(cluster,"ko")~"nfam")) %>%
  left_join(dataP,by = "cluster") %>% 
  mutate(direction = case_when(fc >= 0~ "CTR",
                               fc < 0~ "CRC")) %>%
  group_by(direction,origin) %>% 
  group_by(cluster) %>% 
  mutate(n_zeros = sum(abs == 0))

data$pop_CRC = paste(data$Group,data$Study)

# create individual datasets for CRC and CTR
crc_data = data %>% 
  filter(Group == "CRC" & origin == "nfam") %>% 
  mutate(join = paste(cluster,Study)) %>% 
  group_by(join) %>% 
  mutate(mean_abs_per_cluster_per_group_per_pop = mean(abs)) %>% 
  dplyr::select(!c(Sample_ID,abs)) %>% 
  unique()

ctr_data = data %>% 
  filter(Group == "CTR"  & origin == "nfam") %>% 
  mutate(join = paste(cluster,Study)) %>% 
  group_by(join) %>% 
  mutate(mean_abs_per_cluster_per_group_per_pop = mean(abs)) %>% 
  dplyr::select(!c(Sample_ID,abs)) %>% 
  unique()

# join for calculating difference of abs 
d_plot = merge(crc_data,ctr_data,by = "join")
d_plot$diff_crc_ctr = d_plot$mean_abs_per_cluster_per_group_per_pop.x - d_plot$mean_abs_per_cluster_per_group_per_pop.y

# calculate mean diff of abs per pop for ordering plot
d_plot = d_plot %>% 
  group_by(cluster.x) %>% 
  mutate(mean_diff = mean(diff_crc_ctr))
d_plot$cluster = fct_reorder(d_plot$cluster.x, desc(d_plot$diff_crc_ctr))

# rearrange factor values
d_plot = d_plot %>% 
  mutate(Study = case_when(Study.x == "AT.CRC"~"AT",
                           Study.x == "CN.CRC"~"CN",
                           Study.x == "DE.CRC"~"DE",
                           Study.x == "FR.CRC"~"FR",
                           Study.x == "US.CRC"~"US"))

d_plot$Direction = d_plot$direction.x

# plot
p1 = ggplot(d_plot,aes(y = cluster, x =  diff_crc_ctr))+
  geom_point(aes(shape = Study, color = Direction))+
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
  dplyr::select(cluster,pval.x,direction.x) %>% 
  unique()

p2 = ggplot(d_prev)+
  geom_bar(aes(y = pval.x, x = cluster,fill = direction.x),stat = "identity")+
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
  dplyr::select(cluster,n_zeros.x,direction.x) %>% 
  unique()

p3 = ggplot(d_prev)+
  geom_bar(aes(y = (575-n_zeros.x) / 575, x = cluster,fill = direction.x),stat = "identity")+
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

####
# CRC AUC comparison (Fig 5C)
####

data = read.table("linear_regression.70Train_30Test.lambda.1se.tab",header = T,sep = ",")
data$origin = factor(data$origin,levels = c("KOs","nfams","KOs + nfams"))

data = data %>% 
  mutate(origin_ = case_when(origin == "nfams"~"FESNov\n fams",
                            origin == "KOs + nfams"~"KOs + \nFESNov fams",
                            origin == "KOs"~"KOs"))

data$origin_ = factor(data$origin_,levels = c("KOs", "FESNov\n fams","KOs + \nFESNov fams"))
p1 = ggplot(data,aes(x = origin_,y = auc))+
  geom_boxplot()+
  theme_classic()+
  geom_signif(
    comparisons = list(c("KOs", "KOs + \nFESNov fams")),
    map_signif_level = TRUE
  ) +
  theme(text=element_text(size=25),
           axis.title=element_text(size=25),
  )+
  ylab ("AUC")+
  xlab("Dataset")









