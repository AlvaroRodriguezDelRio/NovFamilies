library(data.table)
library(caret)
library(dplyr)
library(coin)
library(gtools)
library(tidyr)


# load abundance data per novel family / KO
setwd("~/analysis/nov_fams_v2/CRC/DE_analysis/")
data_n = read.csv("abs_per_fam_all.f_fams.h.tab", header = T, stringsAsFactors=F)
data_k = read.csv("kos_abs.tab")

# merge kos and nfam data (same sample order in columns)
data = rbind(data_n,data_k)

samples = read.table("metadata.tsv.f",header = T)

# include colonoscopy into study for correcting in pvalue function for calculating p-values
samples <- samples %>%
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN.CRC', Study,
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))

samples$Study = samples$block
samples = samples[,c(1,6,11)]

# transpose dataframe for having novel families / kos as rows
names = data[,1]
data = data[,-1]
colnames= names(data)
data = data.frame(transpose(data))
names(data) = names
data$sample = colnames

# merge abundance data and metadata
names(samples) = c("sample","Group","Study")
data_m = merge(data, samples, by = "sample")
data_m$Study = as.factor(data_m$Study)
data_m$Group = as.factor(data_m$Group)

#  Wilcoxon test blocked by population and colonoscopy for each novel family / KO
nfams = ncol(data_m)-3
mean_comparison = function(x) p.adjust(pvalue(wilcox_test(x ~ Group | Study, data=data_m)),method = 'fdr',n = nfams)
p_values = as.matrix(lapply(data_m[2:(ncol(data_m)-2)],mean_comparison))

# generalized fold change (Wirbel et al., 2019) calculation for each novel family / KO 
foldch_res = data.frame(cluster = character(),foldc = numeric(),stringsAsFactors=FALSE)
for (i in c(2:(ncol(data_m)-2))){
  cluster_abs = data_m[,c(i,(ncol(data_m)-1))]
  abs_crc = filter(cluster_abs,Group == "CRC")[,1]
  abs_ctr = filter(cluster_abs,Group == "CTR")[,1]
  qcrc <- quantile(log10(abs_crc +1e-20), probs=seq(.1, .9, .05))
  qctr <- quantile(log10(abs_ctr+1e-20), probs=seq(.1, .9, .05))
  change = sum(qctr - qcrc)/length(qctr)
  new_row = c(as.character(names(cluster_abs)[1]),as.numeric(change))
  foldch_res[(nrow(foldch_res) + 1), ] = new_row
}
