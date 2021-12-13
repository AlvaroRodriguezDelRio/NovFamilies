
from ete3 import Tree
import random
from collections import Counter,defaultdict
import re
import math
import pandas as pd
import numpy as np


def load_lineage_info():
    metagenomes = pd.read_csv('/data/jhc/freezer/public/GTDB/GTDB_rev95/reps_from_metagenomes_or_singlecell.txt',
                                    sep='\t', names=['genome', 'rep', 'lineage', 'src'])
    metagenomes[['d','p', 'c', 'o', 'f', 'g', 's']] = metagenomes['lineage'].str.split(';', expand=True)

    no_metagenomes =  pd.read_csv('/data/jhc/freezer/public/GTDB/GTDB_rev95/reps_no_from_metagenomes.txt',
                                    sep='\t', names=['genome', 'rep', 'lineage', 'src'])
    no_metagenomes[['d','p', 'c', 'o', 'f', 'g', 's']] = no_metagenomes['lineage'].str.split(';', expand=True)

    genome2taxonomy =  pd.read_csv('/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab',
                                    sep='\t', names=['genome', 'lineage'])
    genome2taxonomy[['d','p', 'c', 'o', 'f', 'g', 's']] = genome2taxonomy['lineage'].str.split(';', expand=True)

    _tempdfs = []
    _tempdfs2 = []
    _tempdfs3 = []
    for rank in 'dpcofgs':
        rank_df = genome2taxonomy.groupby(rank).agg(n_genomes=('genome', 'nunique'), n_taxa=('lineage', 'nunique'))
        rank_df2 = no_metagenomes.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                        n_taxa=('lineage', 'nunique'))
        rank_df3 = metagenomes.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                        n_taxa=('lineage', 'nunique'))


        _tempdfs.append(rank_df)
        _tempdfs2.append(rank_df2)
        _tempdfs3.append(rank_df3)

    taxa2totalgenomes = pd.concat(_tempdfs).to_dict()
    taxa2ngenomes = taxa2totalgenomes['n_genomes']
    taxa2refseq_genomes = pd.concat(_tempdfs2).to_dict()['n_genomes']
    taxa2metag_genomes = pd.concat(_tempdfs3).to_dict()['n_genomes']

    return taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes,taxa2metag_genomes


# load lca per family
def load_annots(filt_fams):
    fam2tax = defaultdict(lambda:set())

    # take families per genome and save genome tax
    genome2tax = {}
    for line in open("/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab"):
        genome,tax = list(map(str.strip,line.split('\t')))
        genome = genome.replace('.oneLine','')
        genome2tax[genome] = tax
        genome2tax[genome.split('.')[0]] = tax


    genome2fam = defaultdict(lambda:set())
    for line in open("../../../NoISO-fixed.3sp.taxonomy_by_member.tsv"):
        line = line.rstrip()
        fam = line.split('\t')[0]
        member = line.split('\t')[5]
        genome_member = member.split('@')[1]
        if fam not in filt_fams:
            continue
        tax = line.split('\t')[-1]
        genome2fam[genome_member].add(fam)

    # get list of number of families per genomes
    tax2number = defaultdict(lambda:[])
    for genome,fams in genome2fam.items():
        for t in genome2tax[genome].split(';'):
            tax2number[t].append(len(list(fams)))

    # calculate mean number of families per genome
    node2number = Counter()
    for tax,numbers in tax2number.items():
        #print (tax,numbers)
        node2number[tax] = np.mean(numbers)

    return fam2tax,node2number

# load genome id - species id file for renaming leaves in the gtdb file
def load_all_annots():
      genomeid2annot = {}
      name2lineage = defaultdict(lambda:set())

      # load for bacteria
      for line in open("bac120_taxonomy_r202.tsv"):
            line = line.rstrip()
            genome,annot = line.split('\t')
            genomeid2annot[genome] = annot
            for tax in annot.split(';'):
                name2lineage[tax].add(annot)

      # load for arhcaea
      for line in open("ar122_taxonomy_r202.tsv"):
            line = line.rstrip()
            genome,annot = line.split('\t')
            genomeid2annot[genome] = annot
            for tax in annot.split(';'):
                name2lineage[tax].add(annot)

      return genomeid2annot,name2lineage


######################## start ################################


# load final set of fams
filt_fams = set()
for line in open("/data/jhc/cold/MAGs/novel_fams-v2/clustering/filtering/filtered_families.RNAcode.noPfamAcov.BUSTED.noPVOGs.noPfamB.noRefSeq_blastx.txt"):
    line = line.rstrip()
    filt_fams.add(line)

# load tax annot per family (lca) and count number of families per lin.
fam2tax,node2number = load_annots(filt_fams)
# load gtdb genome 2 taxonomy
# load all lineage for all taxa
genomeid2annot,name2lineage = load_all_annots()

# load gtdb trees and combine
tree_bact = Tree("bac120_r202.tree",format = 1,quoted_node_names = True)
tree = Tree("ar122_r202.tree",format = 1,quoted_node_names = True)

root_node_bac = next(n.get_tree_root() for n in tree_bact.traverse())
root_node_arch = next(n.get_tree_root() for n in tree.traverse())
root_node_arch.add_child(root_node_bac)

# rename nodes if they have annot to the order level
for n in tree.traverse():
    if re.search('__',str(n.name)):
        for poss_name in str(n.name).split(':')[-1].split(';'):
            if re.search('o__',poss_name):
                n.name = poss_name.strip()

# quit leaves not to the order level
ex = False
while not ex:
    ex = True
    for n in tree.iter_leaves():
        if not re.search('o__',n.name):
            n.detach()
            ex = False


# get lineage info (number cultivated genomes per branch)
# calculated for each node, can be used for each collapsed level
taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes,taxa2metag_genomes = load_lineage_info()


###
## print info
###

# print tree
tree_out = open("tree_new.nw","w")
tree_out.write(tree.write(format = 1))
tree_out.close()

# print number of genomes table
n_file = open("fraction_uncultivated_cultivated_new.tab",'w')
for node in tree.iter_leaves():
    prop = 0
    if node.name in taxa2metag_genomes and node.name in taxa2refseq_genomes:
        prop = taxa2metag_genomes[node.name]/(taxa2metag_genomes[node.name] + taxa2refseq_genomes[node.name])
    elif node.name in taxa2metag_genomes and node.name not in taxa2refseq_genomes:
        prop = 1
    n_file.write(node.name + '\t' + str(prop) + '\n')
n_file.close()


# print bars table
n_file = open("leave_abundance_new.tab",'w')
for node in tree.iter_leaves():
    number = node2number[node.name]
    n_file.write(node.name + "\t" + str(number) + '\n')
n_file.close()

# print bats table, abundance
n_file = open("leave_abundance.raw.new.tab",'w')
for node in tree.iter_leaves():
    number = node2number[node.name]
    n_file.write(node.name + "\t" + str(number) + '\n')
n_file.close()
