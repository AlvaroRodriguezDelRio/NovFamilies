from ete3 import (Tree, add_face_to_node, RectFace, TextFace, TreeStyle, random_color, StackedBarFace,CircleFace,AttrFace)
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

    genome2taxonomy =  pd.read_csv('/data/jhc/cold/MAGs/novel_fams-v2/build_fam_table/genome2taxonomy',
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



# collapse tree
def collapse_tree(tree,name2lineage,level_trunc):
    prev_dict = {'s':'g',
                 'g':'f',
                 'f':'o',
                 'o':'c',
                 'c':'p',
                 'p':'d'}

    level_t_index = ['d','p','c','o','f','g','s'].index(level_trunc)
    levels2trunc = ['d','p','c','o','f','g','s'][level_t_index+1:]
    levels_prev_trunc = ['d','p','c','o','f','g','s'][:level_t_index+1]
    # eliminate 'first layer', before level of interest
    exit_l = False
    while not exit_l:
        exit_l = True
        for i,node in enumerate(tree.traverse()):
            rank = node.name.split('__')[0]
            if rank in levels2trunc:
                exit_l = False
                if node.up.name.split('__')[0] != prev_dict[rank]: # if this is the only representative of the prev. rank, keep in the tree changing name
                    lineage = list(name2lineage[node.name])[0]
                    new_name = create_dict_annot(lineage)[prev_dict[rank]]
                    node.name = prev_dict[rank] + '__' + new_name
                else:
                    node.detach()
            elif rank == level_trunc: # avoid multiple branches with same names at the level you want to collapse
                if node.up.name == node.name:
                    exit_l = False
                    node.detach()
    return (tree)


# load lca per family
def load_annots(filt_fams):
    fam2tax = defaultdict(lambda:set())

    # take families per genome and save genome tax
    genome2tax = {}
    genome2fam = defaultdict(lambda:set())
    for line in open("NoISO-fixed.3sp.taxonomy_by_member.tsv"):
        line = line.rstrip()
        fam = line.split('\t')[0]
        member = line.split('\t')[5]
        genome_member = member.split('@')[1]
        if fam not in filt_fams:
            continue
        tax = line.split('\t')[-1]
        genome2tax[genome_member] = tax
        genome2fam[genome_member].add(fam)

    # get list of number of families per genomes
    tax2number = defaultdict(lambda:[])
    for genome,fams in genome2fam.items():
        for t in genome2tax[genome].split(';'):
            tax2number[t].append(len(list(fams)))

    # calculate mean number of families per genome
    node2number = Counter()
    for tax,numbers in tax2number.items():
        print (tax,numbers)
        node2number[tax] = np.mean(numbers)

    return fam2tax,node2number

# load genome id - species id file for renaming leaves in the gtdb file
def load_all_annots():
      genomeid2annot = {}
      name2lineage = defaultdict(lambda:set())

      # load for bacteria
      for line in open("/data/jhc/freezer/public/GTDB/GTDB_rev95/bac120_taxonomy_r95.tsv"):
            line = line.rstrip()
            genome,annot = line.split('\t')
            genomeid2annot[genome] = annot
            for tax in annot.split(';'):
                name2lineage[tax].add(annot)

      # load for arhcaea
      for line in open("/data/jhc/freezer/public/GTDB/GTDB_rev95/ar122_taxonomy_r95.tsv"):
            line = line.rstrip()
            genome,annot = line.split('\t')
            genomeid2annot[genome] = annot
            for tax in annot.split(';'):
                name2lineage[tax].add(annot)

      return genomeid2annot,name2lineage

# format tree node names so that they are the tax name in gtdb
def create_dict_annot(tax_line):
    level2tax = {}
    for annot in tax_line.split(';'):
        if annot != "" and annot != 'NULL' and re.search("_",annot):
            level,name = annot.split('__')
            level2tax[level] = name
    return level2tax


def format_tree(tree,genomeid2annot):

    def  fet_most_sp_tax(tax_dict):
        prev_tax = ""
        prev_r = ""
        for x in ['d','p','c','o','f','g','s']:
            if x in tax_dict:
                prev_tax = tax_dict[x]
                prev_r = x
            else:
                return prev_tax,prev_r
        return prev_tax,prev_r

    def get_lca_per_family(genome_annots):
        lca = ""
        lca_pos = 0
        for annot in genome_annots:
            if lca == "":
                lca = create_dict_annot(annot)
                continue

            annot_dict = create_dict_annot(annot)
            for key,item in annot_dict.items():
                if key in lca:
                    if lca[key] != item:
                        del (lca[key])
        return lca


    # get lca per node
    branch2level = {}
    node2name = {}
    name2complete_name = {}
    name2p = {}
    for i,node in enumerate(tree.traverse()):
        annotations = set()
        for leaf in node.iter_leaves():
            taxonomy = genomeid2annot[leaf.name]
            annotations.add(taxonomy)
        lca = get_lca_per_family(annotations)
        lca_name,rank = fet_most_sp_tax(lca)
        branch2level[rank + '__' + lca_name] = rank
        node2name[i] = rank + '__' + lca_name
        node.name =  rank + '__' + lca_name
        if 'p' in lca:
            name2p[rank + '__' + lca_name] = 'p__' + lca['p']

    return tree,branch2level

def calculate_number_nodes_before(node):
    number = 0
    while node.up:
        number += 1
        node = node.up
    return number

######################## start ################################

colors = random_color(num=1000, s=0.3, l=0.9)

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
tree = Tree("bac120_r95.tree",format = 1,quoted_node_names = True)
tree_arch = Tree("ar122_r95.tree",format = 1,quoted_node_names = True)

# merge trees
root_node_bac = next(n.get_tree_root() for n in tree.traverse())
root_node_arch = next(n.get_tree_root() for n in tree_arch.traverse())
root_node_bac.add_child(root_node_arch)

# format tree names & get dict of tax level for each gtdb name
tree,branch2level = format_tree(tree,genomeid2annot)

# collapse tree
tree = collapse_tree(tree,name2lineage,'o')

# get lineage info (number cultivated genomes per branch)
# calculated for each node, can be used for each collapsed level
taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes,taxa2metag_genomes = load_lineage_info()

###
## print info
###

# print tree
tree_out = open("tree.nw","w")
tree_out.write(tree.write(format = 1))
tree_out.close()

# print number of genomes table
n_file = open("fraction_uncultivated_cultivated.tab",'w')
for node in tree.iter_leaves():
    prop = 0
    if node.name in taxa2metag_genomes and node.name in taxa2refseq_genomes:
        prop = taxa2metag_genomes[node.name]/(taxa2metag_genomes[node.name] + taxa2refseq_genomes[node.name])
    elif node.name in taxa2metag_genomes and node.name not in taxa2refseq_genomes:
        prop = 1
    n_file.write(node.name + '\t' + str(prop) + '\n')
n_file.close()

# print bars table
n_file = open("leave_abundance.tab",'w')
for node in tree.iter_leaves():
    number = node2number[node.name]
    n_file.write(node.name + "\t" + str(number) + '\n')
n_file.close()
