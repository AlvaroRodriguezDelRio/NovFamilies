from ete3 import (Tree, add_face_to_node, RectFace, TextFace, TreeStyle, random_color, StackedBarFace,CircleFace,AttrFace,SeqMotifFace,BarChartFace)
from ete3.treeview.faces import SequencePlotFace
import random
from collections import Counter,defaultdict
import re
import math
import pandas as pd
import json


def load_lineage_info():
    metagenomes = pd.read_csv('reps_from_metagenomes_or_singlecell.txt',
                                    sep='\t', names=['genome', 'rep', 'lineage', 'src'])
    metagenomes[['d','p', 'c', 'o', 'f', 'g', 's']] = metagenomes['lineage'].str.split(';', expand=True)

    no_metagenomes =  pd.read_csv('reps_no_from_metagenomes.txt',
                                    sep='\t', names=['genome', 'rep', 'lineage', 'src'])
    no_metagenomes[['d','p', 'c', 'o', 'f', 'g', 's']] = no_metagenomes['lineage'].str.split(';', expand=True)

    genome2taxonomy =  pd.read_csv('genome2taxonomy.r202.tab',
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


def create_dict_annot(tax_line):
    level2tax = {}
    for annot in tax_line.split(';'):
        if annot != "" and annot != 'NULL' and re.search("_",annot):
            level,name = annot.split('__')
            level2tax[level] = name
    return level2tax

# get random color for each biome for representation
def get_route2color(routes,colors):
      route2color = {}
      added = set()
      color = random.sample(colors, 22)
      for i,route in enumerate(routes):
        route2color[route] = color[i]
      return(route2color)


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

      # load for archaea
      for line in open("ar122_taxonomy_r202.tsv"):
          line = line.rstrip()
          genome,annot = line.split('\t')
          genomeid2annot[genome] = annot
          for tax in annot.split(';'):
              name2lineage[tax].add(annot)

      return genomeid2annot,name2lineage

# format tree node names so that they are the tax name in gtdb
def format_tree(tree,genomeid2annot):
    def create_dict_annot(tax_line):
        level2tax = {}
        for annot in tax_line.split(';'):
            if annot != "" and annot != 'NULL' and re.search("_",annot):
                level,name = annot.split('__')
                level2tax[level] = name
        return level2tax

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
    for i,node in enumerate(tree.traverse()):
        annotations = set()
        for leaf in node.iter_leaves():
            taxonomy = genomeid2annot[leaf.name]
            annotations.add(taxonomy)
        lca = get_lca_per_family(annotations)
        lca_name,rank = fet_most_sp_tax(lca)
        node.name = rank + '__' + lca_name

    return tree


# layout function
def layout(node):

    node.img_style['size'] = 0
    node.img_style['hz_line_width'] = 2
    node.img_style['vt_line_width'] = 2

    if not node.up:
        # do not show branch length of the root
        node.img_style['hz_line_width'] = 0

    aligned_column = 0

    if node.is_leaf():

        # bars for representing number of cultivated / uncultivated genomes
        bar_cult = ""
        if node.name in taxa2metag_genomes and node.name in taxa2refseq_genomes:# and node.name in taxa2totalgenomes:
            bar_cult = RectFace((taxa2metag_genomes[node.name]/(taxa2metag_genomes[node.name] + taxa2refseq_genomes[node.name]))*100,20, fgcolor='red', bgcolor='red')
        elif node.name in taxa2metag_genomes and node.name not in taxa2refseq_genomes: # only in uncultivated,ratio = 1
            bar_cult = RectFace(100 ,20, fgcolor='red', bgcolor='red')
        else:
            bar_cult = RectFace(0,20, fgcolor='red', bgcolor='red')
        add_face_to_node(bar_cult, node, column=aligned_column, position='aligned')
        aligned_column += 1

        # space between matrix and uncultivated
        white_bar = RectFace(100,20, fgcolor='white', bgcolor='white')
        add_face_to_node(white_bar, node, column=aligned_column, position='aligned')


        # kos matrix per leave
        for ko,number in taxa2kos[node.name].items():
              aligned_column += 1
              for route,col in route_especific_ko_order[ko].items():
                  color_route = route2color[route]

                  #circle_face = CircleFace(radius=number*3, color=color_route)
                  #circle_face.opacity = 1
                  #add_face_to_node(circle_face, node, column= col, position="aligned")

                  rectface = RectFace(10,10, fgcolor=color_route, bgcolor=color_route)
                  add_face_to_node(rectface, node, column=col + 1, position='aligned')


######################## start ################################

####
# load tree
###


# load gtdb genome 2 taxonomy
genomeid2annot,name2lineage = load_all_annots()

# load gtdb trees and combine
tree = Tree("bac120_r202.tree",format = 1,quoted_node_names = True)
tree_arch = Tree("ar122_r202.tree",format = 1,quoted_node_names = True)
root_node_bac = next(n.get_tree_root() for n in tree.traverse())
root_node_arch = next(n.get_tree_root() for n in tree_arch.traverse())
root_node_bac.add_child(root_node_arch)

# format tree names & get dict of tax level for each gtdb name
tree = format_tree(tree,genomeid2annot)

# collapse tree
collapse_level = 'o'
tree = collapse_tree(tree,name2lineage,collapse_level)


####
# load conserved kos
####

colors = random_color(num=22, s=0.6, l=0.5)
colors_bk = random_color(num=22, s=0.2, l=0.9)


# load taxonomies per family
fam2tax = defaultdict(lambda:set())
for line in open("taxonomy_by_member.tsv"):
    line = line.rstrip()
    fam = line.split('\t')[0]
    tax = line.split('\t')[-1]
    for t in tax.split(';'):
        if t.split('__')[-1] != '':
            fam2tax[fam].add(t)


# load filtered family list
filtered_families = set()
for line in open("ffams.txt"):
    line = line.rstrip()
    filtered_families.add(line)

code2fam = {}
for line in open('fam2codes.tsv'):
    fam,code = list(map(str.strip,line.split('\t')))
    code2fam[code] = fam

# load annotations per lineage
taxa2kos = defaultdict(lambda:Counter())
kos_present = set()
kos_counter = Counter()
taxa_with_kos = set()
for line in open('kpath_confidence_per_nfam.tab'):
    fam_code,annot,acc,param,desc =  list(map(str.strip,line.split('\t')))
    if float(acc) >= 0.9:
        fam = code2fam[fam_code]
        if fam in filtered_families:
            kos_counter[annot] += 1
            for tax in fam2tax[fam]:
                taxa2kos[tax][annot] += 1
                if re.search("__",tax):
                    if tax.split('__')[0] == collapse_level and tax.split('__')[1] != "":
                        kos_present.add(annot)

# load interesting kos
route_set = set()
route_especific_ko_order = defaultdict(lambda:{})
pos = 0
ko2enzyme_name = {}
for line in open('KEGG_maps_selected.tab'):
    route, ko, tr = list(map(str.strip,line.split('\t')))
    ko = ko.replace('map','')
    if ko in kos_present: # if any prediction in neigh
        pos += 1
        route_set.add(route)
        route_especific_ko_order[ko][route] = pos
        ko2enzyme_name[ko] = tr


# discard leaves in the tree with no kos assigned
exit_loop = False
while not exit_loop:
    exit_loop = True
    for node in tree.traverse():
        if node.is_leaf() and len(list(set(taxa2kos[node.name].keys()) & set(route_especific_ko_order.keys()))) == 0:
            node.detach()
            exit_loop = False

        # eliminate prev level leaves after removing leaves from the level we want
        if node.is_leaf() and node.name.split('__')[0] != collapse_level:
            node.detach()
            exit_loop = False


# load genomes to know the percentage of uncultivated genomes in each branch
taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes,taxa2metag_genomes = load_lineage_info()

#####
# layout
#####

route2color = get_route2color(route_set,colors)

tree.dist = 0
ts = TreeStyle()
ts.root_opening_factor = 0.9
ts.show_leaf_name = False#True

# create headers for heatmap and barplot with number of kos
ko2col_ene = {}
number_kos_not_present = 0
for ko in route_especific_ko_order:
    for route,pos in route_especific_ko_order[ko].items():
          color_route = route2color[route]
          txt_face = TextFace(ko,  fsize=9, bold=True)
          txt_face.rotation = 90
          txt_face.vt_align = 0
          ts.aligned_header.add_face(txt_face, column=pos + 1)# +2 for adding space between num genomes and matrix, delete otherwise (also in matrix coordinates in layout function)

          bar_h = max(int(math.log2(kos_counter[ko])*20),1)
          bar = RectFace(10,bar_h, fgcolor=color_route, bgcolor=color_route)
          bar.vt_align = 2

          ts.aligned_header.add_face(bar, column=pos + 1)


# create tree layout and export
ts.layout_fn = layout
ts.scale =  60
tree.render('tree.kpath_conf_09.pdf', tree_style=ts)
