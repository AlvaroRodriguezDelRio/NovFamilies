import json
import sys
from itertools import groupby
from collections import defaultdict

def readlines(file_):
    for line in open(file_):
        yield list(map(str.strip,line.split('\t')))


def calculate_number_neg_between(pos,neg_array,positions_annot):
    number_neg = 0
    pos = int(pos)
    positions_annot = [x for x in positions_annot.split(',')]

    if int(pos)>0:
        for i in range(1,pos):
            for j in neg_array[str(i)]:
                if str(j) in positions_annot:
                    number_neg += 1
    elif pos <0:
        for i in range(pos+1,0):
            for j in neg_array[str(i)]:
                if str(j) in positions_annot:
                    number_neg += 1
    return number_neg

with open("score_per_pos.strand.unknown.tab") as file1:
    for k, g in groupby(file1, key=lambda x:x.split('\t')[0]):
        fam2info = defaultdict(lambda:[])
        neg_strand_pos = defaultdict(lambda:set())
        lines = []
        for line in g:
            if len(line.split('\t')) < 3 or line == '':
                continue
            fam,pos,db,annot,number_genes_with_annot,number_genes_per_pos,nseqs_fam,positions_annot,positions_annot_contrary,num_contrary_strand = list(map(str.strip,line.split('\t')))

            # gather contrary strand coordinates
            if int (num_contrary_strand) > 0:
                for i in positions_annot_contrary.split(','):
                    neg_strand_pos[pos].add(int(i))

            if db != 'unknown':
                lines.append(line)


        fam2info['fam'] = fam
        added_cogs = set()
        for line in lines:
            fam,pos,db,annot,number_genes_with_annot,number_genes_per_pos,nseqs_fam,positions_annot,positions_annot_contrary,num_contrary_strand = list(map(str.strip,line.split('\t')))
            if annot.split('@')[0]+pos not in added_cogs:
                added_cogs.add(annot.split('@')[0] + pos)
                score = int(number_genes_with_annot)/int(nseqs_fam)
                number_neg_between = calculate_number_neg_between(pos,neg_strand_pos,positions_annot)
                number_genes_with_annot = len(positions_annot.split(','))

                cog_dict = {'n':annot.split('@')[0],'score':float(score),"mean_num_in_opposite_strand":float(int(num_contrary_strand)/int(number_genes_with_annot)),"mean_num_pos_opposite_strand_between":float(int(number_neg_between)/int(number_genes_with_annot)),"pos":int(pos)}
                fam2info[db].append(cog_dict)
        sys.stdout.write(json.dumps(fam2info, separators=(',',':'))+"\n")
