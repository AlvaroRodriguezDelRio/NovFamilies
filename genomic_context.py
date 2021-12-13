import sys
import re

def report_neighs(gene2scaff,contig,beg,window_to_report):
        target_gene = gene2scaff[contig][beg]
        target_gene_name = target_gene[0] + "|" + beg + "|" + target_gene[1] + "|" + target_gene[2]
        scaff_starts_sorted = list(map(str,sorted(map(int,gene2scaff[contig].keys()))))
        position_gene = scaff_starts_sorted.index(beg)
        prev_position = position_gene - window_to_report
        while (prev_position < 0):
                prev_position += 1

        post_position = position_gene + window_to_report
        while post_position > (len(scaff_starts_sorted)-1):# last pos in array is len -1
                post_position -= 1

        to_report = []
        for pos in range(prev_position,post_position + 1):
                beg = scaff_starts_sorted[pos]
                gene_info = gene2scaff[contig][beg]
                gene_n = gene_info[0]
                end = gene_info[1]
                strand = gene_info[2]
                gene_n = gene_n.replace(">","")
                neigh_name = gene_n + "|" + beg + "|" + end + "|" + strand
                to_report.append(neigh_name)
        print ("%s\t%s" %(target_gene_name.replace('>',''),','.join(to_report)))

gene2scaff = dict()
for line in open(sys.argv[1]):
        if re.search(">",line):
                gene = line.split(" ")[0]
                contig = "_".join(gene.split("_")[:-1])
                beg = int(line.split("#")[1].replace(" ",""))
                end = line.split("#")[2].replace(" ","")
                strand = line.split("#")[3].replace(" ","")

                if contig not in gene2scaff:
                        gene2scaff[contig] = dict()
                gene2scaff[contig][beg] = [gene,end,strand]


for contig in gene2scaff:
        scaff_starts_sorted = sorted(gene2scaff[contig].keys())
        neighs_array = []
        for gene_pos in scaff_starts_sorted:
                gene = gene2scaff[contig][gene_pos]
                target_gene = gene2scaff[contig][gene_pos]
                target_gene_name = target_gene[0] + "|" + str(gene_pos) + "|" + target_gene[1] + "|" + target_gene[2]
                neighs_array.append(target_gene_name)
        print ('\t'.join([contig,','.join(neighs_array)]))
