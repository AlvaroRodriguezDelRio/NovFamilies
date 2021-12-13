
#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pickle5

famdata = pickle5.load(open('/data/jhc/cold/MAGs/novel_fams-v2/build_fam_table/families_gtdb202.pkl', 'rb'))
famdata.reset_index(inplace=True)
famdata.rename(columns={'index': 'fam'}, inplace=True)

valid_fams = set([line.strip() for line in open("filtered_families.RNAcode.noPfamAcov.BUSTED.noPVOGs.noPfamB.noRefSeq_blastx.txt")])

def load_lineage_info():
    refseq_rep_genomes = pd.read_csv('reps_no_from_metagenomes.txt',
                                    sep='\t', names=['genome', 'rep', 'lineage', 'src'])
    refseq_rep_genomes[['d','p', 'c', 'o', 'f', 'g', 's']] = refseq_rep_genomes['lineage'].str.split(';', expand=True)

    genome2taxonomy =  pd.read_csv('/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab',
                                    sep='\t', names=['genome', 'lineage'])

    genome2taxonomy[['d','p', 'c', 'o', 'f', 'g', 's']] = genome2taxonomy['lineage'].str.split(';', expand=True)

    _tempdfs = []
    _tempdfs2 = []
    for rank in 'dpcofgs':
        rank_df = genome2taxonomy.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                    n_taxa=('lineage', 'nunique'),
                                                )
        rank_df2 = refseq_rep_genomes.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                        n_taxa=('lineage', 'nunique'),
                                                    )

        _tempdfs.append(rank_df)
        _tempdfs2.append(rank_df2)

    taxa2totalgenomes = pd.concat(_tempdfs)
    taxa2ngenomes = taxa2totalgenomes.to_dict()['n_genomes']
    taxa2refseq_genomes = pd.concat(_tempdfs2).to_dict()['n_genomes']

    return taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes
taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes = load_lineage_info()


def get_sinapomorphies(data, min_spec=0.9, min_cov=0.9, min_genomes=10):
    by_phylum = []
    visited_fams = set()
    for rank in "dpcofgs":
        fam_selection = data[(data["%s_specificity" %rank] >= min_spec) &                             (data["%s_coverage" %rank] >= min_cov) &
 (data["n_genomes"]>min_genomes) &                             (data["fam"].isin(valid_fams)) &                             (~data["fam"].isin(visited_fams))
                             ]
        visited_fams.update(fam_selection['fam'].to_list())
        print(len(visited_fams), len(fam_selection))
        fam_selection.reset_index(inplace=True)
        sel = fam_selection.groupby('%s_mostcommon'%rank, as_index=False).agg(n_fams = ('fam', 'count'),
                                                            avg_score = ('%s_score'%rank, 'mean'),
                                                            avg_specificity = ('%s_specificity'%rank, 'mean'),
                                                            avg_coverage = ('%s_coverage'%rank, 'mean'),
                                                            avg_genomes = ('n_genomes', 'mean'),
                                                            avg_taxa = ('n_taxa', 'mean'),

                                                            avg_duprate = ('duprate', 'mean'),
                                                            avg_emapper_hits = ('emapper_hits', 'mean'),
                                                            linage_sample = ('lineages', lambda x: str(list(x)[0])[:1000]),
                                                            fams = ('fam', set))

        if not sel.empty:
             sel[['rank', 'name']] = sel['%s_mostcommon'%rank].str.split('__', 1, expand=True)
             sel['total_genomes'] = sel['%s_mostcommon'%rank].apply(lambda x: taxa2ngenomes[x])
             sel['refseq_genomes'] = sel['%s_mostcommon'%rank].apply(lambda x: taxa2refseq_genomes.get(x, 0))

        by_phylum.append(sel)

    return pd.concat(by_phylum, axis=0)


cols = ['rank', 'name', 'n_fams', 'avg_score', 'avg_specificity', 'avg_coverage',
        'avg_genomes', 'total_genomes', 'refseq_genomes', 'avg_taxa', 'avg_duprate', 'avg_emapper_hits', 'fams', 'linage_sample']

sinapos = get_sinapomorphies(famdata[famdata['emapper_hits']==0], min_spec=1, min_cov=0.9, min_genomes=10)

famcounter = set()
for x in sinapos['fams']:
    famcounter.update(x)
print(len(famcounter))
sinapos[cols].sort_values('avg_taxa', ascending=False).head(50)
sinapos[cols].to_excel('rep_fams_with_all_filters.sp1.new.tab.xls')
