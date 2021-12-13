import pandas as pd
from collections import Counter
import pickle
import json
import sys

class JsonEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Counter):
            return dict(obj)
        elif isinstance(obj, set):
            return list(obj)
        else:
            return super(JsonEncoder, self).default(obj)


def tsv_to_feather():
    header_names = ['fam', 'famsize', 'famspecies', 'src', 'genome', 'member', 'lineage']
    family_data = pd.read_csv('NoISO-fixed.3sp.taxonomy_by_member.tsv', sep='\t',
        header=0,
        names=header_names)

    family_data.to_feather("NoISO-fixed.3sp.taxonomy_by_member.feather")
    return family_data

def load_lineage_info():

    refseq_genomes = pd.read_csv('genome2gtdb202_refseq.tsv',
                                    sep='\t', names=['genome', 'lineage'])
    refseq_genomes[['d','p', 'c', 'o', 'f', 'g', 's']] = refseq_genomes['lineage'].str.split(';', expand=True)

    all_genomes = pd.read_csv('/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab',
                                    sep='\t', names=['genome', 'lineage'])
    all_genomes[['d','p', 'c', 'o', 'f', 'g', 's']] = all_genomes['lineage'].str.split(';', expand=True)

    _tempall = []
    _temprefseq = []
    for rank in 'dpcofgs':
        ranks_all = all_genomes.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                    n_taxa=('lineage', 'nunique'),
                                                )
        ranks_refseq = refseq_genomes.groupby(rank).agg(n_genomes=('genome', 'nunique'),
                                                        n_taxa=('lineage', 'nunique'),
                                                    )

        _tempall.append(ranks_all)
        _temprefseq.append(ranks_refseq)

    taxa2totalgenomes = pd.concat(_tempall).to_dict()
    taxa2ngenomes = taxa2totalgenomes['n_genomes']
    taxa2refseq_genomes = pd.concat(_temprefseq).to_dict()['n_genomes']
    return taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes

def build_family_table():
    # GTDBiso@GB_GCA_003141675@PLND01000034.1_2@d__Bacteria|p__Chloroflexota  4       4       GEM     3300020139_31   GEM@3300020139_31@3300020139_31_00357@d__Bacteria|p__Bipolaricaulota    d__Bacteria;p__Bipol
    genome2lineage = {}
    for line in open('/data/jhc/cold/MAGs/combined/genome2taxonomy.r202.tab'):
        _g, _l = map(str.strip, line.split('\t'))
        if _g.startswith('GB_') or _g.startswith('RS_'):
            _g = _g.split('.')[0]
        genome2lineage[_g] = _l


    emapper_hits = set([line.split()[0] for line in open('microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.all_members.mmseqs.emapper.seed_orthologs')])

    all_fams = {}
    for ln, line in enumerate(open('NoISO-fixed.3sp.taxonomy_by_member.tsv')):
        if ln > 0 and ln % 100000 == 0:
            print(ln, file=sys.stderr)')

        famname, _, famspecies, src, genome, member, lineage95 = map(str.strip, line.split('\t'))


        src, genome, gene, _ = member.split('@')

        if genome == 'MGYG-HGUT-01456.oneLine':
            genome = 'MGYG-HGUT-01456'

        lineage = genome2lineage[genome]
        d, p, c, o, f, g, s = lineage.split(';')


        # Initialize family
        if famname not in all_fams:
            all_fams[famname] = {'nspecies': int(famspecies),
                            'members': set(),
                            'genomes': set(),
                            'emapper_hits': 0,
                            'sources': Counter(),
                            'lineages': Counter(),
                            'd': Counter(),
                            'p': Counter(),
                            'c': Counter(),
                            'o': Counter(),
                            'f': Counter(),
                            'g': Counter(),
                            's': Counter(),
                            'd_taxa': Counter(),
                            'p_taxa': Counter(),
                            'c_taxa': Counter(),
                            'o_taxa': Counter(),
                            'f_taxa': Counter(),
                            'g_taxa': Counter(),
                            's_taxa': Counter(),
                            }

        # Update family info
        # In which taxa terms is observed?
        if genome not in all_fams[famname]['genomes']:
            all_fams[famname]['d'].update([d])
            all_fams[famname]['p'].update([p])
            all_fams[famname]['c'].update([c])
            all_fams[famname]['o'].update([o])
            all_fams[famname]['f'].update([f])
            all_fams[famname]['g'].update([g])
            all_fams[famname]['s'].update([s])
            all_fams[famname]['genomes'].add(genome)

        if lineage not in all_fams[famname]['lineages']:
            all_fams[famname]['d_taxa'].update([d])
            all_fams[famname]['p_taxa'].update([p])
            all_fams[famname]['c_taxa'].update([c])
            all_fams[famname]['o_taxa'].update([o])
            all_fams[famname]['f_taxa'].update([f])
            all_fams[famname]['g_taxa'].update([g])
            all_fams[famname]['s_taxa'].update([s])
            all_fams[famname]['lineages'].update([lineage])

        all_fams[famname]['sources'].update([src])
        all_fams[famname]['members'].add(member)
        if member in emapper_hits:
            all_fams[famname]['emapper_hits'] += 1

    for fam, famdata in all_fams.items():
        famdata['n_genomes'] = len(famdata['genomes'])
        famdata['n_members'] = len(famdata['members'])
        famdata['n_taxa'] = len(famdata['lineages'])
        famdata['duprate'] = len(famdata['members']) / len(famdata['genomes'])

        #famdata['besthit'] = # emapper
        #famdata['refseq'] = # refseq
        clade_counter = []
        for rank in 'dpcofgs':
            taxa_name, taxa_occ = famdata[rank].most_common()[0]

            famdata['%s_mostcommon'%rank] = taxa_name

            # number of genomes from that taxa where family was observed
            famdata['n_%s'%rank] = taxa_occ

            cov = taxa_occ / taxa2ngenomes.get(taxa_name, 1)
            famdata['%s_coverage' %rank] = cov

            spec = taxa_occ / len(famdata['genomes'])
            famdata['%s_specificity' %rank] = spec

            famdata['%s_score'%rank] = min(cov, spec)

            famdata['%s_total_genomes'%rank] = taxa2ngenomes.get(taxa_name, 1)
            famdata['%s_refseq_genomes'%rank] = taxa2refseq_genomes.get(taxa_name, 0)

            #famdata['%s_terms' %rank] = list(famdata[rank].keys())
            for term, count in famdata["%s_taxa"%rank].items():
                if term is None:
                    continue
                term_spec = count / famdata['n_taxa']
                term_cov = count / taxa2totalgenomes['n_taxa'][term]
                clade_counter.append({'term': term,
                                      'ntaxa': count,
                                      'specificity': term_spec,
                                      'coverage': term_cov,
                                      'score': min(term_cov, term_spec),
                                      'ngenomes': famdata[rank][term],
                                      'rank': rank})

        famdata['clade_counter'] = clade_counter
        famdata['name'] = fam
        famdata['code'] = fam2code[fam]

        print(json.dumps(famdata, cls=JsonEncoder))


taxa2totalgenomes, taxa2ngenomes, taxa2refseq_genomes = load_lineage_info()
fam2code = {i.rstrip().split('\t')[0]:i.rstrip().split('\t')[1] for i in open('/data/jhc/cold/MAGs/novel_fams-v2/clustering/microbial_genomes-v1.clustering.folded.parsed.NoISO.3sp.fam2codes.tsv')}
build_family_table()
