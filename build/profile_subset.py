import pandas as pd, re, sys
base = '/lustre/home/vc362/sfari/D_redo_Nov25/'
gtf = base + '13_genetics/WholeTargeted_relaxedsqanti_10reads10samples_CPAT_CDS.gtf'
cls = base + '11_sqanti_final/relaxed_filter/sqantifiltered_monoexonicfiltered_10reads10samples_classification.txt'

g = pd.read_csv(gtf, sep='\t', comment='#', header=None, usecols=[2, 8])
g['tx'] = g[8].str.extract(r'transcript_id "([^"]+)"')[0]
coding = set(g.loc[g[2] == 'CDS', 'tx'])
ingtf = set(g['tx'])
c = pd.read_csv(cls, sep='\t', low_memory=False)
whole = [x for x in c.columns if re.fullmatch(r'Whole\d+', x)]
c['whole_reads'] = c[whole].sum(axis=1)
c['in_gtf'] = c['isoform'].isin(ingtf)
c['gtf_coding'] = c['isoform'].isin(coding)
c['pc_gene'] = c['associated_gene'].isin(set(open(base + '00_meta/protein_coding_ENSGs.txt').read().split()))
print('n whole cols', len(whole))
print(c.groupby('in_gtf').agg(n=('isoform', 'size'), reads=('whole_reads', 'sum')))
print(pd.crosstab(c['structural_category'], c['in_gtf']))
print(pd.crosstab(c['coding'], c['in_gtf']))
print(pd.crosstab(c['pc_gene'], c['in_gtf']))
print(pd.crosstab(c['monomulti'], c['in_gtf']))
print('gtf coding tx', len(coding))
a = c[c.associated_gene.str.startswith('ENSG00000141505')]
print(a[['isoform', 'structural_category', 'associated_transcript', 'coding', 'in_gtf', 'gtf_coding', 'whole_reads', 'CDS_genomic_start', 'CDS_genomic_end']].sort_values('whole_reads', ascending=False).head(25).to_string())
