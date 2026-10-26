"""Build per-gene long-read coverage files for the LRBrainCoverage app.

For one SQANTI filtered set (CPAT CDS GTF + SQANTI classification with per-sample counts):
  coverage(base, region, group) = sum of CPM of transcripts whose structure includes the base
                                  with that region label / sum of CPM of all transcripts of the gene
Region labels come from each transcript's own CPAT ORF: CDS, 5'UTR, 3'UTR, or non-coding
(transcripts with no CPAT ORF above the coding threshold). Only Whole transcriptome samples are used.

Usage: python build_coverage.py <gtf> <classification> <sample_ages.tsv> <mane_genes.tsv> <outdir>
"""
import sys, os, re
import numpy as np
import pandas as pd

gtf_path, cls_path, ages_path, mane_path, outdir = sys.argv[1:6]
os.makedirs(f'{outdir}/genes', exist_ok=True)

# sample groups (same groups as the original app)
ages = pd.read_csv(ages_path, sep='\t')
pre, post = ages.period == 'Prenatal', ages.period == 'Postnatal'
groups = {
    'total_exp': ages['sample'].notna(),
    'prenatal_exp': pre,
    'postnatal_exp': post,
    'prenatal1sttrimester_exp': pre & (ages.age < 14),
    'prenatal2ndtrimester_exp': pre & (ages.age >= 14) & (ages.age < 27),
    'prenatal3rdtrimester_exp': pre & (ages.age >= 27),
    'postnatalchild_exp': post & (ages.age < 18),
    'postnataladult_exp': post & (ages.age >= 18) & (ages.age < 60),
    'postnatalelderly_exp': post & (ages.age >= 60),
}
groups = {k: ages.loc[v, 'sample'].tolist() for k, v in groups.items()}
print({k: len(v) for k, v in groups.items()}, flush=True)

# GTF features
gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, usecols=[0, 2, 3, 4, 6, 8],
                  names=['chr', 'type', 'start', 'end', 'strand', 'attr'])
gtf = gtf[gtf.type.isin(['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'])]
gtf['tx'] = gtf.attr.str.extract(r'transcript_id "([^"]+)"')[0]
gtf = gtf.drop(columns='attr')
coding_tx = set(gtf.loc[gtf.type == 'CDS', 'tx'])
feat = pd.concat([
    gtf[gtf.tx.isin(coding_tx) & (gtf.type != 'exon')],
    gtf[~gtf.tx.isin(coding_tx) & (gtf.type == 'exon')],
])
feat['region'] = feat.type.map({'CDS': 'CDS', 'five_prime_UTR': "5'UTR", 'three_prime_UTR': "3'UTR", 'exon': 'non-coding'})

# check that CDS + UTR features tile the exons of coding transcripts
ln = (gtf.end - gtf.start + 1).groupby([gtf.tx, gtf.type == 'exon']).sum().unstack()
ln = ln.loc[ln.index.isin(coding_tx)]
print('coding tx', len(ln), 'exon length == CDS+UTR length:', (ln[True] == ln[False]).mean().round(4),
      'median diff', (ln[True] - ln[False]).median(), flush=True)

# counts -> CPM per sample, summed per group
whole = [s for s in ages['sample']]
cls = pd.read_csv(cls_path, sep='\t', usecols=['isoform', 'associated_gene', 'structural_category'] + whole)
cls = cls[cls.isoform.isin(set(feat.tx))].copy()
cls['gene_id'] = cls.associated_gene.str.replace(r'\.\d+$', '', regex=True)
cpm = cls[whole] / cls[whole].sum() * 1e6
w = pd.DataFrame({k: cpm[v].sum(axis=1) for k, v in groups.items()})
w.index = cls.isoform
cls['nreads'] = cls[whole].sum(axis=1)
print('transcripts with counts', len(cls), 'of', feat.tx.nunique(), flush=True)

# genes in the app: those with a MANE Select transcript
mane = pd.read_csv(mane_path, sep='\t')
cls = cls[cls.gene_id.isin(mane.gene_id)]
gene_of = cls.set_index('isoform').gene_id
feat = feat[feat.tx.isin(gene_of.index)].copy()
feat['gene_id'] = feat.tx.map(gene_of)
G = list(groups)
gene_tot = w.loc[cls.isoform].groupby(cls.gene_id.values).sum()

# sweep line: +w at start, -w at end+1, cumulative sum per gene x region
wf = w.loc[feat.tx].to_numpy()
ev = pd.concat([
    pd.DataFrame(wf, columns=G).assign(gene_id=feat.gene_id.values, region=feat.region.values, pos=feat.start.values),
    pd.DataFrame(-wf, columns=G).assign(gene_id=feat.gene_id.values, region=feat.region.values, pos=feat.end.values + 1),
])
ev = ev.groupby(['gene_id', 'region', 'pos'], sort=True)[G].sum().reset_index()
ev[G] = ev.groupby(['gene_id', 'region'])[G].cumsum()
ev['end'] = ev.groupby(['gene_id', 'region']).pos.shift(-1) - 1
ev = ev.dropna(subset=['end'])
ev[G] = ev[G].to_numpy() / gene_tot.loc[ev.gene_id, G].to_numpy()
ev[G] = ev[G].fillna(0).round(5).clip(lower=0) + 0.0  # 0/0 when a gene has no reads in a group; +0.0 drops -0.0
ev = ev[(ev[G] > 0).any(axis=1)].rename(columns={'pos': 'start'})
ev['end'] = ev.end.astype(int)

chrom = feat.groupby('gene_id').chr.first()
name = mane.set_index('gene_id').gene_name
for gid, d in ev.groupby('gene_id'):
    d = d.assign(chr=chrom[gid])[['chr', 'start', 'end', 'region'] + G]
    d.to_csv(f'{outdir}/genes/{name[gid]}.tsv.gz', sep='\t', index=False)

# gene summary for the Abundance tab
cat = cls.structural_category.where(cls.structural_category.isin(['FSM', 'ISM', 'NIC', 'NNC']), 'other')
summ = pd.DataFrame({'gene_id': cls.gene_id, 'cat': cat, 'nreads': cls.nreads})
n = summ.groupby(['gene_id', 'cat']).size().unstack(fill_value=0).add_prefix('n')
r = summ.groupby(['gene_id', 'cat']).nreads.sum().unstack(fill_value=0).add_prefix('n').add_suffix('reads')
tot = summ.groupby('gene_id').agg(nreads=('nreads', 'sum'), ntranscripts=('cat', 'size'))
out = mane[['gene_id', 'gene_name', 'transcript_id', 'strand']].merge(tot.join(n).join(r), left_on='gene_id', right_index=True)
out = out[out.gene_name.isin([name[g] for g in ev.gene_id.unique()])]
out.to_csv(f'{outdir}/transcripts.txt', sep='\t', index=False)
print('genes written', len(out), flush=True)
