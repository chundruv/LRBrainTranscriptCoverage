"""Classify isoforms by splice-site motif on the assigned strand and on the opposite strand.

Classes
  canonical            all junctions GT-AG / GC-AG / AT-AC
  strand_flip          no junction canonical, every junction canonical on the opposite strand
                       (CT-AC -> GT-AG, CT-GC -> GC-AG, GT-AT -> AT-AC): strand assigned wrongly
  mixed                some junctions canonical, some not: strand is right, the odd junction is local
  noncanonical         no junction canonical and not rescued by flipping

Published set = isoforms in the CPAT CDS GTF with >= t reads in >= t Whole samples.
Usage: python strand_check.py <classification> <junctions> <gencode.gtf.gz> <stable11.tsv> <cpat.gtf> <t> <out_prefix>
"""
import sys, gzip
import pandas as pd

cls_path, junc_path, gencode_path, st11_path, set_gtf, t, out = sys.argv[1:8]
t = int(t)
CANON = {'GTAG', 'GCAG', 'ATAC'}
comp = str.maketrans('ACGT', 'TGCA')
rc = lambda s: s.translate(comp)[::-1]
flip = lambda s: rc(s[2:4]) + rc(s[0:2])  # motif read from the opposite strand
opp = {'+': '-', '-': '+'}

# junctions (SQANTI3 *_junctions.txt; may or may not carry a header)
cols = ['isoform', 'chrom', 'strand', 'junction_number', 'start', 'end', 'transcript_coord', 'junction_category',
        'start_site_category', 'end_site_category', 'diff_start', 'diff_end', 'bite_junction', 'splice_site', 'canonical']
j = pd.read_csv(junc_path, sep='\t', header=None, names=cols, usecols=range(15), dtype=str)
j = j[j.isoform != 'isoform']
j[['start', 'end']] = j[['start', 'end']].astype(int)
j['splice_site'] = j.splice_site.str.upper()
j['canon'] = j.splice_site.isin(CANON)
j['canon_flip'] = j.splice_site.map(flip).isin(CANON)

# annotated GENCODE introns
ex = pd.read_csv(gencode_path, sep='\t', comment='#', header=None, usecols=[0, 2, 3, 4, 6, 8],
                 names=['chrom', 'type', 'start', 'end', 'strand', 'attr'])
ex = ex[ex.type == 'exon']
ex['tx'] = ex.attr.str.extract(r'transcript_id "([^"]+)"')[0]
ex = ex.sort_values(['tx', 'start'])
intr = pd.DataFrame({'chrom': ex.chrom, 'strand': ex.strand, 'start': ex.end + 1,
                     'end': ex.groupby('tx').start.shift(-1) - 1, 'tx': ex.tx}).dropna()
annot = set(zip(intr.chrom, intr.start, intr.end.astype(int), intr.strand))
j['annot_same'] = [k in annot for k in zip(j.chrom, j.start, j.end, j.strand)]
j['annot_opp'] = [k in annot for k in zip(j.chrom, j.start, j.end, j.strand.map(opp))]

# per isoform
iso = j.groupby('isoform').agg(chrom=('chrom', 'first'), strand=('strand', 'first'), n_junc=('canon', 'size'),
                               n_canon=('canon', 'sum'), n_canon_flip=('canon_flip', 'sum'),
                               n_ctac=('splice_site', lambda s: (s == 'CTAC').sum()),
                               n_annot_same=('annot_same', 'sum'), n_annot_opp=('annot_opp', 'sum'),
                               motifs=('splice_site', lambda s: ','.join(s)))
iso['class'] = 'noncanonical'
iso.loc[iso.n_canon == iso.n_junc, 'class'] = 'canonical'
iso.loc[(iso.n_canon == 0) & (iso.n_canon_flip == iso.n_junc), 'class'] = 'strand_flip'
iso.loc[(iso.n_canon > 0) & (iso.n_canon < iso.n_junc), 'class'] = 'mixed'

# intron chain on the opposite strand already present as another isoform?
chain = j.sort_values(['isoform', 'start']).groupby('isoform').apply(
    lambda d: d.chrom.iloc[0] + ':' + ';'.join(f'{a}-{b}' for a, b in zip(d.start, d.end)))
iso['chain'] = chain
have = set(iso.chain + iso.strand)
iso['opp_chain_exists'] = [(c + opp[s]) in have for c, s in zip(iso.chain, iso.strand)]

# classification fields
c = pd.read_csv(cls_path, sep='\t', usecols=['isoform', 'exons', 'structural_category', 'associated_gene', 'subcategory',
                                             'polyA_motif_found', 'perc_A_downstream_TTS', 'whole_nreads', 'whole_nsamples'])
iso = iso.join(c.set_index('isoform'), how='inner')
iso['polyA_motif_found'] = iso.polyA_motif_found.astype(str).str.upper() == 'TRUE'

# published set
g = pd.read_csv(set_gtf, sep='\t', comment='#', header=None, usecols=[2, 8])
in_gtf = set(g.loc[g[2] == 'transcript', 8].str.extract(r'transcript_id "([^"]+)"')[0])
c['published'] = c.isoform.isin(in_gtf) & (c.whole_nreads >= t) & (c.whole_nsamples >= t)
iso['published'] = iso.index.isin(set(c.isoform[c.published]))

# STable 11 isoforms
st = pd.read_csv(st11_path, sep='\t')
iso['in_stable11'] = iso.index.isin(set(st.Isoform))

# summaries
with open(out + '_summary.txt', 'w') as f:
    p = lambda *a: print(*a, file=f)
    pub = c[c.published]
    p('published set: transcripts', len(pub), 'multi-exon', int((pub.exons > 1).sum()))
    p(pub.structural_category.value_counts().to_string())
    st_iso = iso[iso.in_stable11]
    iso = iso[iso.published]
    j = j[j.isoform.isin(iso.index)]
    p('\nmulti-exon isoforms in published set', len(iso), ' whole reads', int(iso.whole_nreads.sum()))
    s = iso.groupby('class').agg(n=('n_junc', 'size'), reads=('whole_nreads', 'sum'),
                                 two_exon=('exons', lambda e: (e == 2).mean()),
                                 all_junc_annot_opp=('n_annot_opp', lambda x: (x == iso.loc[x.index, 'n_junc']).mean()),
                                 any_junc_annot_opp=('n_annot_opp', lambda x: (x > 0).mean()),
                                 opp_chain_exists=('opp_chain_exists', 'mean'),
                                 polyA_motif=('polyA_motif_found', 'mean'),
                                 median_percA=('perc_A_downstream_TTS', 'median'))
    p(s.round(3).to_string())
    p('\nstructural category by class'); p(pd.crosstab(iso.structural_category, iso['class']).to_string())
    p('\nexons in strand_flip'); p(iso[iso['class'] == 'strand_flip'].exons.clip(upper=6).value_counts().sort_index().to_string())
    nc = j[~j.canon].merge(iso[['class']], left_on='isoform', right_index=True)
    p('\nnon-canonical junctions by motif x class (top 12 motifs)')
    top = nc.splice_site.value_counts().index[:12]
    p(pd.crosstab(nc.splice_site, nc['class']).loc[top].to_string())
    p('\nCT-AC junctions in mixed isoforms: matches annotated intron same strand / opposite strand / novel')
    m = nc[(nc['class'] == 'mixed') & (nc.splice_site == 'CTAC')]
    p(len(m), int(m.annot_same.sum()), int(m.annot_opp.sum()), m.junction_category.value_counts().to_dict())
    p('\nSTable 11 isoforms by class (all, then in published set)')
    p(st_iso['class'].value_counts().to_string()); p(st_iso[st_iso.published]['class'].value_counts().to_string())
    stj = st.assign(cls=st.Isoform.map(iso['class']))
    p('\nSTable 11 junction rows: CTAC by class'); p(pd.crosstab(stj.SpliceSite == 'CTAC', stj.cls).to_string())
    p('\nexample ONT4.7218.10662'); p(st_iso.loc[st_iso.index == 'ONT4.7218.10662'].T.to_string())

iso[iso['class'] != 'canonical'].drop(columns='chain').to_csv(out + '_noncanonical_isoforms.tsv.gz', sep='\t')
st_iso.drop(columns='chain').to_csv(out + '_stable11_isoforms.tsv.gz', sep='\t')
