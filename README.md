# Long-read brain transcript coverage visualisation

This is a shiny app to visualise the coverage of regions expressed across different transcript isoforms in the frontal cortex across development

The y-axis on the plots is similar to a PEXT score, however unlike the PEXT score (which is calucalted using short-read RNAseq from GTeX) we have the full-length transcript isoforms. In our case, the score represents how often a base is included in transcripts found in the frontal cortex

### Data (`data_sqanti/`)

The app uses two SQANTI3 isoform sets (relaxed filter, mono-exonic isoforms removed), with ORFs predicted by CPAT:

- `10r10s/` – isoforms with ≥10 reads in ≥10 samples (`WholeTargeted_relaxedsqanti_10reads10samples_CPAT_CDS.gtf`)
- `2r2s/` – isoforms with ≥2 reads in ≥2 samples (`WholeTargeted_relaxedsqanti_2reads2samples_CPAT_CDS.gtf`)

Mis-stranded isoforms are excluded: those with no canonical junction on the assigned strand but all junctions canonical on the opposite strand, plus isoforms with a CT-AC junction in an otherwise canonical isoform. That is 887 isoforms in 10r10s and 7,040 in 2r2s; the lists are at `13_genetics/strand_filtered/removed_isoforms_*.tsv` on isca. Per-sample read counts come from the matching `sqantifiltered_monoexonicfiltered_*_classification.txt`. Only the 47 whole-transcriptome samples are used. Each `genes/<GENE>.tsv.gz` stores intervals with constant coverage. For each developmental group, coverage at a base and region (CDS, 5'UTR, 3'UTR, non-coding) is the summed CPM of isoforms that contain the base with that label, divided by the summed CPM of all isoforms of the gene. Labels come from each isoform's own CPAT ORF, so one exon can be CDS in some isoforms and UTR in others; see the "How to read this plot" tab. `mane_track.tsv.gz` holds the MANE Select v1.4 gene models. Files are built with `build/build_coverage.py`.

For visualisation of the individual transcript isoforms see [isoforms.com](isoforms.com) and https://www.epigenomicslab.com/brainisoform/ for more details, and view the UCSC tracks at [https://genome.ucsc.edu/s/rb520/relaxed](https://genome.ucsc.edu/s/rb520/relaxed)

To run the app, you can run locally using the `app.R` (no need to download data files) or visit https://chundruvk.shinyapps.io/LRBrainCoverage/

## Citation

Bamford RA, Leung SK, Chundru VK. *et al* An atlas of expressed transcripts in the prenatal and postnatal human cortex. *bioRxiv*. 2024 https://www.biorxiv.org/content/10.1101/2024.05.24.595768
