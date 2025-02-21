# Long-read brain transcript coverage visualisation

This is a shiny app to visualise the coverage of regions expressed across different transcript isoforms in the frontal cortex across development

The y-axis on the plots is similar to a PEXT score, however unlike the PEXT score (which is calucalted using short-read RNAseq from GTeX) we have the full-length transcript isoforms. In our case, the score represents how often a base is included in transcripts found in the frontal cortex

For visualisation of the individual transcript isoforms see [isoforms.com](isoforms.com) and https://www.epigenomicslab.com/brainisoform/ for more details, and view the UCSC track hub at https://genome.ucsc.edu/s/rb520/IsoformTrackHub

To run the app, you can run locally using the `app.R` (no need to download data files) or visit https://chundruvk.shinyapps.io/LRBrainCoverage/

## Citation

Bamford RA, Leung SK, Chundru VK. *et al* An atlas of expressed transcripts in the prenatal and postnatal human cortex. *bioRxiv*. 2024 https://www.biorxiv.org/content/10.1101/2024.05.24.595768
