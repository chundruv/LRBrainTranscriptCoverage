library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)
library(R.utils)
library(data.table)
library(shinydisconnect)
library(shinycssloaders)

# Data are read from the local checkout when present, otherwise from GitHub
base <- if (dir.exists('data_sqanti')) 'data_sqanti' else 'https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/data_sqanti'

datasets <- c("SQANTI relaxed filter, \u226510 reads in \u226510 samples" = '10r10s',
              "SQANTI relaxed filter, \u22652 reads in \u22652 samples" = '2r2s')
transcripts <- lapply(datasets, function(d) fread(paste0(base, '/', d, '/transcripts.txt'), data.table = F))
names(transcripts) <- datasets
mane_track <- fread(paste0(base, '/mane_track.tsv.gz'), data.table = F)

region_cols <- c("CDS" = "#0072B2", "5'UTR" = "#E69F00", "3'UTR" = "#009E73", "non-coding" = "#999999")
group_names <- c('total_exp' = "Total expression",
                 'prenatal_exp' = "Prenatal expression",
                 'postnatal_exp' = "Postnatal expression",
                 'prenatal1sttrimester_exp' = "1st trimester expression",
                 'prenatal2ndtrimester_exp' = "2nd trimester expression",
                 'prenatal3rdtrimester_exp' = "3rd trimester expression",
                 'postnatalchild_exp' = "Child expression",
                 'postnataladult_exp' = "Adult expression",
                 'postnatalelderly_exp' = "Older adult expression")

jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each($els, function(idx, el) {
    var $el = $(el);
    var $proxy = $("#" + $el.data("proxyClick"));
    $el.keydown(function (e) { if (e.keyCode == 13) { $proxy.click(); } });
  });
});
'

plot_note <- HTML("<p style='font-size:0.9em; margin:6px 10px;'>
  Bar height = fraction of the gene's long-read expression from isoforms that include the base.
  Colours split that fraction by how <b>each isoform's own predicted ORF</b> (CPAT) uses the base, so one exon can be
  CDS in some isoforms and UTR in others. The bottom track is the GENCODE/MANE Select annotation for comparison
  (thick = CDS, thin = UTR, green = start codon, vermillion = stop codon). See <b>How to read this plot</b>.</p>")

help_text <- HTML("
<h4>What the bars show</h4>
<p>For every base of the gene, the y-axis is the fraction of the gene's long-read expression that comes from transcript
isoforms containing that base. Expression is the sum, across the samples in the selected developmental group, of each
isoform's read count normalised to counts per million (CPM) within its sample. A value of 1 means every detected isoform of
the gene includes the base; 0.2 means isoforms carrying 20% of the gene's expression include it. This is the same idea as the
pext score, but measured directly on full-length isoforms rather than inferred from short reads.</p>

<h4>What the colours mean</h4>
<p>The coding sequence is <b>not</b> taken from GENCODE or MANE. Every long-read isoform, including those that match an
annotated transcript (FSM), was given its own open reading frame by CPAT (best ORF; coding probability &gt; 0.364).
At each base, the bar is split according to what that base is <i>in each isoform that contains it</i>:</p>
<ul>
<li><span style='color:#0072B2'><b>CDS</b></span> &ndash; inside the isoform's predicted ORF</li>
<li><span style='color:#E69F00'><b>5'UTR</b></span> &ndash; upstream of the predicted start codon</li>
<li><span style='color:#009E73'><b>3'UTR</b></span> &ndash; downstream of the predicted stop codon</li>
<li><span style='color:#999999'><b>non-coding</b></span> &ndash; the isoform has no ORF above the CPAT coding threshold</li>
</ul>
<p>The coloured segments are stacked, so their total height equals the overall inclusion of the base.</p>

<h4>Why the same exon can be both CDS and UTR</h4>
<p>A stacked bar with both blue and orange (or green) means the base is coding in some isoforms and untranslated in others.
The bar shows the share of the gene's expression in each state. This is expected whenever isoforms of one gene use
different ORFs. Common causes are:</p>
<ul>
<li><b>Alternative first exon or transcription start site.</b> If an isoform lacks the exon carrying the annotated start
codon, its best ORF begins at a downstream ATG. Annotated coding exons upstream of that ATG then become 5'UTR in that isoform.</li>
<li><b>Alternative splicing that shifts the frame or creates an earlier stop codon.</b> Annotated coding exons downstream
of the new stop become 3'UTR in that isoform.</li>
<li><b>A different ORF on the same exon structure.</b> CPAT reports one best ORF per isoform. Occasionally this differs from
the annotated ORF even when the exon structure matches.</li>
<li><b>Non-coding isoforms</b> (for example retained-intron transcripts) contribute grey at bases they share with
coding isoforms.</li>
</ul>
<p>Where the colours disagree with the MANE track at the bottom, the isoforms expressed in cortex predict a protein that
starts or ends differently from the annotated one. Where they agree, the brain isoforms use the annotated ORF.</p>

<h4>Example: ASGR1</h4>
<p>ASGR1 is expressed at low levels in cortex (174 reads across 5 isoforms in the &ge;10 reads/&ge;10 samples set).
None of these isoforms contains the MANE exon carrying the annotated start codon (chr17:7178561-7178563; ASGR1 is on the
minus strand). They start instead from alternative first exons at roughly chr17:7178085-7178188. Without the annotated
ATG, CPAT's best ORF in every isoform begins at an internal ATG in the fourth MANE coding exon (chr17:7176881).
The first MANE coding exon is absent from these isoforms (no bar). The second and third coding exons, and the start of
the fourth, lie upstream of the predicted ATG and are therefore 5'UTR in all cortex isoforms. This shows as orange over
exons that the MANE track draws as thick (CDS). Downstream of chr17:7176881 the isoforms use the annotated reading frame, and the
bars are blue. In the permissive &ge;2 reads/&ge;2 samples set, only ~1% of ASGR1 expression includes the annotated start
codon, which is why a small blue segment appears there. In short, the orange bars over annotated coding exons are a
property of the isoforms expressed in brain, not a plotting error. They indicate that cortex ASGR1 transcripts are
predicted to encode an N-terminally truncated protein.</p>

<h4>Gene model track</h4>
<p>The bottom track shows the MANE Select transcript (GENCODE annotation): thick black = CDS, thin black = UTR,
grey = intron, green = start codon, vermillion = stop codon. Hover over a feature for its coordinates.</p>

<h4>Datasets</h4>
<p>Both datasets come from the same SQANTI3 run (relaxed filter, mono-exonic isoforms removed), with CDS predicted by CPAT.
They differ in how much read support an isoform needs:</p>
<ul>
<li><b>&ge;10 reads in &ge;10 samples</b> &ndash; high-confidence isoforms, recommended default.</li>
<li><b>&ge;2 reads in &ge;2 samples</b> &ndash; permissive, includes rarer isoforms. Expect more mixed CDS/UTR calls.</li>
</ul>
<p>Isoforms assigned to the wrong strand were removed from both sets. These are isoforms in which no junction is canonical
on the assigned strand but every junction is canonical (GT-AG, GC-AG, AT-AC) on the opposite strand. Isoforms carrying a
CT-AC junction within an otherwise canonical isoform were also removed.</p>
<p>Only the 47 whole-transcriptome cortex samples contribute to the coverage values. The Targeted panel samples are not used.</p>
")

ui <- fluidPage(
    tags$head(tags$script(HTML(jscode))),
    tags$head(tags$style(HTML(".control-label{font-weight:bold;}"))),
    disconnectMessage(text = "App disconnected. Please see https://github.com/chundruv/LRBrainTranscriptCoverage. Run the app.R locally, no need to download data files",
                      refresh = "", background = "#646464e6", size = 35, width = "full", top = "center",
                      colour = "white", overlayColour = "#999", overlayOpacity = 0.4),
    titlePanel(windowTitle = "Long read brain coverage",
        fluidRow(
            column(9, "Long read brain coverage"),
            column(3, list(tags$a(img(src = "https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png", width = "50", height = "50"),
                                  href = "https://github.com/chundruv/LRBrainTranscriptCoverage", target = '_blank'),
                           tags$a(img(src = "https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/logos/paradigm.png", width = "150", height = "50"),
                                  href = "https://paradigmgenomics.org", target = '_blank'),
                           tags$a(img(src = "https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/logos/hdr-cpg.png", width = "100", height = "50"),
                                  href = "https://www.epigenomicslab.com/brainisoform/", target = '_blank')))
        )
    ),
    `data-proxy-click` = "submit",
    page_sidebar(
        sidebar = sidebar(width = 300,
            radioButtons("dataset", "Isoform set", choices = datasets, selected = '10r10s'),
            selectizeInput('gene', label = 'Please type gene name', choices = NULL),
            numericInput("width", label = "Binwidth", value = 10),
            checkboxGroupInput("regions", label = "Show coverage by predicted region",
                               choices = names(region_cols), selected = names(region_cols)),
            checkboxGroupInput("Groups1", label = "Broad developmental groups",
                               choices = setNames(names(group_names)[1:3], group_names[1:3]),
                               selected = names(group_names)[1:3]),
            checkboxGroupInput("Groups2", label = "Fine-scale prenatal developmental groups",
                               choices = setNames(names(group_names)[4:6], group_names[4:6])),
            checkboxGroupInput("Groups3", label = "Fine-scale postnatal developmental groups",
                               choices = setNames(names(group_names)[7:9], group_names[7:9])),
            radioButtons("rb", "Scale", choiceNames = list("Linear", "Log"), choiceValues = list("linear", "log")),
            actionButton("submit", "Submit!", style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")
        ),
        navset_card_tab(
            full_screen = TRUE,
            nav_panel("Coverage",
                      plot_note,
                      plotlyOutput("plot", height = "100%") %>% withSpinner(color = "#0dc5c1") %>% bslib::as_fill_carrier()),
            nav_panel("How to read this plot", help_text),
            nav_panel("Abundance", htmlOutput("text")),
            nav_panel("About",
                      HTML("Samples (whole-transcriptome frontal cortex, n = 47):<br>
                       Prenatal 1st trimester (&lt;14 weeks post-conception) = 14<br>
                       Prenatal 2nd trimester (14-26 weeks) = 10<br>
                       Prenatal 3rd trimester (27-28 weeks) = 2<br>
                       Postnatal child (0-7 years) = 6<br>
                       Postnatal adult (23-57 years) = 10<br>
                       Postnatal older adult (62-83 years) = 5<br><br>
                       <p>Isoforms: SQANTI3 relaxed filter with mono-exonic isoforms removed, at two read-support thresholds
                       (&ge;10 reads in &ge;10 samples, or &ge;2 reads in &ge;2 samples). ORFs predicted with CPAT.</p>
                       <p>See <a href='https://www.epigenomicslab.com/brainisoform/' target='_blank'>https://www.epigenomicslab.com/brainisoform/</a>
                       for further details and visualisations of this data.</p>
                       <p>For details of the samples and methods, and citation information, see
                       <a href='https://www.biorxiv.org/content/10.1101/2024.05.24.595768' target='_blank'>our bioRxiv preprint</a>.</p>"))
        )
    )
)

# one coverage panel per developmental group
coverage_panel <- function(x, grp, regions, log_scale, width, showleg) {
    p <- plot_ly()
    for (r in regions) {
        d <- x[x$region == r, ]
        p <- p %>% add_trace(x = d$midpos, y = d[[grp]], type = 'bar', width = width,
                             marker = list(color = region_cols[[r]]), name = r, legendgroup = r, showlegend = showleg,
                             hoverinfo = "text",
                             hovertext = paste0('<b>', r, '</b>: ', format(round(d[[grp]], 4), nsmall = 4),
                                                '<br>', d$chr, ':', d$minpos, '-', d$maxpos))
    }
    yaxis <- if (log_scale) list(type = 'log', range = c(-5.1, 0.1), nticks = 6, exponentformat = "e", title = "")
             else list(range = c(-0.1, 1.1), title = "", tickvals = c(0, 0.25, 0.5, 0.75, 1))
    p %>% layout(barmode = 'stack', xaxis = list(showgrid = FALSE, showticklabels = FALSE), yaxis = yaxis) %>%
        add_annotations(text = group_names[[grp]], x = 0.01, y = 1, xref = "paper", yref = "paper",
                        xanchor = "left", yanchor = "bottom", showarrow = FALSE, font = list(size = 14))
}

# MANE Select gene model: CDS thick, UTR thin, introns grey, start/stop codons as markers
gene_model <- function(gene, chr) {
    tr <- mane_track[mane_track$gene_name == gene, ]
    p <- plot_ly()
    styles <- list(intron = list(color = "grey", width = 2), UTR = list(color = "black", width = 6),
                   CDS = list(color = "black", width = 14))
    for (ty in names(styles)) {
        d <- tr[tr$type == ty, ]
        if (nrow(d) == 0) next
        p <- p %>% add_trace(type = "scatter", mode = "lines", x = c(rbind(d$start, d$end, NA)), y = rep(c(1, 1, NA), nrow(d)),
                             line = styles[[ty]], showlegend = FALSE, hoverinfo = "text",
                             text = c(rbind(paste0('MANE ', ty, ': ', chr, ':', d$start, '-', d$end),
                                            paste0('MANE ', ty, ': ', chr, ':', d$start, '-', d$end), NA)))
    }
    codons <- list(start_codon = c("#009E73", "Start codon"), stop_codon = c("#D55E00", "Stop codon"))
    for (ty in names(codons)) {
        d <- tr[tr$type == ty, ]
        if (nrow(d) == 0) next
        p <- p %>% add_trace(type = "scatter", mode = "markers", x = (d$start + d$end) / 2, y = rep(1, nrow(d)),
                             marker = list(color = codons[[ty]][1], size = 12, symbol = "line-ns", line = list(width = 4, color = codons[[ty]][1])),
                             showlegend = FALSE, hoverinfo = "text",
                             text = paste0(codons[[ty]][2], ': ', chr, ':', d$start, '-', d$end))
    }
    p %>% layout(xaxis = list(showgrid = FALSE, showticklabels = TRUE, title = ""),
                 yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, range = c(0.5, 1.5)))
}

server <- function(input, output, session) {

    observeEvent(input$dataset, {
        updateSelectizeInput(session, 'gene', choices = transcripts[[input$dataset]]$gene_name,
                             selected = isolate(input$gene), server = TRUE)
    })

    # intervals -> per-base rows -> bins
    dataset <- eventReactive(input$submit, {
        x <- fread(paste0(base, '/', input$dataset, '/genes/', input$gene, '.tsv.gz'), data.table = F)
        idx <- rep(seq_len(nrow(x)), x$end - x$start + 1)
        y <- x[idx, ]
        y$pos <- y$start + sequence(x$end - x$start + 1) - 1
        y <- y %>% select(-start, -end) %>%
            complete(pos, region = names(region_cols), fill = as.list(setNames(rep(0, length(group_names)), names(group_names)))) %>%
            mutate(chr = x$chr[1])
        y$bins <- cut_interval(y$pos, length = input$width)
        y %>% group_by(bins, region) %>%
            summarise(chr = first(chr), minpos = min(pos), maxpos = max(pos), midpos = mean(pos),
                      across(all_of(names(group_names)), mean), .groups = 'drop') %>%
            as.data.frame()
    }, ignoreNULL = TRUE)

    plot <- eventReactive(input$submit, {
        x <- dataset()
        grps <- unique(c(input$Groups1, input$Groups2, input$Groups3))
        regions <- names(region_cols)[names(region_cols) %in% input$regions]
        p <- lapply(seq_along(grps), function(i)
            coverage_panel(x, grps[i], regions, input$rb == 'log', input$width, showleg = i == 1))
        p1 <- subplot(p, nrows = length(grps), shareX = TRUE, shareY = TRUE)

        gn <- transcripts[[input$dataset]] %>% filter(gene_name == input$gene)
        strand <- c('+' = '+ (\u2192)', '-' = '- (\u2190)')[gn$strand]
        subplot(p1, gene_model(input$gene, x$chr[1]), nrows = 2, heights = c(0.9, 0.1), shareX = TRUE) %>%
            layout(hovermode = 'x unified', legend = list(orientation = 'h'), margin = list(t = 50),
                   title = list(text = paste0("Gene: ", input$gene, ';\tStrand: ', strand, "\nMANE Select: ", gn$transcript_id)))
    }, ignoreNULL = TRUE)

    output$plot <- renderPlotly({plot()})

    text <- eventReactive(input$submit, {
        gn <- transcripts[[input$dataset]] %>% filter(gene_name == input$gene)
        HTML("Isoform set = ", names(datasets)[datasets == input$dataset],
             "<br>Number of reads (whole-transcriptome samples) = ", gn$nreads,
             "<br>Number of unique transcript isoforms = ", gn$ntranscripts,
             "<br>Number of full-splice match transcripts = ", gn$nFSM, " (", gn$nFSMreads, ' reads)',
             "<br>Number of incomplete-splice match transcripts = ", gn$nISM, " (", gn$nISMreads, ' reads)',
             "<br>Number of novel in catalog transcripts = ", gn$nNIC, " (", gn$nNICreads, ' reads)',
             "<br>Number of novel not in catalog transcripts = ", gn$nNNC, " (", gn$nNNCreads, ' reads)',
             "<br>Number of other transcripts (antisense, intergenic, fusion, genic) = ", gn$nother, " (", gn$notherreads, ' reads)',
             "<br><p> See <a href='https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories' target='_blank'>here</a> for details on the classification groups </p>")
    }, ignoreNULL = TRUE)
    output$text <- renderText({text()})
}

shinyApp(ui = ui, server = server)
