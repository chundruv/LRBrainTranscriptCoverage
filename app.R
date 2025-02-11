library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(data.table)


# Load the GTF file and the bigWig data (replace with your actual paths)
gtf_gr <- import.gff(con = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz', format = "gtf")
mane<-fread('http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/mane/current/MANE.GRCh38.v1.3.transcripts_by_gene.tsv', stringsAsFactors = F, data.table=F)
gtf_gr<-gtf_gr[(gtf_gr$type=='gene') | (gtf_gr$transcript_id%in%mane$MANE_Select_Ensembl_id)]

# Define UI
ui <- fluidPage(
    titlePanel("Long read brain pext score"),
    theme = shinythemes::shinytheme("united"),
    # Sidebar for inputs
    page_sidebar(
        sidebar = sidebar(width = 250,
            textInput("gene", label = "Gene Name", value = "MYCBP2"),
            checkboxGroupInput(
                "Groups1", 
                label = "Broad developmental groups", 
                choices = c(
                    "Total expression" = 'total_exp',
                    "Prenatal expression" = 'prenatal_exp',
                    "Postnatal expression" = 'postnatal_exp'
                ),
                selected = c('total_exp', 'prenatal_exp', 'postnatal_exp')
            ),
            checkboxGroupInput(
                "Groups2", 
                label = "Fine-scale prenatal developmental groups", 
                choices = c(
                    "1st trimester expression" = 'prenatal1sttrimester_exp',
                    "2nd trimester expression" = 'prenatal2ndtrimester_exp',
                    "3rd trimester expression" = 'prenatal3rdtrimester_exp'
                ),
                selected = c()
            ),
            checkboxGroupInput(
                "Groups3", 
                label = "Fine-scale postnatal developmental groups", 
                choices = c(
                    "Child expression" = 'postnatalchild_exp',
                    "Adult expression" = 'postnataladult_exp',
                    "Older adult expression" = 'postnatalelderly_exp'
                ),
                selected = c()
            ),
            actionButton("submit", "Submit!"),
            radioButtons("rb", "Scale", choiceNames = list("Linear", "Log(1-pext)"), choiceValues = list("linear", "log")) 
        ),
        
        # Main panel for plot
        card(
            full_screen = TRUE,
            plotlyOutput("plot")
        )
    )
)

# Define server logic
server <- function(input, output) {


    # Reactive expression to load and process expression data for selected groups
    dataset <- eventReactive(input$submit, {
        x1 <- fread(paste0('https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/data/',input$gene,'.txt.gz'), stringsAsFactors=F, data.table=F)
        names(x1)<-c('chr', 'pos', 'gene', 'total_exp', 'prenatal_exp', 'postnatal_exp', 'prenatal1sttrimester_exp', 
                     'prenatal2ndtrimester_exp', 'prenatal3rdtrimester_exp', 'postnatalchild_exp', 'postnataladult_exp','postnatalelderly_exp')
        x1$bins <- cut_interval(x1$pos, length=10)
        x1
    }, ignoreNULL = T)
    
    dataset1 <- eventReactive(input$submit, {
        tmp<-dataset()%>%group_by(bins)%>%
            summarise(chr=unique(chr), minpos=min(pos), maxpos=max(pos), 
                      midpos=mean(pos),total_exp=mean(total_exp), prenatal_exp=mean(prenatal_exp), postnatal_exp=mean(postnatal_exp), 
                      prenatal1sttrimester_exp=mean(prenatal1sttrimester_exp), prenatal2ndtrimester_exp=mean(prenatal2ndtrimester_exp), 
                      prenatal3rdtrimester_exp=mean(prenatal3rdtrimester_exp), postnatalchild_exp=mean(postnatalchild_exp), 
                      postnataladult_exp=mean(postnataladult_exp), postnatalelderly_exp=mean(postnatalelderly_exp))
        as.data.frame(tmp)
    }, ignoreNULL = T)

    
    # Reactive expression to process gene annotations
    dataset2 <- eventReactive(input$submit, {
        as.data.frame(gtf_gr[gtf_gr$gene_name == input$gene & gtf_gr$type == 'exon'])
    }, ignoreNULL = T)
    
    dataset3 <- eventReactive(input$submit, {
        gtf_exons <- gtf_gr[gtf_gr$gene_name == input$gene & gtf_gr$type == 'exon', ]
        gtf_transcript <- gtf_gr[gtf_gr$gene_name == input$gene & gtf_gr$type == 'transcript', ]
        gtf_introns <- GenomicRanges::setdiff(gtf_transcript, gtf_exons)
        gtf_introns$transcript_id <- unique(gtf_exons$transcript_id)
        as.data.frame(gtf_introns)
    }, ignoreNULL = T)
    
    geneshapes<-eventReactive(input$submit, {
        shapes <- list()
        
        # Add exon shapes
        for(i in 1:nrow(dataset2())) {
            exon <- dataset2()[i, ]
            shapes[[length(shapes) + 1]] <- list(
                type = "scatter",
                mode = "rect",
                x0 = exon$start, x1 = exon$end,
                y0 = 0.7, y1 = 1.3,
                fillcolor = 'black',
                line = list(color = "black"),
                text = paste0('Transcript ID: ', exon$transcript_id, '\nExon ID: ', 
                              exon$exon_id, '\nExon range: ', exon$seqnames, ':', 
                              exon$start, '-', exon$end),
                hoverinfo='text'
            )
        }
        #Intron lines
        for(i in 1:nrow(dataset3())) {
            intron <- dataset3()[i, ]
            shapes[[length(shapes) + 1]] <- list(
                type = "line",
                x0 = intron$start, x1 = intron$end,
                y0 = 1, y1 = 1,
                line = list(color = "grey", width = 3),
                text = paste0('Transcript ID: ', intron$transcript_id, '\nIntron range: ', 
                              intron$seqnames, ':',intron$start, '-', intron$end),
                hoverinfo='text'
            )
        }
        shapes
    }, ignoreNULL = T)
    
    # Render the plot
    
    plot<-eventReactive(input$submit,{
        p<-list()

        for(i in unique(c(input$Groups1,input$Groups2,input$Groups3))){
            j <- switch(i,
                   'total_exp' = "#999999",
                   'prenatal_exp' = "#E69F00",
                   'postnatal_exp' = "#56B4E9",
                   'prenatal1sttrimester_exp' = "#009E73",
                   'prenatal2ndtrimester_exp' = "#CC79A7",
                   'prenatal3rdtrimester_exp' = "#0072B2",
                   'postnatalchild_exp' = "#D55E00",
                   'postnataladult_exp' = "#41ab5d",
                   'postnatalelderly_exp' = "#000000")

                p[[i]] <- plotly_build(plot_ly() %>% add_trace(dataset1(), x = dataset1()$midpos, y = dataset1()[,i], 
                                       color=~I(j), width = 10, type = 'bar', name=i, hoverinfo = "text", hovertext= 
                                       ~paste('<br><b>Position range</b>: ',dataset1()$chr,':',dataset1()$minpos,'-',dataset1()$maxpos,'<br>',
                                             
                                             '<br><b>Group</b>: ',i,'<br>',
                                             
                                             '<br><b>Mean PEXT</b>: ',format(round(dataset1()[,i], 4), nsmall = 4),'<br>'))%>%
                                           layout(hovermode = 'x unified', xaxis = list(
                                           title = "",
                                           showgrid = FALSE,
                                           showticklabels = FALSE  # Hide axis labels and ticks
                                       ),
                                       yaxis = list(type = input$rb,
                                           title = "", range<-c(0,1)
                                       )))
        }

        p1<-subplot(p, nrows = length(c(input$Groups1,input$Groups2,input$Groups3)), shareX = TRUE, shareY = TRUE)%>%layout(hovermode = 'x unified', xaxis = list(
            title = "", 
            showgrid = FALSE,
            showticklabels = FALSE
        ),
        yaxis = list(
            title = ""
        ), legend = list(orientation = 'h'))
        
        
        p2<-plot_ly()
        for (i in 1:length(geneshapes())) {
            shape <- geneshapes()[[i]]
            p2 <- p2 %>% add_trace(
                        type = "scatter",
                        mode = "lines",
                        x = c(shape$x0, shape$x1),
                        y = c(shape$y0, shape$y1),
                        showlegend=F,
                        line = list(color = shape$line$color, width = shape$line$width),
                        text = shape$text,  # Add hover text here
                        hoverinfo = "text"   # Enable hoverinfo to display the 'text'
                    )
        }
        p2<-p2%>%layout(shapes = geneshapes(), xaxis=list(showgrid = FALSE,showticklabels = TRUE,title=""), 
                        yaxis=list(title = "", showgrid = FALSE, showticklabels = FALSE))

        subplot(p1,p2, nrows=2, heights = c(0.9,0.1), shareX = T)%>%
            layout(title = list(text = paste0("Gene: ", input$gene,"\nMANEselect: ", unique(dataset3()$transcript_id))), margin=list(t=50))

    }, ignoreNULL = T)
    
    output$plot<-renderPlotly({plot()})
    
}

# Run the application 
shinyApp(ui = ui, server = server)

