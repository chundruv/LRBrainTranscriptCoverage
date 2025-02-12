library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(ggplot2)
library(R.utils)
library(data.table)
library(readr)

transcripts<-fread('https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/transcripts.txt', stringsAsFactors = F, data.table = F)

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
            actionButton("submit", "Submit!", icon("paper-plane"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
            radioButtons("rb", "Scale", choiceNames = list("Linear", "Log"), choiceValues = list("linear", "log")) 
        ),
        
        # Main panel for plot
        card(
            full_screen = TRUE,
            plotlyOutput("plot")
        ),
        card(
            htmlOutput("text")
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
    
    geneshapes<-eventReactive(input$submit, {
        read_rds(paste0('https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/data/',input$gene,'.RDS'))
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
            k <- switch(i,
                        'total_exp' = "Total expression",
                        'prenatal_exp' = "Prenatal expression",
                        'postnatal_exp' = "Postnatal expression",
                        'prenatal1sttrimester_exp' = "1st trimester expression",
                        'prenatal2ndtrimester_exp' = "2nd trimester expression",
                        'prenatal3rdtrimester_exp' = "3rd trimester expression",
                        'postnatalchild_exp' = "Child expression",
                        'postnataladult_exp' = "Adult expression",
                        'postnatalelderly_exp' = "Older adult expression")

            if(input$rb=='linear'){
                p[[i]] <- plotly_build(plot_ly() %>% add_trace(dataset1(), x = dataset1()$midpos, y = dataset1()[,i], 
                                       color=~I(j), width = 10, type = 'bar', name=k, hoverinfo = "text", hovertext= 
                                       ~paste('<br><b>Position range</b>: ',dataset1()$chr,':',dataset1()$minpos,'-',dataset1()$maxpos,'<br>',
                                             
                                             '<br><b>Group</b>: ',i,'<br>',
                                             
                                             '<br><b>Mean PEXT</b>: ',format(round(dataset1()[,i], 4), nsmall = 4),'<br>'))%>%
                                           layout(hovermode = 'x unified', xaxis = list(
                                           title = "",
                                           showgrid = FALSE,
                                           showticklabels = FALSE  # Hide axis labels and ticks
                                       ),
                                       yaxis = list(range=c(0,1),
                                           title = "", tickmode="array", tickvals=c(0,0.25,0.5,0.75,1)
                                       )))
            }else{
                p[[i]] <- plotly_build(plot_ly() %>% add_trace(dataset1(), x = dataset1()$midpos, y = dataset1()[,i], 
                                                               color=~I(j), width = 10, type = 'bar', name=k, hoverinfo = "text", hovertext= 
                                                                   ~paste('<br><b>Position range</b>: ',dataset1()$chr,':',dataset1()$minpos,'-',dataset1()$maxpos,'<br>',
                                                                          
                                                                          '<br><b>Group</b>: ',i,'<br>',
                                                                          
                                                                          '<br><b>Mean PEXT</b>: ',format(round(dataset1()[,i], 4), nsmall = 4),'<br>'))%>%
                                           layout(hovermode = 'x unified', xaxis = list(
                                               title = "",
                                               showgrid = FALSE,
                                               showticklabels = FALSE  # Hide axis labels and ticks
                                           ),
                                           yaxis = list(type = 'log', range=c(-5,0), nticks=4,
                                                        title = ""
                                           )))
            }
        }

        p1<-subplot(p, nrows = length(c(input$Groups1,input$Groups2,input$Groups3)), shareX = TRUE, shareY = TRUE)%>%
            layout(hovermode = 'x unified', xaxis = list(title = "", 
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

        if(transcripts[which(transcripts$gene_name==input$gene),3]=='+'){
            strand<-'+ (\u2192)'
        }else if(transcripts[which(transcripts$gene_name==input$gene),3]=='-'){
            strand<-'- (\u2190)'
        }else{
            strand<-''
        }
        subplot(p1,p2, nrows=2, heights = c(0.9,0.1), shareX = T)%>%
            layout(title = list(text = paste0("Gene: ", input$gene, ';\tStrand: ', strand, "\nMANEselect: ", 
                                              transcripts[which(transcripts$gene_name==input$gene),2])), margin=list(t=50))

    }, ignoreNULL = T)
    
    output$plot<-renderPlotly({plot()})
    text<-eventReactive(input$submit,{
        gn=transcripts[which(transcripts$gene_name==input$gene),]
        HTML("Number of reads = ", gn$nreads, 
               "<br>Number of unique transcript isoforms = ", gn$ntranscripts, 
               "<br>Number of FSM transcripts = ", gn$nFSM, " (", gn$nFSMreads, ' reads)', 
               "<br>Number of ISM transcripts = ", gn$nISM, " (", gn$nISMreads, ' reads)',
               "<br>Number of NIC transcripts = ", gn$nNIC, " (", gn$nNICreads, ' reads)',
               "<br>Number of NNC transcripts = ", gn$nNNC, " (", gn$nNNCreads, ' reads)',)}, ignoreNULL = T)
    output$text<-renderText({text()})
}

# Run the application 
shinyApp(ui = ui, server = server)

