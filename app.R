library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(ggplot2)
library(R.utils)
library(data.table)
library(readr)
library(shinydisconnect)
library(shinycssloaders)

transcripts<-fread('https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/transcripts.txt', stringsAsFactors = F, data.table = F)
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

# Define UI
ui <- fluidPage(
    tags$head(tags$script(HTML(jscode))),
    tags$head(
        tags$style(HTML(
            ".control-label{
                               font-weight:bold;
                               }
             "))),
    disconnectMessage(text = "App disconnected. Please see https://github.com/chundruv/LRBrainTranscriptCoverage. Run the app.R locally, no need to download data files",
                      refresh = "",
                      background = "#646464e6",
                      size = 35,
                      width = "full",
                      top = "center",
                      colour = "white",
                      overlayColour = "#999",
                      overlayOpacity = 0.4),
    titlePanel(windowTitle = "Long read brain coverage",
        fluidRow(
            column(9, "Long read brain coverage"), 
            column(3, list(tags$a(img(src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png", width="50", height="50"), 
                                  href="https://github.com/chundruv/LRBrainTranscriptCoverage", target='_blank'),
                           tags$a(img(src="https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/logos/paradigm.png", width="150", 
                                      height="50"), href="https://paradigmgenomics.org", target='_blank'),
                           tags$a(img(src="https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/logos/hdr-cpg.png",width="100",
                                      height="50"), href="https://www.epigenomicslab.com/brainisoform/", target='_blank')))
        )
    ),
    #theme = bslib::bs_theme(bootswatch = "lux", font_scale = 0.9),
    `data-proxy-click` = "submit",
    # Sidebar for inputs
    page_sidebar(
        sidebar = sidebar(width = 300,
                          selectizeInput(
                              inputId = 'gene',
                              label = 'Please type gene name',
                              choices = NULL
                          ),
                          numericInput("width", label = "Binwidth", value=10),
                          checkboxGroupInput(
                              "regions", 
                              label = "Show covereage by predicted consequence", 
                              choices = c(
                                  "Show CDS" = 'cds',
                                  "Show UTR" = 'utr',
                                  "Show non-coding" = 'nc'
                              ),
                              selected = c('cds', 'utr', 'nc')
                          ),
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
                          radioButtons("rb", "Scale", choiceNames = list("Linear", "Log"), choiceValues = list("linear", "log")),
                          actionButton("submit", "Submit!",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          
        ),
        # Main panel for plot
        navset_card_tab(
            full_screen = TRUE,
            nav_panel( "Coverage",
            
            plotlyOutput("plot", height = "100%") %>% withSpinner(color="#0dc5c1") %>% bslib::as_fill_carrier()
        ),
        nav_panel( "Abundance",
            htmlOutput("text")
        ),
        nav_panel("About",
                  HTML("Sample sizes:<br>
                        Prenatal 1st trimester (6-14weeks) = 14<br>
                       Prenatal 2nd trimester (14-27weeks) = 10<br>
                       Prenatal 3rd trimester (27-birth) = 3<br>
                       Postnatal child (0-18yrs) = 5<br>
                       Postnatal adult (18-60yrs) = 10<br>
                       Postnatal elderly (60-83yrs) = 5<br><br>
                        <p> See <a href='https://www.epigenomicslab.com/brainisoform/' target='_blank'>https://www.epigenomicslab.com/brainisoform/</a> and </p><p><a href='isoforms.com' target='_blank'>isoforms.com</a> for further details and visulaisations of this data</p>
                       <p>For further details of the samples and methods, and citation information, see <a href='https://www.biorxiv.org/content/10.1101/2024.05.24.595768' target='_blank'>our biorxiv preprint</a></p>")
        )
        )
))

# Define server logic
server <- function(input, output, session) {
    
    updateSelectizeInput(session, 'gene', choices = transcripts$gene_name, server = TRUE)
    # Reactive expression to load and process expression data for selected groups
    dataset <- eventReactive(input$submit, {
        x1 <- fread(paste0('https://github.com/chundruv/LRBrainTranscriptCoverage/raw/refs/heads/main/data/',input$gene,'.txt.gz'), stringsAsFactors=F, data.table=F)
        x1$bins <- cut_interval(x1$pos, length=input$width)
        x1
    }, ignoreNULL = T)
    
    dataset1 <- eventReactive(input$submit, {
        tmp<-dataset()%>%group_by(bins, region)%>%
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
        x<-dataset1()
        annotations<-list()
        p<-list()
        count<-0
        for(i in unique(c(input$Groups1,input$Groups2,input$Groups3))){
            if(count==0){
                showleg<-T
            }else{
                showleg<-F
            }
            count<-1
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
                p[[i]] <- plot_ly() 
                if('cds'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='CDS'),], x = x[which(x$region=='CDS'),]$midpos, y = x[which(x$region=='CDS'),i], 
                                                  color=~I("#56B4E9"),width = 10, type = 'bar', name="CDS", showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='CDS'),]$chr,':',x[which(x$region=='CDS'),]$minpos,'-',x[which(x$region=='CDS'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: CDS<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='CDS'),i], 4), nsmall = 4),'<br>'))
                }
                if('utr'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='UTR'),], x = x[which(x$region=='UTR'),]$midpos, y = x[which(x$region=='UTR'),i], 
                                                  color=~I("#E69F00"), width = 10, type = 'bar', name="UTR",showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='UTR'),]$chr,':',x[which(x$region=='UTR'),]$minpos,'-',x[which(x$region=='UTR'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: UTR<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='UTR'),i], 4), nsmall = 4),'<br>'))
                }
                if('nc'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='non-coding'),], x = x[which(x$region=='non-coding'),]$midpos, y = x[which(x$region=='non-coding'),i], 
                                                  color=~I("#999999"), width = 10, type = 'bar', name="non-coding",showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='non-coding'),]$chr,':',x[which(x$region=='non-coding'),]$minpos,'-',x[which(x$region=='non-coding'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: non-coding<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='non-coding'),i], 4), nsmall = 4),'<br>'))
                }          
                p[[i]] <-plotly_build(p[[i]]%>%
                                          layout(hovermode = 'x unified', barmode = 'stack', xaxis = list(
                                              title = "",
                                              showgrid = FALSE,
                                              showticklabels = FALSE  # Hide axis labels and ticks
                                          ),
                                          yaxis = list(range=c(-0.1,1.1),
                                                       title = "", tickmode="array", tickvals=c(0,0.25,0.5,0.75,1)
                                          ))%>%
                                          add_annotations(
                                              text = ~I(k),
                                              x = 0.1,
                                              y = 1+((length(unique(c(input$Groups1,input$Groups2,input$Groups3)))*0.8-2)/10),
                                              yref = "paper",
                                              xref = "paper",
                                              xanchor = "center",
                                              yanchor = "top",
                                              showarrow = FALSE,
                                              font = list(size = 15)
                                          ))
            }else{
                p[[i]] <- plot_ly() 
                if('cds'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='CDS'),], x = x[which(x$region=='CDS'),]$midpos, y = x[which(x$region=='CDS'),i], 
                                                  color=~I("#56B4E9"), width = input$width, type = 'bar', name="CDS", showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='CDS'),]$chr,':',x[which(x$region=='CDS'),]$minpos,'-',x[which(x$region=='CDS'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: CDS<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='CDS'),i], 4), nsmall = 4),'<br>'))
                }
                if('utr'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='UTR'),], x = x[which(x$region=='UTR'),]$midpos, y = x[which(x$region=='UTR'),i], 
                                                  color=~I("#E69F00"), width = input$width, type = 'bar', name="UTR",showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='UTR'),]$chr,':',x[which(x$region=='UTR'),]$minpos,'-',x[which(x$region=='UTR'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: UTR<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='UTR'),i], 4), nsmall = 4),'<br>'))
                }
                if('nc'%in%input$regions){
                    p[[i]] <- p[[i]]%>% add_trace(x[which(x$region=='non-coding'),], x = x[which(x$region=='non-coding'),]$midpos, y = x[which(x$region=='non-coding'),i], 
                                                  color=~I("#999999"), width = input$width, type = 'bar', name="non-coding",showlegend=showleg, hoverinfo = "text", hovertext= 
                                                      ~paste('<br><b>Position range</b>: ',x[which(x$region=='non-coding'),]$chr,':',x[which(x$region=='non-coding'),]$minpos,'-',x[which(x$region=='non-coding'),]$maxpos,'<br>',
                                                             '<br><b>Group</b>: ',i,'<br>',
                                                             '<br><b>Region</b>: non-coding<br>',
                                                             '<br><b>Mean PEXT</b>: ',format(round(x[which(x$region=='non-coding'),i], 4), nsmall = 4),'<br>'))
                }          
                p[[i]] <-plotly_build(p[[i]]%>%
                                          layout(hovermode = 'x unified', barmode = 'stack', xaxis = list(
                                              title = "",
                                              showgrid = FALSE,
                                              showticklabels = FALSE  # Hide axis labels and ticks
                                          ),
                                          yaxis = list(type = 'log', range=c(-5.1,0.1), nticks=6,exponentformat = "e",
                                                       title = ""
                                          ))%>%
                                          add_annotations(
                                              text = ~I(k),
                                              x = 0.1,
                                              y = 1+((length(unique(c(input$Groups1,input$Groups2,input$Groups3)))-3)/10),
                                              yref = "paper",
                                              xref = "paper",
                                              xanchor = "center",
                                              yanchor = "top",
                                              showarrow = FALSE,
                                              font = list(size = 15)
                                          ))
            }}
        
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
                x = seq(shape$x0, shape$x1),
                y = rep(1,shape$x1-shape$x0+1),
                showlegend=F,
                color=shape$line$color,
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
                                              transcripts[which(transcripts$gene_name==input$gene),2])), 
                   margin=list(t=50))
        
    }, ignoreNULL = T)
    
    output$plot<-renderPlotly({plot()})
    text<-eventReactive(input$submit,{
        gn=transcripts[which(transcripts$gene_name==input$gene),]
        HTML("Number of reads = ", gn$nreads, 
             "<br>Number of unique transcript isoforms = ", gn$ntranscripts, 
             "<br>Number of full-splice match transcripts = ", gn$nFSM, " (", gn$nFSMreads, ' reads)', 
             "<br>Number of incomplete-splice match transcripts = ", gn$nISM, " (", gn$nISMreads, ' reads)',
             "<br>Number of novel in catalog transcripts = ", gn$nNIC, " (", gn$nNICreads, ' reads)',
             "<br>Number of novel not in catalog transcripts = ", gn$nNNC, " (", gn$nNNCreads, ' reads)',
             "<br><p> See <a href='https://github.com/ConesaLab/SQANTI3/wiki/SQANTI3-isoform-classification:-categories-and-subcategories' target='_blank'>here</a> for details on the classification groups </p>")}, ignoreNULL = T)
    output$text<-renderText({text()})
}

# Run the application 
shinyApp(ui = ui, server = server)

