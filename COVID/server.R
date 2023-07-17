library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(jsonlite)
library(plotly)
library(ggplot2)
library(readr)
library(reactable)
library(dplyr)
library(Biostrings)
library(seqinr)

#LoadToEnvironment <- function(RData, env=new.env()) {
#  load(RData, env)
#  return(env)
#}

recent = 30  # define 'recent' days to look back
min_ratio = 0.02 # minimum threshold of prevalent variant/lineage
min_freq = 0.1  # minimum threshold of mutation frequency in a variant/lineage

ref_nuc=readRDS('data/ref_nuc.rds')
ref_prot=readRDS('data/ref_prot.rds')
lineages=readRDS('data/lineages.rds')
lineage_max_date_nuc=readRDS('data/lineage_max_date_nuc.rds')
lineage_min_date_nuc=readRDS('data/lineage_min_date_nuc.rds')
lineage_freq_max_nuc=readRDS('data/lineage_freq_max_nuc.rds')
lineage_freq_max_prot=readRDS('data/lineage_freq_max_prot.rds')

source('retrieve_seq.R')

server <- function(input, output, session) {
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if(!is.null(query$tab)) {
      updateTabsetPanel(session, 'tabs', query$tab)
    }
  })
  
  output$voc_display <- renderReactable({
    
    voc_report<-readRDS(paste0('data/voc_report_',input$Date_voc, '.RDS')) 

    voc = voc_report[voc_report$variantType %in% c('Variant of Concern','Previously Circulating Variant of Concern'),]
    
    # get the earliest VOC designation from multiple institutes
    voc$voc_name=ifelse(is.na(voc$short_name), voc$who_name, voc$short_name)
    
    voc$earliest_designation_date=unlist(lapply(voc$classifications, function(x){
      min(x$dateModified[x$variantType=='VOC'], na.rm = TRUE)
    }))
    
    outbreak_date<-as.character(as.Date(voc$classificationTable$VOC$outbreak$label, format="%d %b %Y"))
    outbreak_date<-pmax(voc$earliest_designation_date, outbreak_date, na.rm = TRUE)
    voc$earliest_designation_date<-outbreak_date
    
    voc_display =voc[,c("voc_name", "pangolin_lineage", "nextstrain_clade",
                        "gisaid_clade",  "location_first_identified",
                        "earliest_designation_date", "variantType")]
    
    reactable(voc_display, 
              defaultSorted = "earliest_designation_date", defaultSortOrder="desc")
    
    #env <- eventReactive(input$Date_voc, {
    #  reactiveFileReader(1000, session, paste0("data/all_variable_", input$Date_voc, '.RData'), LoadToEnvironment)
    #})
    #reactable(env()[['voc_display']])
  })
  
  output$voc_pre <-renderPlotly({
    voc_prevalence<-readRDS(paste0('data/voc_prevalence_',input$Date_voc, '.RDS'))
    plot_ly(voc_prevalence, x = ~date, y = ~proportion, color = ~lineage) %>%
      add_lines() %>% layout(title = 'Daily prevalance of the designated variants of concern')
  })
  
  output$recent_voc <-renderReactable({
    min_ratio<-input$pre_voc
    voc_prevalence<-readRDS(paste0('data/voc_prevalence_',input$Date_voc, '.RDS'))
    voc_prevalence_recent =voc_prevalence[voc_prevalence$date>=as.Date(input$Date_voc)-recent, ]
    prevalent_voc=unique(
      voc_prevalence_recent$lineage[voc_prevalence_recent$proportion>min_ratio])
    reactable(data.frame(prevalent_voc))
  })
  
  output$sub_display <- renderReactable({
    min_ratio<-input$pre_sub
    prevalent_sub_df<-readRDS(paste0('data/voc_sublineage_prevalence_',input$Date_voc, '.RDS'))
    sub_prevalence_recent =prevalent_sub_df[prevalent_sub_df$date>=as.Date(input$Date_voc)-recent, ]
    prevalent_sub=unique(
      sub_prevalence_recent$lineage[sub_prevalence_recent$proportion >min_ratio])
    
    prevalent_sub_df=prevalent_sub_df[prevalent_sub_df$lineage %in% prevalent_sub,]
    
    prevalent_sub_df_unique=unique(prevalent_sub_df[,c("lineage","who_name")])
    
    #prevalent_sub_df_unique<-readRDS(paste0('data/voc_sublineage_recent_',input$Date_voc, '.RDS')) 
    reactable(prevalent_sub_df_unique)
  })
  
  output$sub_pre <-renderPlotly({
    min_ratio<-input$pre_sub
    prevalent_sub_df<-readRDS(paste0('data/voc_sublineage_prevalence_',input$Date_voc, '.RDS'))
    sub_prevalence_recent =prevalent_sub_df[prevalent_sub_df$date>=as.Date(input$Date_voc)-recent, ]
    prevalent_sub=unique(
      sub_prevalence_recent$lineage[sub_prevalence_recent$proportion >min_ratio])
    
    prevalent_sub_df=prevalent_sub_df[prevalent_sub_df$lineage %in% prevalent_sub,]
    
    plot_ly(prevalent_sub_df, x = ~date, y = ~proportion, color = ~lineage) %>%
      add_lines() %>% layout(title = 'Daily prevalance of the recently prevalent sublineages')
  })
  
  output$sub_mut <-renderPlotly({
    min_ratio<-input$pre_sub
    
    prevalent_sub_df<-readRDS(paste0('data/voc_sublineage_prevalence_',input$Date_voc, '.RDS'))
    sub_prevalence_recent =prevalent_sub_df[prevalent_sub_df$date>=as.Date(input$Date_voc)-recent, ]
    prevalent_sub=unique(
      sub_prevalence_recent$lineage[sub_prevalence_recent$proportion >min_ratio])
    
    prevalent_mutation_sub2_s<-readRDS(paste0('data/voc_sublineage_mutation_',input$Date_voc, '.RDS'))
    
    prevalent_sub_df=prevalent_sub_df[prevalent_sub_df$lineage %in% prevalent_sub,]
    
    prevalent_sub_df_unique=unique(prevalent_sub_df[,c("lineage","who_name")])
    
    prevalent_mutation_sub2_s<-prevalent_mutation_sub2_s[prevalent_mutation_sub2_s$lineage %in% prevalent_sub |
                                                           prevalent_mutation_sub2_s$lineage %in% unique(prevalent_sub_df_unique$who_name),]
    
    plot_ly(prevalent_mutation_sub2_s,
            x=~reorder(mutation, codon_num), y=~lineage, z=~prevalence, 
            type = "heatmap", colors = "Greys",
            height= 30*length(unique(prevalent_mutation_sub2_s$lineage))+100,
            width = 12*length(unique(prevalent_mutation_sub2_s$mutation))+100) %>%
      layout(margin=list(l=50, r=50, t=50, b=50, pad = 4),showlegend=FALSE,
             xaxis = list(title = 'mutation'), 
             title = 'Mutation frequency of the recently prevalent sublineages')
  })
  
  output$pdfview <- renderUI({
    tags$iframe(style="height:800px; width:100%", src=paste0(input$Date_WHO, '.pdf'))
  })
  
  output$fit_lineage<-DT::renderDT({
    fitness_lineage=readRDS(paste0('data/fitness_lineage_',input$Date_fit, '.RDS'))
    DT::datatable(fitness_lineage)
  })
  
  output$fit_mutation<-DT::renderDT({
    fitness_mutation=readRDS(paste0('data/fitness_mutation_',input$Date_fit, '.RDS'))
    DT::datatable(fitness_mutation)
  })
  
  #output$fit_mutation<-renderPlotly({
  #  
  #})
  
  updateSelectizeInput(session, 'lineage', choices = lineages, selected="Reference", server = TRUE)
  
  output$seq_out <- renderUI({
    retrieve_sequence(input$lineage, input$des, input$gen, input$seq_type)
  })
  
  blast_result <- eventReactive(input$submit, {
    unlink("/tmp/output_seq.txt")
    write.fasta(sequences=as.list(input$seq_input), names='input_sequence',
                file.out= "/tmp/input_seq.fa")
    if(input$seq_type2=='Protein'){
      system('blastp -task blastp-fast -query /tmp/input_seq.fa -db data/lineage_prot -out /tmp/output_seq.txt -max_target_seqs 5')
    }else{
      system('blastn -task megablast -strand plus -query /tmp/input_seq.fa -db data/lineage_nuc -out /tmp/output_seq.txt -max_target_seqs 5')
    }
    
    read_file("/tmp/output_seq.txt")
  })
  
  output$lin_out <- renderText({
    blast_result()
  })
  
  

}
