library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(jsonlite)
library(plotly)
library(ggplot2)
library(readr)
library(reactable)
library(shinyWidgets)

date_who=readRDS('data/date_who.rds')
date_voc=readRDS('data/date_voc.rds')


ui <- dashboardPage(
  title = "VIVID - Virus Variant Dissector",
  dashboardHeader(title =  span(HTML(paste0(tags$u('Vi'), 'rus ', tags$u('V'), 'ar', tags$u('i'), 'ant ', tags$u('D'), 'issector')),
                                style = "font-size: 16px")),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      #menuItem("Highlights", tabName = "summary", icon = icon("highlighter")),
      menuItem("VOC epidemiology", tabName = "epi", icon = icon("chart-bar")),
      menuItem("WHO weekly update", tabName = "who", icon = icon("newspaper")),
      menuItem("Fitness prediction", tabName = "fitness", icon = icon("virus")),
      menuItem("Lineage to sequence", tabName = "seq1", icon = icon("dna")),
      menuItem("Sequence to lineage", tabName = "seq2", icon = icon("sitemap")),
      menuItem("Resources", tabName = "resources", icon = icon("list")),
      menuItem("Feedback", icon = icon("comment"), 
               href = "mailto:ds@rvacmed.com?subject=Feedback for the COVID Monitor")
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
           #shiny-tab-resources li {
                font-size: 2rem;
                margin-bottom: 0.8rem; 
           }
      "))
    ),
    tabItems(
      
      tabItem(tabName = "epi",
              h1("Monitor the designated variants of concern (VOCs)"),
              selectInput("Date_voc", label = "Report Date", 
                          choices = date_voc, selected=max(date_voc) ),
              h3("Designated VOCs"),
              p("Current VOCs retrieved from ", 
                tags$a(target="_blank", 
                       href="https://outbreak.info/situation-reports#voc",
                       "outbreak.info - VOC Reports"), 
                ". VOCs were designated by outbreak.info, CDC, ECDC, WHO, 
                or Public Health England. Only the earliest designataion dates 
                from these institutes are reported here."),
              reactableOutput("voc_display"),
              
              h3("Prevalance of the VOCs"),
              p("The daily prevalence is calculated as the proportion of 
                the VOC sequences from all sequences collected on that day."),
              plotlyOutput("voc_pre"),
              
              h3("Recently prevalent VOC(s)"),
              p("The prevalent VOC(s) in recent month (defined as 
                peak prevalence > the selected threshold within 30 days from report date)"),
              sliderTextInput(
                inputId = "pre_voc", label = "Prevalence threshold", 
                choices = c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), selected = 0.02, 
                grid = TRUE
              ),
              reactableOutput("recent_voc"),
              
              h3("Recently prevalent pangolin lineages of VOCs"),
              p("The prevalent sublineages in recent month (defined as 
                peak prevalence > the selected threshold within 30 days from report date)"),
              sliderTextInput(
                inputId = "pre_sub", label = "Prevalence threshold", 
                choices = c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), selected = 0.02, 
                grid = TRUE
              ),
              reactableOutput("sub_display"),
              p(""),
              plotlyOutput("sub_pre"),
              p(""),
              p("The spike protein sequence mutations of 
                the recently prevalent pangolin lineage(s)."),
              plotlyOutput("sub_mut")
              
      ),
      
      tabItem(tabName = "who",
              h1("WHO weekly epidemiological update"),
              selectInput("Date_WHO", label = "Report Date", 
                          choices = date_who, selected=max(date_who) ),
              uiOutput("pdfview")
              
      ),
      
      tabItem(tabName = "fitness",
              h1("Fitness prediction"),
              selectInput("Date_fit", label = "Report Date", 
                          choices = date_voc, selected=max(date_voc) ),
              h3("Predicted lineage fitness scores"),
              p("Warning: the algorithm failed to predict in current 'variant soup'.
                We are actively working on new prediction algorithms.",  style = "color: red;"),
              DT::DTOutput("fit_lineage"),
              p(""),
              h3("Predicted individual mutation fitness scores"),
              DT::DTOutput("fit_mutation"),
              p(""),
              #plotlyOutput("sub_mut"),
              p(""),
              h3("Note"),
              p("Predicting the emerging variants or mutations is much-needed but a challenging task. We implemented the 
                best publicly available machine learning algorithm here, with weekly retraining and prediction."),
              p("The algorithm was developed by ",
                tags$a(target="_blank", 
                       href="https://www.science.org/doi/10.1126/science.abm1208",
                       "Obermeyer et al 2022"),
                ". It is able to measure the relative fitness (R", tags$sub("0"), ") of lineages according to
                the spreading curve in population, and estimate the contribution of individual mutations to 
                the overall fitness by fitting a linear model, similar with genome-wide association studies (GWAS), 
                thus the individual mutation scores can also be used to predict the overall fitness of a new lineage. "),
              p("The prediction is generally accurate, except for lineages with many unseen mutations, e.g. the Omicron 
              lineages with much more mutations than Alpha, Beta, and Delta lineages could not be accurately 
                predicted before any omicron lineage appeared. But once the algorithm sees a few such new lineages, 
                it is able to quickly correct the prediction by retraining using the new data, before they become predominant. 
                Therefore, we retain the algorithm every week to keep the prediction up-to-date."
              )

      ),
      
      tabItem(tabName = "seq1",
              h1("Retrieve the DNA, RNA, and protein sequence of the spike gene by lineage"),
              selectizeInput('lineage', 'Select pango lineage(s)', NULL, 
                             multiple = T, selected='Reference'),
              selectizeInput('des', 'Include descendants', c('Yes','No'), selected='No'),
              selectizeInput('gen', 'Sequence generation method', 
                             c('Earliest', 'Latest',
                               'Most frequent DNA', 'Most frequent protein'),
                             selected='Earliest'),
              selectizeInput('seq_type', 'Sequence type', c('Protein', 'mRNA', 'cDNA'),
                             multiple = T, selected='Protein'),
              
              uiOutput('seq_out'),
              
              p(""),
              h3("Note"),
              p("The sequence IDs with designated lineages come from ", 
                tags$a(target="_blank", 
                       href="https://github.com/cov-lineages/pango-designation/blob/master/lineages.csv",
                       "pango-designation github"),
                "."
                ),
              p("If selecting \"Include descendants\", the sequence generation method will 
                consider the lineage and all its descendants, 
                e.g. BA.4 will include sequences from BA.4, BA.4.1, BA.4.1.1, and BA.4.1.2. 
                It's highly recommended to consult the data science team before selecting \"Yes\"."
              ),
              p("The following methods can be used to retrieve spike protein/RNA/DNA sequences:"),
              tags$ul(
                tags$li("\"Earliest\" - the earliest collected genome of the selected lineage 
                        (after removing sequences with non-ATCG bases)."),
                tags$li("\"Latest\" - the latest collected genome of the selected lineage 
                        (after removing sequences with non-ATCG bases)."),
                tags$li("\"Most frequent DNA\" - the most frequent spike cDNA sequence of the selected lineage 
                        (after removing sequences with non-ATCG bases)."),
                tags$li("\"Most frequent protein\" - the most frequent spike protein sequence of the selected lineage 
                (after removing sequences with non-standard amino acids). 
                        No DNA or RNA sequence will be generated if this is selected.")
              ),
              p('If multiple sequences from the same date or with the same frequency, 
                only distinct sequences are reported. Please note that sometimes the DNA sequences are distinct, 
                but the translated protein sequences are the same. In this case, all protein sequences are reported.'),
              p("The sequence frequencies are calculated from the subset of GISAID sequences listed in ",
                tags$a(target="_blank", 
                       href="https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json",
                       "pango-designation github"),
                "(about 10% of all GISAID sequences). We do not use all GISAID sequences because the lineage designation in 
                GISAID database is often lagged, especially for the recent lineages, 
                which may have zero genome assigned."
              ),
              p("The retrieved sequences may be different when using this tool 
                at a different time, because the sequence list and lineage designation are updated frequently.")
              
      ),
      
      tabItem(tabName = "seq2",
              h1("Identify the pango lineage from DNA, RNA or protein sequence"),
              textAreaInput('seq_input', "Input your sequence (only the sequence strings, no fasta header or annotations)", value = "",
                            width = '60%', rows=10),
              selectizeInput('seq_type2', 'Sequence type', c('Protein', 'mRNA', 'cDNA'),
                             multiple = F, selected='Protein'),
              actionButton('submit','Submit'),
              p(''),
              h3("Blast result"),
              verbatimTextOutput('lin_out') %>% withSpinner()
      ),
      
      tabItem(tabName = "resources",
              h1("A list of useful COVID resources"),
              tags$ul(
                tags$li(tags$a(target="_blank", 
                               href="https://www.gisaid.org/", 
                               "GISAID")),
                tags$li(tags$a(target="_blank", 
                               href="https://covid19.who.int/", 
                               "WHO - COVID-19 Dashboard")),
                tags$li(tags$a(target="_blank", 
                               href="https://www.who.int/activities/tracking-SARS-CoV-2-variants", 
                               "WHO - Tracking SARS-CoV-2 Variants")),
                tags$li(tags$a(target="_blank", 
                               href="https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports", 
                               "WHO - COVID-19 Weekly Epidemiological Update")),
                tags$li(tags$a(target="_blank", 
                               href="https://outbreak.info/situation-reports#voc", 
                               "outbreak.info - Variant of Concern Reports")),
                tags$li(tags$a(target="_blank", 
                               href="https://covariants.org/variants", 
                               "CoVariants - Variants and Mutations")),
                tags$li(tags$a(target="_blank", 
                               href="https://covariants.org/shared-mutations", 
                               "CoVariants - Shared mutations across variants")),
                tags$li(tags$a(target="_blank", 
                               href="https://cov-lineages.org/lineage_list.html", 
                               "Pango Lineage List")),
                tags$li(tags$a(target="_blank", 
                               href="https://github.com/cov-lineages/pango-designation/blob/master/pango_designation/alias_key.json", 
                               "Pango Lineage Hierarchy")),
                tags$li(tags$a(target="_blank", 
                               href="https://nextstrain.org/ncov/gisaid/global/6m", 
                               "Nextstrain - Real-time tracking of COVID genomic evolution")),
                tags$li(tags$a(target="_blank", 
                               href="https://clades.nextstrain.org/", 
                               "Nextclade - Clade assignment, mutation calling, and sequence quality checks")),
                tags$li(tags$a(target="_blank", 
                               href="https://cov-spectrum.org/", 
                               "CoV-Spectrum - Explore up-to-date genome data and monitor variants")),
                tags$li(tags$a(target="_blank", 
                               href="https://covid19dashboard.regeneron.com/", 
                               "Regeneron - COVID-19 dashboard")),
              )
              
      )
      
    )
  ),
  footer=dashboardFooter(left = "Copyright © 2022, RVAC Medicines. All rights reserved.",
                         right = tags$a(href="Terms.html", "Terms of Use"))
)



