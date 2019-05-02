library(rsconnect)
library(shiny)
library(plotly)
library(dplyr)

# Inputs ------------------------------------------------------------------


#File inputs

metadata <-
  read.csv(file = "Miseq01_metadata_ITS_table.csv",
           sep = ",",
           stringsAsFactors = F)
abundance <-
  read.csv(file = "Miseq01_abundance_ITS_table.csv",
           sep = ",",
           stringsAsFactors = F)
taxonomy <-
  read.csv(file = "Miseq01_taxonomy_ITS_table.csv",
           sep = ",",
           stringsAsFactors = F)


# User Interface ----------------------------------------------------------


# Define UI for app that draws a histogram 
ui <- fluidPage(
  # set up 2 tabs, one for Fungal Abundance by Sample Type
  # one for Sample Types by Fungal Taxonomy
  tabsetPanel(
    tabPanel("Fungal Taxa by Sample Type", fluid = TRUE,
             
             # Sidebar layout with input and output definitions
             sidebarLayout(
               # Sidebar panel for inputs
               sidebarPanel(
                 # Input: Select a sample_type (later sampleType)
                 selectInput(
                   inputId = "sampleType",
                   label = "Select a sample type to see what fungi are in it:",
                   choices = sort(unique(metadata$sample_type))
                 ),
                 selectInput(
                   inputId = "taxonLevel",
                   label = "Select the fungal taxonomic level that you want to view:",
                   choices = colnames(taxonomy)[3:9]
                 ),
                 p("Bubble size indicates the number of sequences"),
                 p("Hover over the bubbles for more information")
               ),
               
               # Main panel for displaying outputs
               mainPanel(# Output: Plotly
                 plotlyOutput("metaPlot"))
             )
    ),
    tabPanel("Sample Types by Fungal Taxa", fluid = TRUE,
             
             # Sidebar layout with input and output definitions
             sidebarLayout(
               # Sidebar panel for inputs
               sidebarPanel(
                 # Input: Select a sampleType
                 selectInput(
                   inputId = "fungalTaxon",
                   label = "Select a fungal genus to see which sample types it's in:",
                   choices = sort(unique(taxonomy$Genus))
                 ),
                 p("Bubble size indicates the number of sequences"),
                 p("Hover over the bubbles for more information")
               ),
               
               # Main panel for displaying outputs
               mainPanel(# Output: Plotly
                 plotlyOutput("taxPlot"))
             )
    )))




# Code To Run -------------------------------------------------------------
# Run this code every time the user changes the value of the dropdown

server <- function(input, output) {
  
  

# Fungal Taxa By Sample Type

 
  filteredMetadata <- reactive({
    
    # ############# For TESTING
         # input<- list()
         # input$sampleType <- "Mushroom"
         # input$taxonLevel <- "Genus"
    # # 
  ## filter all tables by sample_type
    sample_filter    <- input$sampleType
    taxonomic_filter <- input$taxonLevel
    
    # filter metadata by sample_type
    filteredMeta <- metadata[metadata$sample_type == sample_filter, ]
    
    # filters abundance table to keep rows where match samples in filteredMeta and pulls out id
    filteredAbun <- abundance[abundance$Group %in% filteredMeta$id, ]
    abundIDs     <- filteredAbun[[1]]
    filteredAbun <- filteredAbun[2:ncol(filteredAbun)]
    # removes empty columns from abundance table (i.e. no reads)
    filteredAbun  <-
      filteredAbun[(colSums(filteredAbun, na.rm = T) != 0)]
    
    #filters taxonomy table to only include rows where OTUs match nonZeroAbun
    filteredTax <- taxonomy[taxonomy$OTU %in% colnames(filteredAbun), ]
   
  ## generate totals
  
  # get total reads for each OTU per sample
  totals <- colSums(filteredAbun)
  # get number of samples in which each OTU is present
  sampleCounts <- colSums(filteredAbun != 0)
  # combine into dataframe
  forPlotly <-
    cbind(Total = totals,
          Counts = sampleCounts,
          OTU = colnames(filteredAbun))
  # join with taxonomy
  
  forPlotly <-  forPlotly %>%
    as.data.frame(stringsAsFactors = F) %>%
    left_join(filteredTax, by = c("OTU"="OTU_ID")) 
  # remove unused taxonomy columns
  forPlotly <- forPlotly[colnames(forPlotly) %in% c("Total", "Counts", input$taxonLevel)]  
  # group by user selected taxonomic level and rename column "taxonLevel"
  colnames(forPlotly)[3] <- "taxonLevel"
  forPlotly <- group_by(forPlotly,taxonLevel)
  class(forPlotly$Total)  <- "numeric"
  class(forPlotly$Counts) <- "numeric"

  
  # remove unidentified seqs
  forPlotly <- forPlotly[!is.na(forPlotly$taxonLevel),]
  
  # summarize totals and counts by taxonomic group
  forPlotly <- summarise(forPlotly,
                         Total = sum(Total),
                         Counts = max(Counts))
  class(forPlotly$Total)  <- "numeric"
  class(forPlotly$Counts) <- "numeric"
  
  forPlotly

  })
  
 ## Generate Plot  
    output$metaPlot <- renderPlotly({
      
      p <- plot_ly(
        filteredMetadata(),
        x = ~ Counts,
        y = ~ Total,
        type = 'scatter',
        mode = 'markers',
        size = ~ Total,
        color = ~ taxonLevel,
        colors = 'Paired',
        #Choosing the range of the bubbles' sizes:
        sizes = c(10, 50),
        marker = list(opacity = 0.5, sizemode = 'diameter')
        
      ) %>%
        
        layout(title = 'Fungal Abundance',
               xaxis = list(showgrid = FALSE, title = "Number of Distinct Samples"),
               yaxis = list(showgrid = FALSE, title = "Total Reads", type = "log"),
               showlegend = FALSE)
    })
  
    

# Sample Types By Fungal Taxon --------------------------------------------
filteredTaxdata <- reactive({
  ####### FOR TESTING
   #  input <-list()
   #  input$fungalTaxon <- "Aspergillus"
   # # 
  ## filter by taxon
  # filter taxonomy by user selected taxon
  filteredTax2 <- taxonomy[taxonomy$Genus == input$fungalTaxon,]
  filteredTax2 <- filteredTax2[!is.na(filteredTax2$Genus),] 
  
  # filter abundance by selected taxon, first pulling out "Group" column
  abundanceIDs2 <- abundance$Group
  filteredAbun2 <- abundance[colnames(abundance) %in% c(filteredTax2$OTU)]
  # remove empty rows and corresponding abundanceIDs (= Group column)
  toRemove <- rowSums(filteredAbun2)>0
  filteredAbun2 <- filteredAbun2[toRemove,]
  abundanceIDs2 <- abundanceIDs2[toRemove]
  
  # filter metadata by samples in which taxon is present
  filteredMeta2 <- metadata[metadata$id %in% abundanceIDs2,]
  
  ## Generate Totals
  # abundance totals by sample
  totals2    <- rowSums(filteredAbun2)
  forPlotly2 <- cbind(Total2 = totals2, id = abundanceIDs2) %>% as.data.frame()
  
  # join with sample metadata
  forPlotly2 <- left_join(forPlotly2,filteredMeta2, by = "id")
  
  # generate counts for each sample type
  forPlotly2 <- group_by(forPlotly2,sample_type)
  forPlotly2 <- summarise(forPlotly2, Total = sum(Total2), Counts = n())
  forPlotly2 <- forPlotly2[-nrow(forPlotly2),]
  forPlotly2
  
})



output$taxPlot <- renderPlotly({
  
  p <- plot_ly(
    filteredTaxdata(),
    x = ~ Counts,
    y = ~ Total,
    type = 'scatter',
    mode = 'markers',
    size = ~ Total,
    color = ~ sample_type,
    colors = 'Paired',
    #Choosing the range of the bubbles' sizes:
    sizes = c(10, 50),
    marker = list(opacity = 0.5, sizemode = 'diameter')
    
  ) %>%
    
    layout(title = 'Fungal Abundance',
           xaxis = list(range(0,20), showgrid = FALSE, title = "Number of Distinct Samples"),
           yaxis = list(showgrid = FALSE, title = "Total Reads", type = "log"),
           showlegend = FALSE)
})
}
shinyApp(ui = ui,server = server)
