library(shiny)
library(gmodels)

source("Import_DEG_sets.R") # for list_dfs_2

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20))



# Define UI ----
ui <- fluidPage(
  titlePanel("Good Biological Samples Volcano"),
  
  fluidRow(
    
    column(width = 2.5,
           selectInput("my_comparison",
                       label = "Volcano plots",
                       choices = df_names,
                       width = "40%"),
           
           # Add checkbox for removing the non coding RNAs
           checkboxInput("FilterOut_Rvnc", label = "Remove Rvnc genes", value = T),
           
           fluidRow(
             column(width = 3,
                    textInput("my_GeneID", 
                              label = "Gene ID (Will label gene orange)",
                              value = "Rv...")
             ),
             column(width = 2.5,
                    uiOutput("gene_link")  # New UI output for the link
             )
           ),
           # Add checkbox to toggle gene set points
           checkboxInput("show_gene_set", label = "Show gene sets", value = FALSE),
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           # Dropdown for selecting the gene set within the chosen rda file.
           # Start with no selection
           selectInput("my_GeneSet",
                       label = "Gene Set (Will label genes yellow)",
                       choices = NULL)
    ),
    
    column(width = 5,
           plotlyOutput("volcano_plot",
                        width = 1000, height = 600),
    ),
  )
  
)

# Define server logic ----
server <- function(input, output, session) {
  
  # Gene Link
  output$gene_link <- renderUI({
    req(input$my_GeneID)  # Ensure there's a valid input
    url <- paste0("https://mycobrowser.epfl.ch/genes/", input$my_GeneID)
    tags$a(href = url, target = "_blank", paste0("View Details of ", input$my_GeneID, " on Mycobrowser"))
  })
  
  # When a new gene set source is selected, update the gene set dropdown
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(session, "my_GeneSet",
                      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
                      selected = NULL)
  })
  
  # Volcano Plot
  output$volcano_plot <- renderPlotly({
    
    # Add data for labelling a single gene
    single_gene <- list_dfs_2[[input$my_comparison]] %>% 
      filter(GENE_ID == input$my_GeneID)
    
    # Add data for labelling a gene set
    gene_set <- list_dfs_2[[input$my_comparison]] %>%
      filter(GENE_ID %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Make new data name for plotting
    plot_data <- list_dfs_2[[input$my_comparison]]
    
    # Filter out the Rvnc genes if requested (made new dataframe to plot)
    if (input$FilterOut_Rvnc) {
      plot_data <- plot_data %>% 
        filter(!str_detect(GENE_ID, regex("Rvnc", ignore_case = T)))
    }

    # Make the Volcano Plot
    my_volcano <- plot_data %>%
      
      ggplot(aes(x = LOG2FOLD, y = -log10(AVG_PVALUE), col = DE, label = DE_labels, text = GENE_ID, label2 = GENE_NAME, label3 = PRODUCT)) + 
      geom_point() + 
      
      # Add a differently colored point
      geom_point() +
      # Conditionally add the gene set points (yellow) based on checkbox value
      { if(input$show_gene_set) 
        geom_point(data = gene_set, color = "yellow", aes(col = DE, label = DE_labels, text = GENE_ID))
        else 
          NULL } +
      geom_point(data = single_gene, color = "yellow", aes(col = DE, label = DE_labels, text = GENE_ID)) + 
      
      labs(title = input$my_comparison) + 
      geom_vline(xintercept = c(-1,1), col = "grey", linetype = "dashed") + 
      geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
      scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
    # geom_label_repel(max.overlaps = 10) # Can do geom_text_repel or geom_label_rebel
    
    # Determine the max and min axes values for labeling 
    plot_build <- ggplot_build(my_volcano)
    y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
    x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
    x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
    
    # Add the gene number annotations
    text_up <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant up") %>% nrow()
    text_down <- list_dfs_2[[input$my_comparison]] %>% filter(DE == "significant down") %>% nrow()
    my_volcano_annotated <- my_volcano +
      annotate("text", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
      annotate("text", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
    
    final_plot <- my_volcano_annotated + my_plot_themes 
    final_plot
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)