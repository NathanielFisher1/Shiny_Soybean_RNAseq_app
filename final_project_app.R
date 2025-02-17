# Load packages
library(shiny)  
library(bslib)  
library(ggplot2)  
library(colourpicker)  
library(emojifont)  
library(vroom) 
library(tools)  
library(shinyjs)  
library(dplyr) 
library(DT)  
library(plotly)  
library(tidyverse)  
library(tidyr)  
library(devtools)  
library(SummarizedExperiment)  
library(matrixStats)  
library(gplots)  
library(DESeq2)  
library(ggthemes)  
library(ggbeeswarm) 
library(shinyWidgets) 
library(shinythemes)

# Load normalized counts matrix
testdf <- read.csv('./data/normalizedcountsmatrix_input.csv', header = TRUE, row.names = 1)

# List gene names for the individual gene plotting
searchable <- rownames(testdf) #for last panel

# UI
ui <- fluidPage(tags$style(type="text/css",
                           ".shiny-output-error { visibility: visible; }",
                           ".shiny-output-error:before { visibility: visible; }"),

theme = shinythemes::shinytheme("cyborg"),
#shinythemes::themeSelector(), # I removed this, too distracting

# Panel 1 is title panel
  titlePanel(HTML("<h3>Soybean (<i>Glycine max</i>) RNA-seq Data Visualizer</h3>")),
  tabsetPanel(
    tabPanel("Sample Metadata", 
             sidebarLayout(
              sidebarPanel(
                style = "position: relative; height: 800px; overflow-y: scroll;",
                           tags$h6("These data were gathered by Brown et. al 2015 and can be accessed", 
                                   HTML('<a href="https://rdrr.io/bioc/bigPint/man/se_soybean_cn_sub.html" target="_blank">here</a>.')),
                           tags$br(),
                           tags$h6("Cotyledon samples were collected from soybean plants in triplicate at 3 different time points and RNA sequencing was performed on samples."),
                           tags$br(),
                           tags$h6("Reference for further information on methods:"),
                           tags$h6("Brown AV, Hudson KA (2015) Developmental profiling of gene expression in soybean trifoliate leaves and cotyledons. BMC Plant Biol 15:169"),       
                           tags$br(),
                           tags$div(
                             style = "display: block; margin-left: auto; margin-right: auto;",
                             tags$img(src = "soybean_img.jpg", height = "450px")), 
                           tags$h6("Seikei Zusetsu (1804)", style = "font-size: 14px; margin-top: 5px;")),
             
              # Main panel tabs and outputs
               mainPanel( 
               tabsetPanel(type = "tabs",
                           tabPanel("Summary", tableOutput("summary_table")),
                           tabPanel("Table", DTOutput("table")),
                           tabPanel("Plots", plotOutput("replicate_plot"),plotOutput("growth_stage_plot")))))),
    
    # Next is panel for analyzing trasncript ounts
    tabPanel("Transcript Counts",
             sidebarLayout(sidebarPanel(
               sliderInput("var_slide", "Select genes within or greater than variance percentile of:", min = 0, max = 100, post = '%', value = 50),
               sliderInput("nonzero_slide", "Select genes with number of samples non-zero at least:", min = 0, max = 9, value = 9)),
             
             # Make tabs for panel
             mainPanel(
               tabsetPanel(type = "tabs",
                           
                           # Summary table
                           tabPanel("Summary", tableOutput("counts_summary_table")),
                           
                           # Make tab for plots of filtering based on mising samples and p values
                           tabPanel("Scatter Plots", 
                                    plotOutput("med_v_var"),
                                    plotOutput("med_v_nzeros")),
                           
                           # Make heatmap tab
                           tabPanel("Heatmap", plotOutput("heatmap")),
                           
                           # Make PCA plot tab
                           tabPanel("PCA",
                                    selectInput("pc_select1", "Select PC for X axis", c("PC1", 
                                                                                        "PC2", 
                                                                                        "PC3", 
                                                                                        "PC4", 
                                                                                        "PC5", 
                                                                                        "PC6", 
                                                                                        "PC7", 
                                                                                        "PC8", 
                                                                                        "PC9")),
                                    
                                    selectInput("pc_select2", "Select PC for Y axis", c("PC1", 
                                                                                        "PC2", 
                                                                                        "PC3", 
                                                                                        "PC4", 
                                                                                        "PC5", 
                                                                                        "PC6", 
                                                                                        "PC7", 
                                                                                        "PC8", 
                                                                                        "PC9")),
                                    
                                    selectInput("pca_fill", "Select Metadata to Group by", c("Sample_Name","Replicate", "Growth_Stage")),
                                    
                                     plotOutput("pca_scatter")))))),
    
    # Third panel is for DE analysis
    tabPanel("Differential Expression", 
             
             # Interface to select comparisons for volcano plots and p-value threshold for volcano plots
             sidebarLayout(sidebarPanel(
               radioButtons("var1", "Choose Comparison:", choices = c("Early_vs_Middle", "Early_vs_Late", "Middle_vs_Late")),
               radioButtons("button1", "Choose the column for the x-axis:", choices = c("logFC","logCPM","Likelihood_Ratio","p_value","adjusted_p_value")), #radiobuttons for x
               radioButtons("button2", "Choose the column for the y-axis:", choices = c("logFC","logCPM","Likelihood_Ratio","p_value","adjusted_p_value"), selected = "p_value"), #radiobuttons for y
               colourpicker::colourInput("col1","Select base point color:","#CC79A7"), #input color 1
               colourpicker::colourInput("col2","Select highlight point color:","#009E73"), #input color 2
               sliderInput("slide", "Select the magnitude of the FDR adjusted p-value coloring:", min = -3.5, max = 0, value = -1, step = .05) #slider for p-value              
             ),
               
             # Tab is for summary table of DE results
             mainPanel(
               tabsetPanel(
                 type = 'tabs',
                 tabPanel("DE Summary", DTOutput("de_summary_table")),
                 tabPanel("Plot", plotOutput("volcano")))))),
    
    # Tab is for volcano plot of DE results
    tabPanel("Visualize Inidividual Genes", 
            sidebarLayout(sidebarPanel(
              radioButtons("field1", "Please choose a field to group samples by:", choices = c("Sample_Name","Replicate", "Growth_Stage")),
                selectizeInput(
                  inputId = "searchme", 
                  label = "Search a gene of interest:",
                  multiple = FALSE,
                  choices = c("Search Bar" = "", searchable),
                  options = list(
                    create = FALSE,
                    placeholder = "e.g. Glyma05g05880.1",
                    maxItems = '1',
                    onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
                    onType = I("function (str) {if (str === \"\") {this.close();}}"))),
              radioButtons("field2", "Choose plot type:", choices = c('bar', 'boxplot', 'violin', 'beeswarm')),
              actionButton("plot_button", "Generate Plot", class = "btn-success")
            ),mainPanel(suppressMessages(plotOutput("individual_gene_plot")))))))

# Define server logic
server <- function(input, output, session) {
  cat(file=stderr(), "Server function started\n")
  
  # Read in summary data
  summary_table_data <- read.csv("./data/metadata_input.csv", header = T)
  
  # Display data as a table
  draw_summary_table <- function(inputdf){
    inputdf <- mutate(inputdf, growth_stage = as.factor(growth_stage), replicate = as.factor(replicate))
    final <- data.frame("Column_Name" = colnames(inputdf),
                        "Type" = sapply(inputdf[2,], FUN = function(x) class(x)),
                        "Distinct_Values" = character(3))
    final$Distinct_Values[1] <- paste(unique(inputdf$sample),collapse=', ')
    final$Distinct_Values[2] <- paste(unique(inputdf$growth_stage),collapse=', ')
    final$Distinct_Values[3] <- paste(unique(inputdf$replicate),collapse=', ')
    return(final)
  }
  
  # Make simple barplots for replicates
  make_barplot_replicate <- function(inputdf){
    plot <- ggplot(inputdf,aes(x = replicate)) + geom_bar(fill = "#CC79A7", color = "black", width = 0.7, size = 1)+
      labs(title = "Replicates", x = "Replicate", y = "Count") + theme_minimal() + theme(legend.position = "none") +
      theme(plot.background = element_rect(fill = "#151515", color = 
            axis.text.x = element_text(size = 18, color = "white"),  
            axis.text.y = element_text(size = 18, color = "white"),  
            axis.title.x = element_text(size = 24, color = "white"), 
            axis.title.y = element_text(size = 24, color = "white"), 
            plot.title = element_text(size = 30, color = "white"),
            legend.text = element_text(size = 24, color = "white"),
            legend.title = element_text(size = 24, color = "white"),
            panel.grid.major.x = element_blank(),  
            panel.grid.minor.x = element_blank())   
    return(plot)
    
  }
  
  # Make simple barplot for growth stages
  make_barplot_growth_stage <- function(inputdf){
    plot <- ggplot(inputdf,aes(x = growth_stage)) + geom_bar(fill = "#009E73", color = "black", width = 0.7, size=1)+
      labs(title = "Growth Stages", x = "Growth Stage", y = "Count") + theme_minimal() + theme(legend.position = "none") +
      theme(plot.background = element_rect(fill = "#151515", color = 
            axis.text.x = element_text(size = 18, color = "white"),  
            axis.text.y = element_text(size = 18, color = "white"),  
            axis.title.x = element_text(size = 24, color = "white"), 
            axis.title.y = element_text(size = 24, color = "white"), 
            plot.title = element_text(size = 30, color = "white"),
            legend.text = element_text(size = 24, color = "white"),
            legend.title = element_text(size = 24, color = "white"),
            panel.grid.major.x = element_blank(),  
            panel.grid.minor.x = element_blank())  
    return(plot)
  }
  
  # Making outputs for first tab
    output$summary_table <- renderTable(draw_summary_table(summary_table_data))#output table
    output$table <- renderDT(datatable(summary_table_data,options = list(ordering = TRUE)))
    output$replicate_plot <- renderPlot(make_barplot_replicate(summary_table_data))
    output$growth_stage_plot <- renderPlot(make_barplot_growth_stage(summary_table_data))
    
    # All below is for counts tab
    counts_data <- read.csv("./data/normalizedcountsmatrix_input.csv", header = TRUE,row.names = 1)
    
    # Make summary table
    make_count_summary_table <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf
      num_sample <-ncol(dff)
      num_genes <- nrow(dff)
      
      filtered_df <- as_tibble(dff) %>%
        mutate(row_var = rowVars(as.matrix(dff)),
               percentile = percent_rank(row_var)*100) %>%
        filter(rowSums(. != 0) >= slidernonzero) %>% # Need to add 1 because rownames are nonzero
        filter(percentile >= slidervar)
      
      num_genes_pass <- nrow(filtered_df)
      pct_genes_pass <- num_genes_pass / num_genes*100
      num_genes_fail <- num_genes - num_genes_pass
      pct_genes_fail <- 100 - pct_genes_pass

      final <- data.frame("Number of Samples" = num_sample,
                          "Total Number of Genes" = num_genes,
                          "Number of Genes Passing Filter" = num_genes_pass,
                          "Number of Genes Failing Filter" = num_genes_fail,
                          "Percent Genes Passing Filter" = pct_genes_pass,
                          "Percent Genes Failing Filter" = pct_genes_fail, check.names = FALSE)
      return(final)
    }
    
    # Scatterplots
    make_scatter_med_var <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf
      filtered_df <- as_tibble(dff) %>%
        mutate(row_var = rowVars(as.matrix(dff)),
               percentile = percent_rank(row_var)*100,
               nonzero = rowSums(. != 0),
               row_med = rowMedians(as.matrix(dff))) 
      filtered_df$Threshold <- ifelse(filtered_df$percentile >= slidervar & filtered_df$nonzero >= slidernonzero, 'TRUE', 'FALSE')
      plot <- ggplot(filtered_df, aes(x = log(row_med), y = log(row_var), color = Threshold)) + theme_dark() +
        scale_color_manual(name ="Within Threshold", values = c('FALSE' = "#CC79A7", 'TRUE' = "#009E73"))+
        geom_point() +
        labs(x = "log10(Median Count)", y = "log10(Count Variance)", title = "Median Count vs. Count Variance") + 
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
        
      return(plot)
    }
    
    make_scatter_med_zeros <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf
      filtered_df <- as_tibble(dff) %>%
        mutate(row_var = rowVars(as.matrix(dff)),
               percentile = percent_rank(row_var)*100,
               nonzero = rowSums(. != 0),
               row_med = rowMedians(as.matrix(dff))) 
      filtered_df$Threshold <- ifelse(filtered_df$percentile >= slidervar & filtered_df$nonzero >= slidernonzero, 'TRUE', 'FALSE')
      plot <- ggplot(filtered_df, aes(x = log(row_med), y = 9-nonzero , color = Threshold)) + theme_dark() +
        scale_color_manual(name ="Within Threshold", values = c('FALSE' = "#CC79A7", 'TRUE' = "#009E73"))+
        geom_point() +
        labs(x = "log10(Median Count)", y = "# of Zero Values/row", title = "Median Count vs. # of Zero Values/row") +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
      return(plot)
    }
    
    # Make Heatmap
    make_heatmap <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf # Need to assign variable to reactive object or it does not work 
      filtered_df <- as_tibble(dff, rownames = NA) %>% # Filtering while preserving rownames
        rownames_to_column() %>%
        mutate(row_var = rowVars(as.matrix(dff)), # Defining row
               percentile = percent_rank(row_var)*100) %>%
        filter(rowSums(. != 0) -3 >= slidernonzero) %>% # Need to subtract 3 because I added 3 nonzero rows
        filter(percentile >= slidervar) %>%
        select(-row_var,-percentile)

      input <- as.data.frame(filtered_df[,-1])
      rownames(input) <- filtered_df$rowname

      heat_data <- as.matrix(input)
      par(mar=c(5.1,4.1,4.1,4))
      heatmap.2(heat_data, scale = "column", trace = 'none', adjRow = c(.1,0), margins = c(10,10))
    }
    
    # Make PCA plots
    plot_pca <- function(inputdf,input1_pc,input2_pc, input3_fill) {
      input <- inputdf 
      pca <- prcomp(t(input)) 
      plotting_data <- data.frame('Sample_Name' = colnames(input),
                                  'Growth_Stage' = c(rep("early",3),rep("middle",3),rep("late",3)),
                                  'Replicate' = as.character(rep(c(1,2,3),3)),
                                  pca$x)
      pct_var1 <- round(pca$sdev[as.numeric(str_extract_all(input1_pc, "\\d")[[1]])]^2/sum(pca$sdev^2)*100, digits =1)
      pct_var2 <- round(pca$sdev[as.numeric(str_extract_all(input2_pc, "\\d")[[1]])]^2/sum(pca$sdev^2)*100, digits =1)
      

      
      plot <- ggplot(plotting_data,aes_string(x = input1_pc, y = input2_pc, color = input3_fill)) + 
        geom_point(size = 4) + labs(x = paste(input1_pc,', % Variance Explained:',pct_var1) , y = paste(input2_pc,', % Variance Explained:',pct_var2)) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  
              axis.text.y = element_text(size = 18, color = "white"),  
              axis.title.x = element_text(size = 24, color = "white"), 
              axis.title.y = element_text(size = 24, color = "white"), 
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.minor.x = element_blank())   
      return(plot)
    }
    
  # Make tab outputs
  output$counts_summary_table <- renderTable(make_count_summary_table(counts_data,input$var_slide,input$nonzero_slide))
  output$med_v_var <- renderPlot(make_scatter_med_var(counts_data,input$var_slide,input$nonzero_slide))
  output$med_v_nzeros <- renderPlot(make_scatter_med_zeros(counts_data,input$var_slide,input$nonzero_slide))
  output$heatmap <- renderPlot(
                               make_heatmap(counts_data,input$var_slide,input$nonzero_slide), height = 800)
  output$pca_scatter <- renderPlot(plot_pca(counts_data,input$pc_select1,input$pc_select2,input$pca_fill))
    
  
  # DE Panel 
  
  # Read in Data
  DE_summary_table_data <- read.csv("./data/differentialexpression_input.csv", header = TRUE)

  # Volcano plot
  volcano_plot <-
    function(dataf,comparison, x_name, y_name, slider, color1, color2) { 
      comp = ifelse(comparison == "Early_vs_Middle",'S1_S2',ifelse(comparison == "Early_vs_Late",'S1_S3', 'S2_S3'))
      x_final = ifelse(x_name == "logFC", "logFC", 
                       ifelse(x_name == "logCPM", "logCPM", 
                              ifelse(x_name == "Likelihood_Ratio", "LR", 
                                     ifelse(x_name == 'p_value', "PValue", 
                                            ifelse(x_name == 'adjusted_p_value', 'FDR', NULL)))))
      y_final = ifelse(y_name == "logFC", "logFC", 
                       ifelse(y_name == "logCPM", "logCPM", 
                              ifelse(y_name == "Likelihood_Ratio", "LR", 
                                     ifelse(y_name == 'p_value', "PValue", 
                                            ifelse(y_name == 'adjusted_p_value', 'FDR', NULL)))))
      
      dataf <- na.omit(dataf) # Remove NAs
      xval <- as.numeric(unlist(dataf[paste0(comp,'.',x_final)])) # Make vector of xvals
      yval <- as.numeric(unlist(dataf[paste0(comp,'.',y_final)])) # Make vector of yvals
      dataf$pointcolor <- ifelse(dataf[paste0(comp,'.','FDR')]>1*10^slider,'FALSE','TRUE') # Make new column to store color values
      p1 <- ggplot(dataf,aes(x=xval,y=-log10(yval), color = pointcolor)) + theme_dark()+ # Make plot
        geom_point()+
        scale_color_manual(name =paste0("FDR adjusted p-value < 1 x 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2)) +
        theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.position = "bottom")+
        labs(x = x_name, y= paste0("-log10(",y_name,")")) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.minor.x = element_blank(),
              legend.position = "bottom")   # Remove minor vertical grid lines
      return(p1)
    }
  
  
  # Make Outputs
  output$de_summary_table <- renderDT(datatable(DE_summary_table_data,options = list(ordering = TRUE))%>%
                                        formatRound(c(2:16), 3))
  output$volcano <- renderPlot(volcano_plot(DE_summary_table_data,input$var1, input$button1, input$button2, input$slide, input$col1, input$col2),
                               width = 800, height = 600)
  
# Individual gene panel
  
# Gene name is input$searchme
# Plot type is input$field2
# Category is input$field1


  # Make individual gene plot function
  
  make_final_plot <- function(inputdf, gene, category, plot_type){
    
    req(input$plot_button)  # Wait for button click before plotting
    
    
    for_plot <- data.frame('Sample_Name' = suppressWarnings(colnames(inputdf)),
                           'Growth_Stage' = factor(
                             c(rep("early",3), rep("middle",3), rep("late",3)),
                             levels = c("early", "middle", "late"),  # Factor growth stage
                             ordered = TRUE
                           ),
                           'Replicate' = as.character(rep(c(1,2,3),3)),
                           t(inputdf[gene,]))
    
    # For barplot
    make_barplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_bar(stat = 'identity', fill = 'blue') + labs(title = gene, y = 'Count') + theme_dark() +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
      return(p)
    }
    
    # For boxplot
    make_boxplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_boxplot(fill = "blue") + labs(title = gene, y = 'Count') + theme_dark() +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
      return(p)
    }
    
    # For violinplot
    make_violinplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_violin(fill = "blue") + labs(title = gene, y = 'Count') + theme_dark() +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
      return(p)
    }
    
    # For beeswarmplot
    make_beeswarmplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_beeswarm(color = "blue") + labs(title = gene, y = 'Count') + theme_dark() +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "#151515", color = "#151515"),
              axis.text.x = element_text(size = 18, color = "white"),  # Increase x-axis text size
              axis.text.y = element_text(size = 18, color = "white"),  # Increase y-axis text size
              axis.title.x = element_text(size = 24, color = "white"), # Increase x-axis title size
              axis.title.y = element_text(size = 24, color = "white"), # Increase y-axis title size)
              plot.title = element_text(size = 30, color = "white"),
              legend.text = element_text(size = 24, color = "white"),
              legend.title = element_text(size = 24, color = "white"),
              panel.grid.major.x = element_blank(),  # Remove major vertical grid lines
              panel.grid.minor.x = element_blank())   # Remove minor vertical grid lines
      return(p)
    }
    
    # Plotting based on radio buttons selected
    plot <- if(plot_type == 'bar'){
      make_barplot(for_plot, category, gene)
    } else if (plot_type == 'boxplot'){
      make_boxplot(for_plot, category, gene)
    } else if (plot_type == 'violin'){
      make_violinplot(for_plot, category, gene)
    } else {make_beeswarmplot(for_plot, category, gene)}
    return(plot)
  }
  
  
  # Assign output
  output$individual_gene_plot <- suppressWarnings(renderPlot({
    input$plot_button
    isolate(make_final_plot(counts_data,input$searchme,input$field1,input$field2))
    
  }))
    

  cat(file=stderr(), "Server function completed\n")
  
  
}

# Create Shiny app object

shinyApp(ui = ui, server = server)

