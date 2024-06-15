#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (interactive()){

library(shiny) # for the app
library(bslib)
library(ggplot2) # for plots
library(colourpicker) # you might need to install this
library(emojifont) # for emojis
library(vroom) # for reading in tsv or csv files
library(tools) #just some tools in the toolbelt ;)
library(shinyjs) # for images
library(dplyr) #data manipulation
library(DT) # for sortable table
library(plotly) # for plotting
library(tidyverse) # for data wrangling
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

#just some initial inputs for the search bar
setwd('/Users/nathaniel_fisher/Desktop/BUMS/R')
testdf <- read.csv('normalizedcountsmatrix_input.csv', header = TRUE, row.names = 1)
searchable <- rownames(testdf) # for last box

# Define UI for application that draws a histogram
ui <- fluidPage(tags$style(type="text/css",
                           ".shiny-output-error { visibility: hidden; }",
                           ".shiny-output-error:before { visibility: hidden; }"
),
# setBackgroundColor(
#   color = c("orange", "yellow", "red"),
#   gradient = "linear",
#   direction = "bottom"
# ),
shinythemes::themeSelector(),
  titlePanel("Let's learn about soybeans!"),
  # Create tabs with content
  tabsetPanel(
    tabPanel("Samples", 
             sidebarLayout(sidebarPanel(
               titlePanel(title = fluidRow(icon = emoji("angry"), # Emoji 1
                                           emoji("whale2"),            # Emoji 2
                                           emoji("dolphin"),       # Emoji 3
                                           emoji("whale2"),
                                           emoji("dolphin"),
                                           emoji("whale2"),
                                           emoji("dolphin"),
                                           emoji("whale2"),
                                           emoji("dolphin"),
                                           emoji("whale2"),
                                           emoji("dolphin"))), 
               fileInput(inputId = 'summary_data_file', "Upload summary matrix in CSV format:") #set fileinput
               
             ),
             mainPanel( #main panel tabs and outputs
               tabsetPanel(type = "tabs",
                           tabPanel("Summary", tableOutput("summary_table")),
                           tabPanel("Table", DTOutput("table")),
                           tabPanel("Plots", plotOutput("replicate_plot"),plotOutput("growth_stage_plot")))))
    ),
    tabPanel("Counts",
             sidebarLayout(sidebarPanel(
               titlePanel(title = fluidRow(icon = emoji("angry"), # Emoji 1
                                           emoji("dragon"),            # Emoji 2
                                           emoji("snake"),       # Emoji 3
                                           emoji("dragon"),
                                           emoji("snake"),
                                           emoji("dragon"),
                                           emoji("snake"),
                                           emoji("dragon"),
                                           emoji("snake"),
                                           emoji("dragon"),
                                           emoji("snake")  )), #just adding emoji for fun
               fileInput(inputId = 'counts_file', "Upload normalized counts matrix in CSV format:"),
               sliderInput("var_slide", "Select genes within or greater than variance percentile of:", min = 0, max = 100, post = '%', value = 50),
               sliderInput("nonzero_slide", "Select genes with number of samples non-zero at least:", min = 0, max = 9, value = 9)
               
             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Summary", tableOutput("counts_summary_table")),
                           tabPanel("Scatter Plots", 
                                    plotOutput("med_v_var"),
                                    plotOutput("med_v_nzeros")
                                    ),
                           tabPanel("Heatmap", plotOutput("heatmap")),
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
                                    
                                     plotOutput("pca_scatter")))))
    ),
    tabPanel("DE", 
             sidebarLayout(sidebarPanel(
               titlePanel(title = fluidRow(icon = emoji("angry"), # Emoji 1
                                           emoji("gorilla"),            # Emoji 2
                                           emoji("monkey"),       # Emoji 3
                                           emoji("gorilla"),
                                           emoji("monkey"),
                                           emoji("gorilla"),
                                           emoji("monkey"),
                                           emoji("gorilla"),
                                           emoji("monkey"),
                                           emoji("gorilla"),
                                           emoji("monkey"))),
               fileInput(inputId = 'DE_file', "Upload results from differential expression analysis in CSV format:"),
               radioButtons("var1", "Choose Comparison:", choices = c("Early_vs_Middle", "Early_vs_Late", "Middle_vs_Late")),
               radioButtons("button1", "Choose the column for the x-axis:", choices = c("logFC","logCPM","Likelihood_Ratio","p_value","adjusted_p_value")), #radiobuttons for x
               radioButtons("button2", "Choose the column for the y-axis:", choices = c("logFC","logCPM","Likelihood_Ratio","p_value","adjusted_p_value"), selected = "p_value"), #radiobuttons for y
               colourpicker::colourInput("col1","Select base point color:","purple"), #input color 1
               colourpicker::colourInput("col2","Select highlight point color:","green"), #input color 2
               sliderInput("slide", "Select the magnitude of the FDR adjusted p-value coloring:", min = -3.5, max = 0, value = -1, step = .05) #slider for p-value              
             ),
               
             mainPanel(
               tabsetPanel(
                 type = 'tabs',
                 tabPanel("DE Summary", DTOutput("de_summary_table")),
                 tabPanel("Plot", plotOutput("volcano"))
               
               
             )))
        
    ),
    tabPanel("Visualize Inidividual Genes", 
            sidebarLayout(sidebarPanel(
              titlePanel(title = fluidRow(icon = emoji("angry"), # Emoji 1
                                          emoji("tropical_fish"),            # Emoji 2
                                          emoji("blowfish"),       # Emoji 3
                                          emoji("tropical_fish"),
                                          emoji("blowfish"),
                                          emoji("tropical_fish"),
                                          emoji("blowfish"),
                                          emoji("tropical_fish"),
                                          emoji("blowfish"),
                                          emoji("tropical_fish"),
                                          emoji("blowfish"))),
              radioButtons("field1", "Please choose a field to group samples by:", choices = c("Sample_Name","Replicate", "Growth_Stage")),
              #title = "Search Bar",
              #fluidRow(
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
            ),mainPanel(suppressMessages(plotOutput("individual_gene_plot"))
              
            )

              

            )
    )
  )
  
)

# Define server logic
server <- function(input, output, session) {

  #load summary data
  #   load_summary_table_data <- reactive({ #load data reactive, with validate statement to prevent warnings, return dataframe
  #   req(input$summary_data_file) #getinputfile
  #   ext <- tools::file_ext(input$summary_data_file$name)
  #   switch(ext,
  #          csv = vroom::vroom(input$summary_data_file$datapath, delim = ","),
  #          tsv = vroom::vroom(input$summary_data_file$datapath, delim = "\t"),
  #          validate("Invalid file; Please upload a .csv or .tsv file")
  #   )
  # })
  
  #all below is only for Samples tab
  load_summary_table_data <- reactive({ #load data reactive, with validate statement to prevent warnings, return dataframe
      infile <- input$summary_data_file
      validate(
        need(input$summary_data_file != "", "Please select a data set in .csv format")
      )      
      dttf <- read.csv(infile$datapath, header = TRUE)
      return(dttf)
    })
  
  
  #display data as a table
  draw_summary_table <- function(inputdf){
    inputdf <- mutate(inputdf, growth_stage = as.factor(growth_stage), replicate = as.factor(replicate))
    final <- data.frame("Column_Name" = colnames(inputdf),
                        "Type" = sapply(inputdf[2,], FUN = function(x) class(x)),
                        "Distinct_Values" = character(3))
    final$Distinct_Values[1] <- paste(unique(inputdf$sample),collapse=', ')
    final$Distinct_Values[2] <- paste(unique(inputdf$growth_stage),collapse=', ')
    final$Distinct_Values[3] <- paste(unique(inputdf$replicate),collapse=', ')
    
    #I(list(unique(inputdf$sample),unique(inputdf$growth_stage), unique(inputdf$replicate))))
    return(final)
  }
  
  # plotting functions
  make_barplot_replicate <- function(inputdf){
    plot <- ggplot(inputdf,aes(x = replicate)) + geom_bar(fill = "red", color = "black", width = 0.7, size = 1)+
      labs(title = "Replicates", x = "Replicate", y = "Count") + theme_dark() + theme(legend.position = "none")
    return(plot)
  }
  make_barplot_growth_stage <- function(inputdf){
    plot <- ggplot(inputdf,aes(x = growth_stage)) + geom_bar(fill = "blue", color = "black", width = 0.7, size=1)+
      labs(title = "Growth Stages", x = "Growth Stage", y = "Count") + theme_dark() + theme(legend.position = "none")
    return(plot)
  }
  # making outputs for first tab
    output$summary_table <- renderTable(draw_summary_table(load_summary_table_data()))#output table
    output$table <- renderDT(datatable(load_summary_table_data(),options = list(ordering = TRUE)))
    output$replicate_plot <- renderPlot(make_barplot_replicate(load_summary_table_data()))
    output$growth_stage_plot <- renderPlot(make_barplot_growth_stage(load_summary_table_data()))
    
    # all below is for counts tab
    load_counts_data <- reactive({ #load data reactive, with validate statement to prevent warnings, return dataframe
      infile <- input$counts_file
      validate(
        need(input$counts_file != "", "Please select a data set in .csv format")
      )      
      dttf <- read.csv(infile$datapath, header = TRUE,row.names = 1)
      return(dttf)
    }) 
    
    make_count_summary_table <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf
      num_sample <-ncol(dff)
      num_genes <- nrow(dff)
      
      filtered_df <- as_tibble(dff) %>%
        mutate(row_var = rowVars(as.matrix(dff)),
               percentile = percent_rank(row_var)*100) %>%
        filter(rowSums(. != 0) >= slidernonzero) %>% # need to add 1 because rownames are nonzero
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
    
    #scatterplots
    
    make_scatter_med_var <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf
      filtered_df <- as_tibble(dff) %>%
        mutate(row_var = rowVars(as.matrix(dff)),
               percentile = percent_rank(row_var)*100,
               nonzero = rowSums(. != 0),
               row_med = rowMedians(as.matrix(dff))) 
      filtered_df$Threshold <- ifelse(filtered_df$percentile >= slidervar & filtered_df$nonzero >= slidernonzero, 'TRUE', 'FALSE')
      plot <- ggplot(filtered_df, aes(x = log(row_med), y = log(row_var), color = Threshold)) + theme_dark() +
        scale_color_manual(name ="Within Threshold", values = c('FALSE' = "red", 'TRUE' = "green"))+
        geom_point() +
        labs(x = "log10(Median Count)", y = "log10(Count Variance)", title = "Median Count vs. Count Variance")
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
        scale_color_manual(name ="Within Threshold", values = c('FALSE' = "red", 'TRUE' = "green"))+
        geom_point() +
        labs(x = "log10(Median Count)", y = "# of Zero Values/row", title = "Median Count vs. # of Zero Values/row")
      return(plot)
    }
    
    # make heatmap
    make_heatmap <- function(inputdf,slidervar,slidernonzero){
      dff <- inputdf # for some reason need to assign variable to reactive object or it does not work as well
      filtered_df <- as_tibble(dff, rownames = NA) %>% #filtering while preserving rownames
        rownames_to_column() %>%
        mutate(row_var = rowVars(as.matrix(dff)), #adding row which is
               percentile = percent_rank(row_var)*100) %>%
        filter(rowSums(. != 0) -3 >= slidernonzero) %>% # need to subtract 3 because I added 3 nonzero rows
        filter(percentile >= slidervar) %>%
        select(-row_var,-percentile)

      input <- as.data.frame(filtered_df[,-1])
      rownames(input) <- filtered_df$rowname

      heat_data <- as.matrix(input)
      par(mar=c(5.1,4.1,4.1,4))
      heatmap.2(heat_data, scale = "column", trace = 'none', adjRow = c(.1,0), margins = c(10,10))
    }
    
    plot_pca <- function(inputdf,input1_pc,input2_pc, input3_fill) {
      input <- inputdf # for some reason need to assign variable to reactive object or it does not work as well
      pca <- prcomp(t(input)) 
      plotting_data <- data.frame('Sample_Name' = colnames(input),
                                  'Growth_Stage' = c(rep("early",3),rep("middle",3),rep("late",3)),
                                  'Replicate' = as.character(rep(c(1,2,3),3)),
                                  pca$x)
      pct_var1 <- round(pca$sdev[as.numeric(str_extract_all(input1_pc, "\\d")[[1]])]^2/sum(pca$sdev^2)*100, digits =1)
      pct_var2 <- round(pca$sdev[as.numeric(str_extract_all(input2_pc, "\\d")[[1]])]^2/sum(pca$sdev^2)*100, digits =1)
      

      
      plot <- ggplot(plotting_data,aes_string(x = input1_pc, y = input2_pc, color = input3_fill)) + 
        geom_point() + labs(x = paste(input1_pc,', % Variance Explained:',pct_var1) , y = paste(input2_pc,', % Variance Explained:',pct_var2))
      return(plot)
    }
    # volcano_plot <-
    #   function(dataf, x_name, y_name, slider, color1, color2) { #argument inputs
    #     dataf <- na.omit(dataf) #remove nas
    #     xval <- as.numeric(unlist(dataf[x_name])) #make vector of xvals
    #     yval <- as.numeric(unlist(dataf[y_name])) #make vector of yvals
    #     dataf$pointcolor <- ifelse(dataf$padj>1*10^slider,'FALSE','TRUE') #make new column to store color values
    #     p1 <- ggplot(dataf,aes(x=xval,y=-log10(yval), color = pointcolor)) + theme_bw()+ #make plot
    #       geom_point()+
    #       scale_color_manual(name =paste0("padj < 1 x 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2)) +
    #       theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
    #             legend.position = "bottom",legend.text = element_text(size = 25))+
    #       labs(x = x_name, y= paste0("-log10(",y_name,")"))
    #     return(p1)
    #   }
    
    
    
    
  output$counts_summary_table <- renderTable(make_count_summary_table(load_counts_data(),input$var_slide,input$nonzero_slide))
  output$med_v_var <- renderPlot(make_scatter_med_var(load_counts_data(),input$var_slide,input$nonzero_slide))
  output$med_v_nzeros <- renderPlot(make_scatter_med_zeros(load_counts_data(),input$var_slide,input$nonzero_slide))
  output$heatmap <- renderPlot(
                               make_heatmap(load_counts_data(),input$var_slide,input$nonzero_slide), height = 800)
  output$pca_scatter <- renderPlot(plot_pca(load_counts_data(),input$pc_select1,input$pc_select2,input$pca_fill))
    
  
  #DE Panel is below
  
  load_DE_summary_table_data <- reactive({ #load data reactive, with validate statement to prevent warnings, return dataframe
    infile <- input$DE_file
    validate(
      need(input$DE_file != "", "Please select a data set in .csv format")
    )      
    dttf <- read.csv(infile$datapath, header = TRUE)
    return(dttf)
  })
  
  
  #for volcano plot
  volcano_plot <-
    function(dataf,comparison, x_name, y_name, slider, color1, color2) { #argument inputs
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
      
      dataf <- na.omit(dataf) #remove nas
      xval <- as.numeric(unlist(dataf[paste0(comp,'.',x_final)])) #make vector of xvals
      yval <- as.numeric(unlist(dataf[paste0(comp,'.',y_final)])) #make vector of yvals
      dataf$pointcolor <- ifelse(dataf[paste0(comp,'.','FDR')]>1*10^slider,'FALSE','TRUE') #make new column to store color values
      p1 <- ggplot(dataf,aes(x=xval,y=-log10(yval), color = pointcolor)) + theme_dark()+ #make plot
        geom_point()+
        scale_color_manual(name =paste0("FDR adjusted p-value < 1 x 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2)) +
        theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.position = "bottom")+
        labs(x = x_name, y= paste0("-log10(",y_name,")"))
      return(p1)
    }
  
  
  
  output$de_summary_table <- renderDT(datatable(load_DE_summary_table_data(),options = list(ordering = TRUE))%>%
                                        formatRound(c(2:16), 3))
  output$volcano <- renderPlot(volcano_plot(load_DE_summary_table_data(),input$var1, input$button1, input$button2, input$slide, input$col1, input$col2),
                               width = 800, height = 600)
  
  #individual gene panel
  
# gene name is input$searchme
# plot type is input$field2
# category is input$field1


  
  make_final_plot <- function(inputdf, gene, category, plot_type){
    for_plot <- data.frame('Sample_Name' = suppressWarnings(colnames(inputdf)),
                           'Growth_Stage' = c(rep("early",3),rep("middle",3),rep("late",3)),
                           'Replicate' = as.character(rep(c(1,2,3),3)),
                           t(inputdf[gene,]))
    
    make_barplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_bar(stat = 'identity', fill = 'blue') + labs(title = gene, y = 'Count') + theme_dark()
      return(p)
    }
    
    
    make_boxplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_boxplot(fill = "blue") + labs(title = gene, y = 'Count') + theme_dark()
      return(p)
    }
    
    
    make_violinplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_violin(fill = "blue") + labs(title = gene, y = 'Count') + theme_dark()
      return(p)
    }
    
    
    make_beeswarmplot <- function(inputdf, category, gene){
      p <- ggplot(inputdf, aes_string(x = category, y = gene)) +
        geom_beeswarm(fill = "blue") + labs(title = gene, y = 'Count') + theme_dark()
      return(p)
    }
    
    plot <- if(plot_type == 'bar'){
      make_barplot(for_plot, category, gene)
    } else if (plot_type == 'boxplot'){
      make_boxplot(for_plot, category, gene)
    } else if (plot_type == 'violin'){
      make_violinplot(for_plot, category, gene)
    } else {make_beeswarmplot(for_plot, category, gene)}
    
    
    return(plot)
  }
  
  
  
  output$individual_gene_plot <- suppressWarnings(renderPlot({
    input$plot_button
    isolate(make_final_plot(load_counts_data(),input$searchme,input$field1,input$field2))
    
  }))
    
  
#     eventReactive(make_final_plot(load_counts_data(),input$searchme,input$field1,input$field2)){runif(input$plot_button)}
#   
#   observeEvent(input$plot_button)
# output$individual_gene_plot <- renderPlot(observeEvent(input$plot_button, {
#   make_final_plot(load_counts_data(),input$searchme,input$field1,input$field2)
# }))
  
#   renderPlot( reactive({ if (input$plot_button >0){
#   make_final_plot(load_counts_data(),input$searchme,input$field1,input$field2)
# }
#   
# }))
#   
# write function that takes gene input, category, plot type and then plots it 
  
  
  
}

# Create Shiny app object
shinyApp(ui = ui, server = server)



}