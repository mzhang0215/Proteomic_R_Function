#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinyWidgets)
library(plotly)
library(dplyr)

load(file = "protein_display_dataset_new2.Rdata")

#protein_display_dataset <- as.data.frame(fread("protein_display_dataset.tsv", header =TRUE),
#                                         stringsAsFactors=F)

button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;
/* Change the text size to 15 pixels. */
font-size: 15px;
}"


# Define UI for application that shows protein in volcano plot
shinyUI(fluidPage(
  
  #Navbar structure for UI
  navbarPage("100 Plus proteomics", theme = shinytheme("lumen"),
             tabPanel("Protein Browser", icon = icon("brain"),
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(width = 3, 
                                     # tags$head(
                                     #   tags$style(HTML("
                                     #                   .my_table .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                                     #                   vertical-align: center;
                                     #                   border: 2px solid black; 
                                     #                   }
                                     #                   "))
                                     #   ),
                                     
                                     titlePanel("Quick Search"),
                                     #shinythemes::themeSelector(),
                                     fluidRow(column(12,
                                                     # Select which Gender(s) to plot
                                                     selectizeInput(
                                                       inputId = "protein_symbol", label = "",
                                                       selected = "MAPT",
                                                       choices = protein_display_dataset$protein,
                                                       width = "800px"
                                                     )
                                                     )),
                                     hr(),
                                     fluidRow(column(12,
                                                     h3("Full name"),
                                                     tags$style("#full_name {font-size:18px;
                                                                color:gray;}"),
                                                     textOutput("full_name")
                                                     )),

                                     fluidRow(column(12,
                                                     h3("Cell-type marker"),
                                                     tags$style("#celltype_marker {font-size:18px;
                                                                color:gray;}"),
                                                     textOutput("celltype_marker")
                                                     )),
                                     
                                     hr(),
                                     fluidRow(column(12,
                                                     h3("Pearson correlation with abundance"),
                                                     tableOutput('correlation_table'),
                                                     h5("Pearson correlation between protein abundance and Age, Amyloid stage, and Braak stage in AD, ND, and Centenarian groups.")
                                                     )),
                                     hr(),
                                     fluidRow(column(12,
                                                     h3("Abundance test in paired groups"),
                                                     tableOutput('pairTest_table'),
                                                     h5("T-test for abundances between the paired groups, i.e., AD vs. ND, AD vs. CEN, and ND vs. CEN. The upper right shows the p-values, and the lower left shows the t-statistics.")
                                                     ))
                      ),
                      mainPanel(width = 9,
                                fluidRow(column(12,
                                                plotlyOutput(outputId = "volcano_plot", width = "95%", height = "600px")
                                  )),
                                hr(),
                                fluidRow(column(12,
                                                plotlyOutput(outputId = "scatter_plot", width = "95%", height = "400px")
                                  )),
                                hr(),
                                fluidRow(column(12,
                                                plotlyOutput(outputId = "single_line_plot", width = "95%", height = "400px")
                                ))
                                )
                      )
                      ),
             tabPanel("Abundance vs. Braak stage & Age", icon = icon("random"),
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(width = 3, 
                                     titlePanel("Protein[s] selection"),
                                     fluidRow(column(12,
                                                     selectizeInput(inputId = "proteins",
                                                                    label = "",
                                                                    selected = c("MAPT", "APP", "SIRT2"),
                                                                    choices = protein_display_dataset$protein,
                                                                    width = "800px", 
                                                                    multiple = TRUE)
                                     ))
                        ),
                        mainPanel(width = 9,
                                  fluidRow(
                                    column(12,
                                           plotlyOutput(outputId = "line_plot", width = "95%", height = "800px")
                                  ))
                        )
                      )
             ),
             tabPanel("Pairwise correlations", icon = icon("fa-sharp fa-solid fa-chart-line"),
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     titlePanel("Protein selection"),
                                     fluidRow(column(12,
                                                     selectizeInput(
                                                       inputId = "proteins2", label = "",
                                                       selected = c("MAPT", "APP", "SIRT2"),
                                                       choices = protein_display_dataset$protein,
                                                       width = "800px",
                                                       multiple = TRUE
                                                     )
                                     )),
                                     hr(),
                                     fluidRow(column(12,
                                                     h3("Correlation method"),
                                                     selectizeInput(inputId = "cormethod",
                                                                    label = "",
                                                                    selected = "pearson",
                                                                    choices = c("pearson", "spearman"),
                                                                    width = "800px"),
                                                     h5("X-cross indicates the correlation is not significant.")
                                     ))
                        ),
                        mainPanel(width = 9,
                                  fluidRow(
                                    column(12,
                                           plotlyOutput(outputId = "protein_corr_plot",
                                                        width = "1000px",
                                                        height = "800px")
                                    ))
                        )
                      )
             ),
             tabPanel("Abundance vs. APOE-e4", icon = icon("fa-solid fa-dna"),
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     titlePanel("Protein selection"),
                                     fluidRow(column(12,
                                            selectizeInput(
                                              inputId = "protein_APOE", label = "",
                                              selected = "MAPT",
                                              choices = protein_display_dataset$protein,
                                              width = "800px"
                                              )
                                     )),
                                     hr(),
                                     fluidRow(column(12,
                                                     h3("Test between APOE-e4+ and APOE-e4- in AD group"),
                                                     tableOutput('APOE_t_table'),
                                                     h5("T-test between APOE-e4 positive and APOE-e4 negative samples in AD group.")
                                     ))
                        ),
                        mainPanel(width = 9,
                                  fluidRow(
                                    column(12,
                                           plotlyOutput(outputId = "protein_APOE_plot", width = "95%", height = "90%")
                                  ))
                        )
                      )
             ),
             tabPanel("Sample Characteristic", icon = icon("chart-bar"),
                      tags$style(button_color_css),
                      # Sidebar layout with a input and output definitions
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     titlePanel("Sample characteristic"),
                                     selectInput(inputId = "variable1",
                                                 label = "Select variable1:",
                                                 choices = c("Age", "Braak stage", "Amyloid stage", "Sex", "APOE allele"),
                                                 selected = "Age",
                                                 width = "220px"
                                     ),
                                     selectInput(inputId = "variable2",
                                                 label = "Select variable2:",
                                                 choices = c("None", "Age", "Braak stage", "Amyloid stage", "Sex", "APOE allele"),
                                                 selected = "None",
                                                 width = "220px"
                                     )
                                     
                          
                        ),
                        mainPanel(width = 6,
                                  plotOutput(outputId = "distribution_plot"))
                      )
                      ),
             navbarMenu("More", icon = icon("info-circle"),
                        tabPanel("100 Plus proteomics", fluid = TRUE),
                        tabPanel("100 Plus study", fluid = TRUE)
                        )
  )
))
