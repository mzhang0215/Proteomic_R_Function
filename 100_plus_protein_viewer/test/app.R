# ## Only run examples in interactive R sessions
# if (interactive()) {
#   
#   library("shiny")
#   library("shinyWidgets")
#   
#   
#   # simple use
#   
#   ui <- fluidPage(
#     multiInput(
#       inputId = "id", label = "Fruits :",
#       choices = c("Banana", "Blueberry", "Cherry",
#                   "Coconut", "Grapefruit", "Kiwi",
#                   "Lemon", "Lime", "Mango", "Orange",
#                   "Papaya"),
#       selected = "Banana", width = "350px"
#     ),
#     verbatimTextOutput(outputId = "res")
#   )
#   
#   server <- function(input, output, session) {
#     output$res <- renderPrint({
#       input$id
#     })
#   }
#   
#   shinyApp(ui = ui, server = server)
#   
#   
#   # with options
#   
#   ui <- fluidPage(
#     multiInput(
#       inputId = "id", label = "Fruits :",
#       choices = c("Banana", "Blueberry", "Cherry",
#                   "Coconut", "Grapefruit", "Kiwi",
#                   "Lemon", "Lime", "Mango", "Orange",
#                   "Papaya"),
#       selected = "Banana", width = "400px",
#       options = list(
#         enable_search = FALSE,
#         non_selected_header = "Choose between:",
#         selected_header = "You have selected:"
#       )
#     ),
#     verbatimTextOutput(outputId = "res")
#   )
#   
#   server <- function(input, output, session) {
#     output$res <- renderPrint({
#       input$id
#     })
#   }
#   
#   shinyApp(ui = ui, server = server)
#   
# }

library(shiny)
library(dplyr)
data(mtcars)

ui <- fluidPage(
  titlePanel("MTCARS"),
  selectInput("Columns","Columns",
              names(mtcars), multiple = TRUE),
  verbatimTextOutput("dfStr")
)

server <- function(input, output) {
  Dataframe2 <- reactive({
    mtcars[,input$Columns] 
  })
  output$dfStr <- renderPrint({
    str(Dataframe2())
  })
}

shinyApp(ui, server)