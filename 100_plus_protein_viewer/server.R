#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(plotly)
#library(ggthemr)
library(ggcorrplot)

load(file = "protein_display_dataset_new2.Rdata")
load(file = "metadata_imp2.Rdata")
load(file = "protein_intensity_combat.Rdata")

load(file = "exps_at_stage_ADND.Rdata")
load(file = "exps_at_stage_CHC.Rdata")

load(file = "exps_at_Astage_ADND.Rdata")
load(file = "exps_at_Astage_CHC.Rdata")

load(file = "exps_at_age_NDCHC.Rdata")
load(file = "exps_at_age_AD.Rdata")
load(file = "group_wise_test.Rdata")

Astrocytes = "#FEB24C"
Endothelial = "#88419D"
Microglia = "#9ECAE1"
ExNeurons = "#ADDD8E"
InNeurons = "#339900"
OPC = "#A50F15"
Oligodendrocytes = "#DD3497"
Unknown = "#008080"
No_unique_enrichment = "#969696"
NS_DEPs = "gray88"

color_AD <- "#FF3300"
color_ND <- "blue"
color_CHC <- "#FFCC00"
color_ADND <- "purple"
color_NDCHC <- "green"

point_size = 15
tb_margin = 100
height0 = 900
width0 = 900

NDCHC_age0 <- min(as.numeric(colnames(exps_at_age_NDCHC)))
NDCHC_age1 <- max(as.numeric(colnames(exps_at_age_NDCHC)))

AD_age0 <- min(as.numeric(colnames(exps_at_age_AD)))
AD_age1 <- max(as.numeric(colnames(exps_at_age_AD)))

age_min <- min(NDCHC_age0, AD_age0)
age_max <- max(NDCHC_age1, AD_age1)

AD_index <- which(metadata$GROUP == "AD")
ND_index <- which(metadata$GROUP == "ND")
CHC_index <- which(metadata$GROUP == "CHC")

tmp <- plotly:::Schema
tmp$traces$scatter$attributes$line$dash$values <- rep(tmp$traces$scatter$attributes$line$dash$values, 10)
assignInNamespace("Schema", tmp, ns = "plotly")

alpha_value <- c(rep(1, sum(protein_display_dataset$celltype != "")), 
                 rep(0.3, sum(protein_display_dataset$celltype == "")))

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  cdata <- session$clientData
  
  output$full_name <- renderText({
    protein_display_dataset[which(protein_display_dataset$protein == input$protein_symbol),]$NAME
  })
  
  output$celltype_marker <- renderText({
    if(length(which(protein_display_dataset$protein == input$protein_symbol)) == 0)
      return("")
    celltype <- protein_display_dataset[which(protein_display_dataset$protein == input$protein_symbol),]$celltype
    if(celltype == ""){
      return(paste(input$protein_symbol, "is not a cell-type marker protein."))
    }
    if(celltype == "Astrocytes"){
      return(paste(input$protein_symbol, "is a marker of Astrocytes."))
    }
    if(celltype == "Endothelial"){
      return(paste(input$protein_symbol, "is a marker of Endothelial Cell."))
    }
    if(celltype == "Microglia"){
      return(paste(input$protein_symbol, "is a marker of Microglia."))
    }
    if(celltype == "ExNeurons"){
      return(paste(input$protein_symbol, "is a marker of ExNeurons."))
    }
    if(celltype == "InNeurons"){
      return(paste(input$protein_symbol, "is a marker of InNeurons."))
    }
    if(celltype == "OPC"){
      return(paste(input$protein_symbol, "is a marker of OPC."))
    }
    if(celltype == "Oligodendrocytes"){
      return(paste(input$protein_symbol, "is a marker of Oligodendrocytes."))
    }
    if(celltype == "Unknown"){
      return(paste(input$protein_symbol, "is a marker of Unknown-type cell."))
    }
  })
  
  output$correlation_table <- renderTable({
    if(length(which(protein_display_dataset$protein == input$protein_symbol)) != 0){
      data_matrix <- subset(protein_display_dataset, protein == input$protein_symbol, 
             select=c("braak_rho", "braak_p", "amyloid_rho", "amyloid_p", "age_rho", "age_p"))
      data_matrix <- data.frame(rho = t(data_matrix[c(1, 3, 5)]), p=t(data_matrix[c(2, 4, 6)]))
      data_matrix <- cbind(c("Braak stage", "Amlyoid stage", "Age"), data_matrix)
      colnames(data_matrix) <- c("", "coefficient", "p-value")
      return(data_matrix)
    }
  },
  bordered = TRUE,
  spacing = 'm',
  width = "100%"
  )
  
  output$pairTest_table <- renderTable({
    if(input$protein_symbol %in% names(group_wise_test)){
      data_matrix <- group_wise_test[[input$protein_symbol]]
      return(data_matrix)
    }
  },
  bordered = TRUE,
  spacing = 'm',
  width = "100%"
  )
  
  output$APOE_t_table <- renderTable({
    if(length(which(protein_display_dataset$protein == input$protein_APOE)) != 0){
      data_matrix <- subset(protein_display_dataset, protein == input$protein_APOE, 
                            select=c("APOE_t", "APOE_FDR"))
      data_matrix <- cbind(c(input$protein_APOE), data_matrix)
      colnames(data_matrix) <- c("", "t statistic", "p-value (after FDR)")
      return(data_matrix)
    }
  },
  bordered = TRUE,
  spacing = 'm',
  width = "100%"
  )

  
  x_lfc <- list(
    title = "Flod change (log2)",
    dtick = 2,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  x_slope <- list(
    title = "Linear slope",
    dtick = 0.25,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  x_slope_age <- list(
    title = "Linear slope",
    dtick = 0.02,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  y_FDR <- list(
    title = "FDR (-log10)",
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  y_2 <- list(
    title = "",
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  
  f <- list( size = 15, color = "black")
  title_a <- list(
    text = "Braak stage related protein, ANOVA test",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = -4,
    y = 1,
    yshift=25,
    showarrow = FALSE
  )
  
  title_b <- list(
    text = "Braak stage related protein, LM+F test",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = -0.6,
    y = 1,
    yshift=25,
    showarrow = FALSE
  )
  
  title_c <- list(
    text = "Age related protein, LM+F test",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = -0.07,
    y = 1,
    yshift=25,
    showarrow = FALSE
  )
  
  output$volcano_plot <- renderPlotly({
    if (length(which(protein_display_dataset$protein == input$protein_symbol)) == 0){
      subplot(
        plot_ly(protein_display_dataset, x = ~lfc, y = ~(-log10(ANOVA_adj_p)), text = ~tips, key = ~protein, color = ~marker, 
                colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                marker = list(opacity = alpha_value),
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "FDR (-log10): %{y:.2f}<br>",
                  "%{xaxis.title.text}: %{x:.2f}",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~marker, showlegend = TRUE) %>%
          add_markers(size = I(point_size)) %>%
          layout(xaxis = x_lfc, yaxis = y_FDR, annotations = title_a, legend = list(x = 0, y = -0.2, orientation = 'h')
                 ),
        
        plot_ly(protein_display_dataset, x = ~slope, y = ~(-log10(lm_adj_p)), text = ~tips, key = ~protein, color = ~marker, 
                colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                marker = list(opacity = alpha_value),
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "FDR (-log10): %{y:.2f}<br>",
                  "%{xaxis.title.text}: %{x:.2f}",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~marker, showlegend = FALSE) %>%
          add_markers(size = I(point_size)) %>%
          layout(xaxis = x_slope, yaxis = y_2, annotations = title_b),
        
        plot_ly(protein_display_dataset, x = ~effect_size, y = ~(-log10(age_adj_p)), text = ~tips, key = ~protein, color = ~marker,  
                colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                marker = list(opacity = alpha_value),
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "FDR (-log10): %{y:.2f}<br>",
                  "%{xaxis.title.text}: %{x:.2f}",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~marker, showlegend = FALSE) %>%
          add_markers(size = I(point_size)) %>%
          layout(xaxis = x_slope_age, yaxis = y_2, annotations = title_c),
        
        titleX = TRUE,
        titleY = TRUE
      )

    }
    else{
      protein <- protein_display_dataset[which(protein_display_dataset$protein == input$protein_symbol), ]
      a_lfc <- list(
        x = protein$lfc,
        y = -log10(protein$ANOVA_adj_p),
        text = protein$protein,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 2,
        ax = 20,
        ay = -30
      )
      a_lm <- list(
        x = protein$slope,
        y = -log10(protein$lm_adj_p),
        text = protein$protein,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 2,
        ax = 20,
        ay = -30
      )
      a_age <- list(
        x = protein$effect_size,
        y = -log10(protein$age_adj_p),
        text = protein$protein,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 2,
        ax = 20,
        ay = -30
      )
      
      subplot(plot_ly(protein_display_dataset, x = ~lfc, y = ~(-log10(ANOVA_adj_p)), text = ~tips, key = ~protein, color = ~marker,  
                      colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                      marker = list(opacity = alpha_value),
                      hovertemplate = paste(
                        "<b>&nbsp;%{text}</b>",
                        "FDR (-log10): %{y:.2f}<br>",
                        "%{xaxis.title.text}: %{x:.2f}",
                        "<extra></extra>"
                      ), height=600,
                      legendgroup = ~marker, showlegend = TRUE) %>%
                add_markers(size = I(point_size)) %>%
                layout(xaxis = x_lfc, yaxis = y_FDR, annotations = list(a_lfc, title_a), legend = list(x = 0, y = -0.2, orientation = 'h')
                       ),
              
              plot_ly(protein_display_dataset, x = ~slope, y = ~(-log10(lm_adj_p)), text = ~tips, key = ~protein, color = ~marker,
                      colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                      marker = list(opacity = alpha_value),
                      hovertemplate = paste(
                        "<b>&nbsp;%{text}</b>",
                        "FDR (-log10): %{y:.2f}<br>",
                        "%{xaxis.title.text}: %{x:.2f}",
                        "<extra></extra>"
                      ), height=600,
                      legendgroup = ~marker, showlegend = FALSE) %>%
                add_markers(size = I(point_size)) %>%
                layout(xaxis = x_slope, yaxis = y_2, annotations = list(a_lm, title_b)),
              
              plot_ly(protein_display_dataset, x = ~effect_size, y = ~(-log10(age_adj_p)), text = ~tips, key = ~protein, color = ~marker,   
                      colors = c(Astrocytes, Endothelial, ExNeurons, InNeurons, Microglia, OPC, Oligodendrocytes, Unknown, No_unique_enrichment),
                      marker = list(opacity = alpha_value),
                      hovertemplate = paste(
                        "<b>&nbsp;%{text}</b>",
                        "FDR (-log10): %{y:.2f}<br>",
                        "%{xaxis.title.text}: %{x:.2f}",
                        "<extra></extra>"
                      ), height=600,
                      legendgroup = ~marker, showlegend = FALSE) %>%
                add_markers(size = I(point_size)) %>%
                layout(xaxis = x_slope_age, yaxis = y_2, annotations = list(a_age, title_c)),
              
              titleX = TRUE,
              titleY = TRUE
              )
    }
  })
  
  x_braak <- list(
    title = "Braak stage",
    dtick = 1,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    range = list("0/I", "VI"),
    linewidth = 0.5
  )
  x_Abeta <- list(
    title = "NIA Amyloid stage",
    dtick = 1,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  x_group <- list(
    title = "Group",
    dtick = 1,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  x_age <- list(
    title = "Age [y]",
    dtick = 10,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )
  y <- list(
    title = "Protein abundance (log2)",
    zeroline = FALSE,
    linecolor = toRGB("black"),
    linewidth = 0.5
  )

  x_age2 <- list(
    title = "Age [y]",
    dtick = 5,
    zeroline = FALSE,
    linecolor = toRGB("black"),
    range = list(age_min, age_max),
    linewidth = 0.5
  )
  
  
  pre_protein = "MAPT"
  
  output$scatter_plot <- renderPlotly({
    protein = ""
    if (length(which(protein_display_dataset$protein == input$protein_symbol)) == 0){
      protein <- pre_protein
    }
    else{
      protein <- input$protein_symbol
      pre_protein <<- input$protein_symbol
    }
    
    plot_matrix <- data.frame(braak = metadata$Braak_text,
                              Abeta = metadata$amyloid,
                              group = metadata$group_n,
                              age = as.numeric(metadata$age),
                              diagnosis = as.character(metadata$diagnosis),
                              intensity = as.numeric(protein_intensity_combat[which(rownames(protein_intensity_combat) == protein),]))
    
    subplot(
      plot_ly(plot_matrix, x = ~Abeta, y = ~intensity, color = ~group, text = ~diagnosis, type = "box", boxpoints = "all", pointpos = 0,
              colors = c(color_AD, color_ND, color_CHC), source = "B",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ),
              legendgroup = ~group, showlegend = TRUE) %>%
        layout(xaxis = x_Abeta, yaxis = y, boxmode = "group", 
               legend = list(x = 0.3, y = -0.2, orientation = 'h')),
      
      plot_ly(plot_matrix, x = ~braak, y = ~intensity, color = ~group, text = ~diagnosis, type = "box", boxpoints = "all", pointpos = 0,
              colors = c(color_AD, color_ND, color_CHC), source = "B",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ),
              legendgroup = ~group, showlegend = FALSE) %>%
        layout(xaxis = x_braak, yaxis = y_2, boxmode = "group"
               ),
      
      plot_ly(plot_matrix, x = ~age, y = ~intensity, color = ~group, text = ~diagnosis,
              colors = c(color_AD, color_ND, color_CHC), source = "B",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ),
              legendgroup = ~group, showlegend = FALSE) %>%
        add_markers(size = I(point_size))%>%
        layout(xaxis = x_age, yaxis = y_2),
      
      titleX = TRUE,
      titleY = TRUE
    )
    
  })
  
  output$single_line_plot <- renderPlotly({
    protein = ""
    if (!(input$protein_symbol %in% protein_display_dataset$protein)){
      protein <- pre_protein
    }
    else{
      protein <- input$protein_symbol
      pre_protein <<- input$protein_symbol
    }
    
    exps_at_Astage_ADND_sub <- exps_at_Astage_ADND[which(rownames(exps_at_Astage_ADND) == protein), ]
    amyloid  <- colnames(exps_at_Astage_ADND_sub)
    amyloid <- rep(amyloid, nrow(exps_at_Astage_ADND_sub))
    intens <- as.vector(t(exps_at_Astage_ADND_sub))
    plot_df_amyloid_ADND <- data.frame(intens = intens,
                                       amyloid = amyloid)
    plot_df_amyloid_ADND$group <- "AD and ND groups"
    exps_at_Astage_CHC_sub <- exps_at_Astage_CHC[which(rownames(exps_at_Astage_CHC) == protein), ]
    amyloid  <- colnames(exps_at_Astage_CHC_sub)
    amyloid <- rep(amyloid, nrow(exps_at_Astage_CHC_sub))
    intens <- as.vector(t(exps_at_Astage_CHC_sub))
    plot_df_amyloid_CHC <- data.frame(intens = intens,
                                      amyloid = amyloid)
    plot_df_amyloid_CHC$group <- "Centenarian group"
    plot_df_amyloid <- rbind(plot_df_amyloid_ADND, plot_df_amyloid_CHC)
    plot_df_amyloid$amyloid <- factor(plot_df_amyloid$amyloid, levels = c("0", "1", "2", "3"))
    
    exps_at_stage_ADND_sub <- exps_at_stage_ADND[which(rownames(exps_at_stage_ADND) == protein), ]
    braak  <- colnames(exps_at_stage_ADND_sub)
    braak <- rep(braak, nrow(exps_at_stage_ADND_sub))
    intens <- as.vector(t(exps_at_stage_ADND_sub))
    plot_df_Braak_ADND <- data.frame(intens = intens,
                                     braak = braak)
    plot_df_Braak_ADND$group <- "AD and ND groups"
    exps_at_stage_CHC_sub <- exps_at_stage_CHC[which(rownames(exps_at_stage_CHC) == protein), ]
    braak  <- colnames(exps_at_stage_CHC_sub)
    braak <- rep(braak, nrow(exps_at_stage_CHC_sub))
    intens <- as.vector(t(exps_at_stage_CHC_sub))
    plot_df_Braak_CHC <- data.frame(intens = intens,
                                    braak = braak)
    plot_df_Braak_CHC$group <- "Centenarian group"
    plot_df_Braak <- rbind(plot_df_Braak_ADND, plot_df_Braak_CHC)
    plot_df_Braak$braak <- factor(plot_df_Braak$braak, levels = c("0/I", "II", "III", "IV", "V", "VI"))
    
    exps_at_age_NDCHC_sub <- exps_at_age_NDCHC[which(rownames(exps_at_age_NDCHC) == protein), ]
    age  <- as.numeric(colnames(exps_at_age_NDCHC_sub))
    age <- rep(age, nrow(exps_at_age_NDCHC_sub))
    intens <- as.vector(t(exps_at_age_NDCHC_sub))
    plot_df_age_NDCHC <- data.frame(intens = intens,
                                    age = age)
    plot_df_age_NDCHC$group <- "ND group"
    plot_df_age_NDCHC[which(plot_df_age_NDCHC$age >= 100), ]$group <- "Centenarian group"
    exps_at_age_AD_sub <- exps_at_age_AD[which(rownames(exps_at_age_AD) == protein), ]
    age  <- as.numeric(colnames(exps_at_age_AD_sub))
    age <- rep(age, nrow(exps_at_age_AD_sub))
    intens <- as.vector(t(exps_at_age_AD_sub))
    plot_df_age_AD <- data.frame(intens = intens,
                                 age = age)
    plot_df_age_AD$group <- "AD group"
    plot_df_age <- rbind(plot_df_age_NDCHC, plot_df_age_AD)
    
    y_min <- min(min(plot_df_amyloid$intens), min(plot_df_Braak$intens), min(plot_df_age$intens)) - 0.05
    y_max <- max(max(plot_df_amyloid$intens), max(plot_df_Braak$intens), max(plot_df_age$intens)) + 0.05
    
    y_ratio <- list(
      title = "LFQ intensity ratio",
      range = list(y_min, y_max), 
      zeroline = FALSE,
      linecolor = toRGB("black"),
      linewidth = 0.5
    )
    
    y_ratio2 <- list(
      title = "",
      range = list(y_min, y_max), 
      zeroline = FALSE,
      linecolor = toRGB("black"),
      linewidth = 0.5
    )
    
    #print(plot_df_age)
    subplot(
      plot_ly(plot_df_amyloid, x = ~amyloid, y = ~intens, type = 'scatter', mode = 'lines+markers', color = ~group, colors = c(color_ADND, color_CHC), 
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b>",
                "(%{x}, %{y:.2f})",
                "<extra></extra>"
              ), 
              legendgroup = ~group, showlegend = TRUE) %>%
        layout(xaxis = x_Abeta, yaxis = y_ratio, shapes = list(hline())#, legend = list(x = 0, y = -0.2, orientation = 'h')
        ),
      
      plot_ly(plot_df_Braak, x = ~braak, y = ~intens, type = 'scatter', mode = 'lines+markers',  color = ~group, colors = c(color_ADND, color_CHC), 
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b>",
                "(%{x}, %{y:.2f})",
                "<extra></extra>"
              ), 
              legendgroup = ~group, showlegend = FALSE) %>%
        layout(xaxis = x_braak, yaxis = y_ratio2, shapes = list(hline(x1 = "VI"))#, legend = list(x = 0, y = -0.1, orientation = 'h')
        ),
      
      plot_ly(subset(plot_df_age, age <100), x = ~age, y = ~intens, type = 'scatter', mode = 'lines+markers', color = ~group, colors = c(color_AD, color_CHC, color_ND), 
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b>",
                "(%{x}, %{y:.2f})",
                "<extra></extra>"
              ), 
              legendgroup = ~group, showlegend = TRUE) %>%
        add_trace(data = subset(plot_df_age, age >=100), x = ~age, y = ~intens, mode = 'lines+markers', showlegend = FALSE) %>%
        layout(xaxis = x_age2, yaxis = y_ratio2, shapes = list(hline(x0=age_min, x1=age_max), vline(y0=y_min, y1=y_max)), legend = list(x = 0.3, y = -0.2, orientation = 'h')
               ),
      
      titleX = TRUE,
      titleY = TRUE,
      nrows = 1
    )
  })
  
  output$protein_corr_plot <- renderPlotly({
    prot_df <- data.frame(t(protein_intensity_combat[which(rownames(exps_at_age_NDCHC) %in% input$proteins2), ]))
    #print(metadata)
    
    cor_e = cor(prot_df, use = "all.obs", method=input$cormethod)
    cor_p = cor_pmat(prot_df, use = "all.obs", method=input$cormethod, exact=FALSE)
    corr.plot <- ggcorrplot(
      cor_e, hc.order = TRUE, type = "lower", 
      outline.col = "white", show.diag = FALSE,
      p.mat = cor_p, sig.level = 0.05#,
      #title = "Pairwise correlations in ALL samples"
      #colors = c("#6D9EC1", "white", "#E46726")
    ) #+ labs(subtitle = "Pairwise correlations in ALL samples")
    
    cor_e_AD = cor(prot_df[AD_index, ], use = "all.obs", method=input$cormethod)
    cor_p_AD = cor_pmat(prot_df[AD_index, ], use = "all.obs", method=input$cormethod, exact=FALSE)
    corr.plot.AD <- ggcorrplot(
      cor_e_AD, hc.order = TRUE, type = "lower", 
      outline.col = "white", show.diag = FALSE,
      p.mat = cor_p_AD, sig.level = 0.05#,
      #title = "Pairwise correlations in AD group"
      #colors = c("#6D9EC1", "white", "#E46726")
    ) #+ labs(subtitle = "Pairwise correlations in AD samples")
    
    cor_e_ND = cor(prot_df[ND_index, ], use = "all.obs", method=input$cormethod)
    cor_p_ND = cor_pmat(prot_df[ND_index, ], use = "all.obs", method=input$cormethod, exact=FALSE)
    corr.plot.ND <- ggcorrplot(
      cor_e_ND, hc.order = TRUE, type = "lower", 
      outline.col = "white", show.diag = FALSE,
      p.mat = cor_p_ND, sig.level = 0.05#,
      #title = "Pairwise correlations in ND group"
      #colors = c("#6D9EC1", "white", "#E46726")
    ) #+ labs(subtitle = "Pairwise correlations in ND samples")
    
    cor_e_CHC = cor(prot_df[CHC_index, ], use = "all.obs", method=input$cormethod)
    cor_p_CHC = cor_pmat(prot_df[CHC_index, ], use = "all.obs", method=input$cormethod, exact=FALSE)
    corr.plot.CHC <- ggcorrplot(
      cor_e_CHC, hc.order = TRUE, type = "lower", 
      outline.col = "white", show.diag = FALSE,
      p.mat = cor_p_CHC, sig.level = 0.05#,
      #title = "Pairwise correlations in Centenarian group"
      #colors = c("#6D9EC1", "white", "#E46726")
    ) #+ labs(subtitle = "Pairwise correlations in Centenarian samples")
    #print(cor_p_CHC)
    
    subplot(ggplotly(corr.plot), ggplotly(corr.plot.AD), ggplotly(corr.plot.ND), ggplotly(corr.plot.CHC), 
            margin = 0.05,
            titleX = TRUE,
            titleY = TRUE,
            nrows = 2) %>% 
      layout(annotations = list( 
              list(x = 0 , y = 1.03, text = "Pairwise correlations in ALL samples", showarrow = F, xref='paper', yref='paper'), 
              list(x = 0.845 , y = 1.03, text = "Pairwise correlations in AD group", showarrow = F, xref='paper', yref='paper'),
              list(x = 0 , y = 0.46, text = "Pairwise correlations in ND group", showarrow = F, xref='paper', yref='paper'),
              list(x = 0.93 , y = 0.46, text = "Pairwise correlations in Centenarian group", showarrow = F, xref='paper', yref='paper')
              ) 
            ) 
    
  })
  
  output$protein_APOE_plot <- renderPlotly({
    protein = ""
    if (!(input$protein_APOE %in% protein_display_dataset$protein)){
      protein <- pre_protein
    }
    else{
      protein <- input$protein_APOE
      pre_protein <<- input$protein_APOE
    }
    
    plot_matrix <- data.frame(braak = metadata$Braak_text,
                              Abeta = metadata$amyloid,
                              group = metadata$group_n,
                              age = as.numeric(metadata$age),
                              diagnosis = as.character(metadata$diagnosis),
                              APOE = as.character(metadata$hasAPOEe4),
                              intensity = as.numeric(protein_intensity_combat[which(rownames(protein_intensity_combat) == protein),]))
    #print(protein)
    #print(plot_matrix)
    subplot(
      plot_ly(plot_matrix, x = ~group, y = ~intensity, color = ~APOE, text = ~diagnosis, type = "box", boxpoints = "all", pointpos = 0,
              colors = c("#1F77B4", "#FF7F0E"), 
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ), height=1000,
              legendgroup = ~APOE, showlegend = TRUE) %>%
        layout(xaxis = x_group, yaxis = y, boxmode = "group", 
               legend = list(x = 0.2, y = -0.1, orientation = 'h')),
      
      plot_ly(data = plot_matrix, x = ~Abeta, y = ~intensity, color = ~group, text = ~diagnosis, type = 'scatter',
              mode = 'markers', symbol = ~APOE, symbols = c('circle', 'o'),
              colors = c(color_AD, color_ND, color_CHC), marker = list(size = 10), source = "C",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ), height=1000,
              legendgroup = ~group, 
              showlegend = TRUE)%>%
        layout(xaxis = x_Abeta, yaxis = y, legend = list(x = 0.2, y = -0.1, orientation = 'h')),
      
      plot_ly(data = plot_matrix, x = ~braak, y = ~intensity, color = ~group, text = ~diagnosis, type = 'scatter',
              mode = 'markers', symbol = ~APOE, symbols = c('circle', 'o'),
              colors = c(color_AD, color_ND, color_CHC), marker = list(size = 10), source = "C",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ), height=1000,
              legendgroup = ~group, 
              showlegend = FALSE) %>%
        layout(xaxis = x_braak, yaxis = y_2),
      
      plot_ly(data = plot_matrix, x = ~age, y = ~intensity, color = ~group, text = ~diagnosis,type = 'scatter',
              mode = 'markers', symbol = ~APOE, symbols = c('circle', 'o'),
              colors = c(color_AD, color_ND, color_CHC), marker = list(size = 10), source = "C",
              hovertemplate = paste(
                "<b>&nbsp;%{text}</b><br>",
                "X: %{x:.2f}<br>",
                "Y: %{y:.2f}",
                "<extra></extra>"
              ), height=1000,
              legendgroup = ~group, 
              showlegend = FALSE) %>%
        layout(xaxis = x_age, yaxis = y_2),
      
      margin = 0.05,
      
      titleX = TRUE,
      titleY = TRUE,
      nrows = 2
    )
    
  })
  
  title_rbraak_ADND <- list(
    text = "In AD and ND groups",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = -0.5,
    y = 1,
    yshift=25, 
    showarrow = FALSE
  )
  title_rbraak_CEN <- list(
    text = "In Centenarian group",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = -0.5,
    y = 1,
    yshift=25, 
    showarrow = FALSE
  )
  title_rage_NDCHC <- list(
    text = "In ND and Centenarian groups",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = 64,
    y = 1,
    yshift=25,
    showarrow = FALSE
  )
  title_rage_AD <- list(
    text = "In AD group",
    font = f,
    xref = "x",
    yref = "paper",
    xanchor = 'left',
    x = 64,
    y = 1,
    yshift=25,
    showarrow = FALSE
  )
  
  hline <- function(y = 1, x0 = "0/I", x1 = "VI", color = "gray") {
    list(
      type = "line", 
      x0 = x0, 
      x1 = x1, 
      #xref = "x",
      y0 = y, 
      y1 = y, 
      line = list(color = color)
    )
  }
  
  vline <- function(y0 = 0, y1 = 1, x = 100, color = "red") {
    list(
      type = "line", 
      x0 = x, 
      x1 = x, 
      y0 = y0, 
      y1 = y1, 
      line = list(color = color, dash="dash", width=3)
    )
  }  
  
  output$line_plot <- renderPlotly({
    exps_at_stage_ADND_sub <- exps_at_stage_ADND[which(rownames(exps_at_stage_ADND) %in% input$proteins), ]
    braak  <- colnames(exps_at_stage_ADND_sub)
    braak <- rep(braak, nrow(exps_at_stage_ADND_sub))
    intens <- as.vector(t(exps_at_stage_ADND_sub))
    protein <- rep(rownames(exps_at_stage_ADND_sub), each = ncol(exps_at_stage_ADND))
    plot_df_Braak_ADND <- data.frame(protein = protein,
                          intens = intens,
                          braak = braak)
    plot_df_Braak_ADND$protein <- factor(plot_df_Braak_ADND$protein, levels = input$proteins)
    plot_df_Braak_ADND$braak <- factor(plot_df_Braak_ADND$braak, levels = c("0/I", "II", "III", "IV", "V", "VI"))
    
    exps_at_stage_CHC_sub <- exps_at_stage_CHC[which(rownames(exps_at_stage_CHC) %in% input$proteins), ]
    braak  <- colnames(exps_at_stage_CHC_sub)
    braak <- rep(braak, nrow(exps_at_stage_CHC_sub))
    intens <- as.vector(t(exps_at_stage_CHC_sub))
    protein <- rep(rownames(exps_at_stage_CHC_sub), each = ncol(exps_at_stage_CHC))
    plot_df_Braak_CHC <- data.frame(protein = protein,
                          intens = intens,
                          braak = braak)
    plot_df_Braak_CHC$protein <- factor(plot_df_Braak_CHC$protein, levels = input$proteins)
    plot_df_Braak_CHC$braak <- factor(plot_df_Braak_CHC$braak, levels = c("0/I", "II", "III", "IV", "V", "VI"))
    
    exps_at_age_NDCHC_sub <- exps_at_age_NDCHC[which(rownames(exps_at_age_NDCHC) %in% input$proteins), ]
    age  <- as.numeric(colnames(exps_at_age_NDCHC_sub))
    age <- rep(age, nrow(exps_at_age_NDCHC_sub))
    intens <- as.vector(t(exps_at_age_NDCHC_sub))
    protein <- rep(rownames(exps_at_age_NDCHC_sub), each = ncol(exps_at_age_NDCHC))
    plot_df_age_NDCHC <- data.frame(protein = protein,
                              intens = intens,
                              age = age)
    plot_df_age_NDCHC$protein <- factor(plot_df_age_NDCHC$protein, levels = input$proteins)
    
    exps_at_age_AD_sub <- exps_at_age_AD[which(rownames(exps_at_age_AD) %in% input$proteins), ]
    age  <- as.numeric(colnames(exps_at_age_AD_sub))
    age <- rep(age, nrow(exps_at_age_AD_sub))
    intens <- as.vector(t(exps_at_age_AD_sub))
    protein <- rep(rownames(exps_at_age_AD_sub), each = ncol(exps_at_age_AD))
    plot_df_age_AD <- data.frame(protein = protein,
                              intens = intens,
                              age = age)
    plot_df_age_AD$protein <- factor(plot_df_age_AD$protein, levels = input$proteins)
    
    y_min <- min(min(plot_df_Braak_ADND$intens), min(plot_df_Braak_CHC$intens), min(plot_df_age_NDCHC$intens), min(plot_df_age_AD$intens)) - 0.05
    y_max <- max(max(plot_df_Braak_ADND$intens), max(plot_df_Braak_CHC$intens), max(plot_df_age_NDCHC$intens), max(plot_df_age_AD$intens)) + 0.05
    
    y_ratio <- list(
      title = "LFQ intensity ratio",
      range = list(y_min, y_max), 
      zeroline = FALSE,
      linecolor = toRGB("black"),
      linewidth = 0.5
    )
    
    y_ratio2 <- list(
      title = "",
      range = list(y_min, y_max), 
      zeroline = FALSE,
      linecolor = toRGB("black"),
      linewidth = 0.5
    )
    
    #print(plot_df_age)
      subplot(
        plot_ly(plot_df_Braak_ADND, x = ~braak, y = ~intens, text = ~protein, type = 'scatter', mode = 'lines+markers', #linetype = ~protein, 
                color = ~protein,
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "(%{x}, %{y:.2f})",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~protein, showlegend = TRUE) %>%
          layout(xaxis = x_braak, yaxis = y_ratio, annotations = title_rbraak_ADND, shapes = list(hline()), height = height0, width = width0, #margin = list(l = 0, r = 0, t = 0, b = tb_margin),
                 legend = list(x = 0, y = -0.1, orientation = 'h')
          ),
        
        plot_ly(plot_df_Braak_CHC, x = ~braak, y = ~intens, text = ~protein, type = 'scatter', mode = 'lines+markers', #linetype = ~protein, 
                color = ~protein,
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "(%{x}, %{y:.2f})",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~protein, showlegend = FALSE) %>%
          layout(xaxis = x_braak, yaxis = y_ratio2, annotations = title_rbraak_CEN, height = height0, width = width0, #margin = list(l = 0, r = 0, t = 0, b = tb_margin),
                 shapes = list(hline(x1 = "VI"))
          ),
        
        plot_ly(plot_df_age_AD, x = ~age, y = ~intens, text = ~protein, type = 'scatter', mode = 'lines+markers', #linetype = ~protein, 
                color = ~protein,
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "(%{x}, %{y:.2f})",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~protein, showlegend = FALSE) %>%
          layout(xaxis = x_age2, yaxis = y_ratio, annotations = title_rage_AD, height = height0, width = width0, #margin = list(l = 0, r = 0, t = 0, b = 0),
                 shapes = list(hline(x0=age_min, x1=age_max)#, vline(y0=y_min, y1=y_max)
                                                                                           )),
        
        plot_ly(subset(plot_df_age_NDCHC, age <100), x = ~age, y = ~intens, text = ~protein, type = 'scatter', mode = 'lines+markers', #linetype = ~protein, 
                color = ~protein,
                hovertemplate = paste(
                  "<b>&nbsp;%{text}</b>",
                  "(%{x}, %{y:.2f})",
                  "<extra></extra>"
                ), height=600,
                legendgroup = ~protein, showlegend = FALSE) %>%
          add_trace(data = subset(plot_df_age_NDCHC, age >=100), x = ~age, y = ~intens, text = ~protein, mode = 'lines+markers', #linetype = ~protein, 
                    color = ~protein) %>%
          layout(xaxis = x_age2, yaxis = y_ratio2, annotations = title_rage_NDCHC, height = height0, width = width0, #margin = list(l = 0, r = 0, t = 0, b = 0),
                 shapes = list(hline(x0=age_min, x1=age_max), vline(y0=y_min, y1=y_max)
                 )),
        
        margin = 0.065,
        titleX = TRUE,
        titleY = TRUE,
        nrows = 2
      )
    })
  
  output$distribution_plot <- renderPlot({
    x <- NULL
    y <- "None"
    if(input$variable1 == "Age")
      x <- "age"
    if(input$variable1 == "Braak stage")
      x <- "fBraak"
    if(input$variable1 == "Amyloid stage")
      x <- "amyloid"
    if(input$variable1 == "Sex")
      x <- "sex"
    if(input$variable1 == "APOE allele")
      x <- "APOE"
    
    if(input$variable2 == "Age")
      y <- "age"
    if(input$variable2 == "Braak stage")
      y <- "fBraak"
    if(input$variable2 == "Amyloid stage")
      y <- "amyloid"
    if(input$variable2 == "Sex")
      y <- "sex"
    if(input$variable2 == "APOE allele")
      y <- "APOE"
    
    #ggthemr("pale", layout = "clear", spacing = 1.6,
     #       text_size = 12, type = "outer", line_weight = 0.5)
    
    if((y == "None") | (x == y)){
      plot_distr_group <- metadata[,c(x, "group_n")]
      colnames(plot_distr_group) <- c("variable", "group")
      
      if(input$variable1 == "Age"){
        ggplot(plot_distr_group, aes(x=variable, fill=group)) + 
          geom_histogram(binwidth=0.5, alpha=1, position="stack") +
          xlab(input$variable1) + ylab("Number of cases") +
          scale_fill_manual(name=NULL, values = c(color_AD, color_ND, color_CHC))+
          theme_light() +
          theme(legend.position="bottom")

      }
      else{
        ggplot(plot_distr_group, aes(x=as.factor(variable), fill=group)) + 
          geom_histogram(binwidth=0.5, alpha=1, position="stack", stat="count") +
          xlab(input$variable1) + ylab("Number of cases") +
          scale_fill_manual(name=NULL, values = c(color_AD, color_ND, color_CHC))+
          theme_light() +
          theme(legend.position="bottom")
      }
    }
    else{
      plot_distr_group <- metadata[,c(x, y, "group_n")]
      colnames(plot_distr_group) <- c("variable1", "variable2", "group")
      ggplot(plot_distr_group, aes(x=variable1, y=variable2, color=group)) + 
        geom_point(position=position_jitter(h=0.1, w=0.1), size = 3, alpha = 0.5) +
        xlab(input$variable1) + ylab(input$variable2) +
        scale_color_manual(name=NULL, values = c(color_AD, color_ND, color_CHC))+
        theme_light() +
        theme(legend.position="bottom")
      
    }
  })
  
  
  observeEvent(event_data("plotly_click", source = "A"), {
    updateSearchInput(
      session = session,
      inputId = "protein_symbol",
      value = event_data("plotly_click")$key,
      trigger = TRUE
    )
  }, ignoreInit = TRUE)
  
})
