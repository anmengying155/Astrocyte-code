library(shiny)
library(ggplot2)
library(dplyr)


load("D:/File/document/单细胞测序代码/AUCell_tutorial/plot_df.RData") 

ui <- fluidPage(
  titlePanel("选择感兴趣的 gene sets"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("selected_genes", "请选择通路（gene set）:",
                         choices = unique(plot_df$gene_set),
                         selected = unique(plot_df$gene_set)[1:3]) # 默认选前3个
    ),
    mainPanel(
      plotOutput("dotplot")
    )
  )
)

server <- function(input, output) {
  output$dotplot <- renderPlot({
    req(input$selected_genes)
    
    # 过滤数据
    plot_data <- plot_df %>% filter(gene_set %in% input$selected_genes)
    
    ggplot(plot_data, aes(x = gene_set, y = mean_auc, color = mean_auc, size = log10_p)) +
      geom_point() +
      facet_wrap(~orig.ident + celltype, scales = "free_y") +
      scale_color_gradient(low = "#ffebeb", high = "#ff6767") +
      theme_minimal(base_size = 13) +
      labs(
        title = "Selected Gene Sets",
        x = "Gene Set", y = "Mean AUCell",
        color = "Mean AUCell", size = "-log10(p)"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 12))
  })
}

shinyApp(ui = ui, server = server)
