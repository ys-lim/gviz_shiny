source("/mnt/gtklab01/yongshan/results/intermediate_files/gviz_shiny/browser.R")

ui <- fluidPage(
  theme = shinytheme("lumen"),
  titlePanel(h1("RNA-seq Visualisation")),
  mainPanel(width = 12,
            browser_ui("browser")
  )
)

server <- function(input,output,session) {
  callModule(browser_server, "browser")
}

shinyApp(ui, server)
