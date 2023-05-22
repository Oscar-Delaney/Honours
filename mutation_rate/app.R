library(shiny)
library(shinyjs)
library(rsconnect)
source("mutation_rate_evol.R")

intro_text <- "Population size is the initial number of individuals in the population, which is also the carrying capacity.
Initial mutation rate(s) is in log_10 format. You can input a single value, or a sequence of values in the format of R code, such as c(-2, -1.9, -1.2) or seq(-2, -1, 0.2).
Initial fitness score is the initial net number of beneficial mutations in the population, where fitness is (1+s) to the power of this quantity, so initially 1.
Mutation rate mutability is the proportion of genomic mutations that affect the mutation rate.
Mutation rate jump size represents the magnitude of mutation rate change (increase or decrease).
s is the magnitude of fitness advantage/disadvantage of a mutation.
Baseline reproduction rate is the number of offspring per individual per generation initially.
Maximum reproduction rate is the upper limit of reproduction rate in the population. In the deterministic model this can safely be set to Inf, but in the stochastic model Inf gives errors for a large enough number of generations.
If Deterministic model is checked, the model uses an infinite population size.
Minimum detectable population size is the smallest population size that is not rounded to 0. This is only relevant in the deterministic model, and values below about 1e-20 lead to noticeably slower runtimes.
The simulator generates three types of output: an average mutation rate plot, an average fitness score plot, and a heatmap of the final population distribution."


basic_content <- wellPanel(
    numericInput("size", "Population size", value = "1e6"),
    numericInput("generations", "Number of generations", value = "1e3"),
    textInput("init_mu", "Initial mutation rate(s):", "c(-2, -1.9, -1.2)"),
    numericInput("init_w", "Initial fitness score", value = 0),
    numericInput("a", "Mutation rate mutability", value = 0.1),
    numericInput("p_mu_up", "P(mutation rate increases)", value = 0.5),
    numericInput("jump", "Mutation rate jump size", value = 1e-1),
    numericInput("p_w_up", "P(mutation increases fitness)", value = 0.1),
    numericInput("s", "s", value = 1e-2),
    numericInput("r", "Baseline reproduction rate", value = 1),
    numericInput("r_max", "Maximum reproduction rate", value = "1e5"),
    checkboxInput("det", "Deterministic model", value = FALSE),
    numericInput("sensitivity", "Minimum detectable population size:", value = "1e-10")
)

# Define UI for application
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Mutation Rate Evolution Simulation"),
    p(
        actionButton("run", "Run Simulation"),
        actionButton("reset_all", "Reset all parameters to defaults")
        ),
        div(
            id = "everything",
            mainPanel(
                tabsetPanel(
                    tabPanel("Basic", fluidRow(column(6, basic_content),
                        column(6, p(intro_text)))),
                    tabPanel("Graph (Mutation rate)",
                        p("Mutation rate over time:"),
                        plotOutput("muPlot", height = "800", width = "100%")),
                    tabPanel("Graph (Fitness)",
                        p("Net beneficial mutations over time:"),
                        plotOutput("wPlot", height = "800", width = "100%")),
                    tabPanel("Graph (Heatmap)",
                        p("Heatmap of mutation rate and fitness in endpoint population:"),
                        plotOutput("heatmap", height = "800", width = "100%"))
                )
            )
        )
    )

# Define server logic
server <- function(input, output) {
   observeEvent(input$run, {
      results <- evolve(
         size = input$size,
         init_mu = eval(parse(text = input$init_mu)),
         init_w = input$init_w,
         a = input$a,
         p_mu_up = input$p_mu_up,
         jump = input$jump,
         p_w_up = input$p_w_up,
         generations = input$generations,
         s = input$s,
         r = input$r,
         r_max = input$r_max,
         sensitivity = input$sensitivity,
         det = input$det
      )

      pop <- results$pop
      stats <- results$stats

      output$muPlot <- renderPlot({
         plot_univariate(stats, "mu", "log_10 Average Mutation Rate")
      })

      output$wPlot <- renderPlot({
         plot_univariate(stats, "w", "Average Fitness Score")
      })

      output$heatmap <- renderPlot({
         plot_heatmap(pop)
      })
   })
   # add code to reset all parameters to defaults
    observeEvent(input$reset_all, {
        shinyjs::reset("everything")
    })
}

# Run the application
shinyApp(ui = ui, server = server)

# To deploy the app to shinyapps.io use deployApp()