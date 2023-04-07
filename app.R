library(shiny)
library(rsconnect)

# Load the stochastic.R file
source("stochastic.R")

# UI
ui <- fluidPage(
  titlePanel("Stochastic Simulation of Antibiotic Resistance"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("rep", "Number of Runs", value = 10, min = 1, step = 1),
      checkboxInput("pharmacokinetic", "Pharmacokinetic Model", FALSE),
      selectInput("stewardship", "Antibiotic Stewardship Strategy",
                  choices = c("Cycling" = "cycl",
                              "Combination" = "comb",
                              "Drug 1 Only" = "1_only",
                              "Drug 2 Only" = "2_only")),
      numericInput("time", "Simulation Time (hours)", value = 50, min = 1),
      numericInput("freq", "Bottleneck Frequency (hours)", value = 10, min = 1),
      sliderInput("N0", "log_10 of Nutrient influx as bottleneck", value = 4,
        min = 1, max = 8, step = 1),
      sliderInput("D", "log_10 of Dilution fraction at bottleneck", value = -1,
        min = -5, max = 0, step = 1),
      sliderInput("m1", "mutation rate for A1 resistance", value = -2,
        min = -9, max = -2, step = 1),
      sliderInput("m2", "mutation rate for A2 resistance", value = -2,
        min = -9, max = -2, step = 1),
      sliderInput("HGT", "recombination rate", value = -4,
        min = -9, max = -3, step = 1),
      sliderInput("cost1", "fitness cost of A1 resistance", value = 0,
        min = 0, max = 0.1, step = 0.01),
      sliderInput("cost2", "fitness cost of A2 resistance", value = 0,
        min = 0, max = 0.1, step = 0.01),
      sliderInput("theta", "drug-drug interaction", value = 0,
        min = -1, max = 1, step = 0.1),
      actionButton("run_simulation", "Run Simulation")
    ),
    mainPanel(
      plotOutput("simulation_plot")
    )
  )
)

server <- function(input, output) {
  # Create a reactive expression for the simulation results
  simulation_result <- eventReactive(input$run_simulation, {
    simulate(
      rep = input$rep,
      pharmacokinetic = input$pharmacokinetic,
      stewardship = input$stewardship,
      time = input$time,
      freq = input$freq,
      N0 = round(10 ^ input$N0),
      D = 10 ^ input$D,
      m1 = 10 ^ input$m1,
      m2 = 10 ^ input$m2,
      HGT = 10 ^ input$HGT,
      psi = c(0.1, 0.1 - input$cost1, 0.1 - input$cost2,
        0.1 - input$cost1 - input$cost2),
      theta = rep(input$theta, 4)
    )
  })

  output$simulation_plot <- renderPlot({
    # Check if the simulation_result has been executed
    if (!is.null(simulation_result())) {
      # Create a plot of the simulation results
      log_plot(summarise(simulation_result()))
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)