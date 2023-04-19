library(shiny)
library(shinyMatrix)
library(shinyjs)
library(rsconnect)

# Load the stochastic.R file
source("stochastic.R")

# Default values
drugs_default <- matrix(c(
    -1, -1, # log_10 mutation rates
    -1, -1, # log_10 drug elimination rates
    10, 10 # drug influx concentrations
    ), nrow = 2, ncol = 3, dimnames = list(c("A1", "A2"),
        c("mutation_rate", "elimination_rate", "influx")))

growth_default <- matrix(c(
    1e2, 0, 0, 0, # init: initial populations
    1, 1, 1, 1, # mu: growth rates
    1e2, 1e2, 1e2, 1e2, # k: half-maximal growth rates
    1, 1, 1, 1 # alpha: resources used per unit growth
    ), nrow = 4, ncol = 4,
    dimnames = list(c("S", "R1", "R2", "R12"),
        c("init", "mu", "k", "alpha")))

pd_default <- matrix(c(
    0.1, 0.1, 0.1, 0.1, # psi
    0.2, 0.2, 0.2, 0.2, # phi1
    1, 100, 1, 100, # zeta1
    1, 1, 1, 1, # kappa1
    0.2, 0.2, 0.2, 0.2, # phi2
    1, 1, 100, 100, # zeta2
    1, 1, 1, 1, # kappa2
    0, 0, 0, 0 # theta
    ), nrow = 4, ncol = 8,
    dimnames = list(c("S", "R1", "R2", "R12"),
    c("psi", "phi1", "zeta1", "kappa1", "phi2", "zeta2", "kappa2", "theta")))

# UI
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Stochastic Simulation of Antibiotic Resistance"),
  div(id = "everything",
  mainPanel(
    tabsetPanel(
      tabPanel("Basic",
               fluidRow(
                 column(6, wellPanel(
                   numericInput("rep", "Number of Runs", value = 10, min = 1, step = 1),
                   checkboxInput("pharmacokinetic", "Pharmacokinetic Model", FALSE),
                   selectInput("stewardship", "Antibiotic Stewardship Strategy",
                               choices = c("Cycling" = "cycl",
                                           "Combination" = "comb",
                                           "Drug 1 Only" = "1_only",
                                           "Drug 2 Only" = "2_only")),
                   selectInput("display", "Uncertainty Display",
                               choices = c("Median + 25th and 75th percentiles" = TRUE,
                                           "Mean + 95% confidence interval" = FALSE))
                 )),
                 column(6, wellPanel(
                   sliderInput("time", "Simulation Time (hours)",
                    value = 50, min = 10, max = 100, step = 10),
                   sliderInput("freq", "Bottleneck Frequency (hours)",
                    value = 10, min = 1, max = 30, step = 1),
                   sliderInput("dt", "log_10 of time between data points (hours)",
                    value = -1, min = -3, max = 0, step = 1),
                   sliderInput("HGT", "log_10 of recombination rate",
                    value = -4, min = -9, max = -3, step = 1)
                 ))
               )
      ),
      tabPanel("Drugs",
               fluidRow(
                 column(12, wellPanel(
                    matrixInput("drugs", value = drugs_default,
                    rows = list(extend = FALSE, names = TRUE),
                    cols = list(extend = FALSE, names = TRUE),
                    class = "numeric")
                )))
      ),
      tabPanel("Growth & Bottlenecks",
         fluidRow(
           column(6, wellPanel(
            matrixInput("growth", value = growth_default,
            rows = list(extend = FALSE, names = TRUE),
            cols = list(extend = FALSE, names = TRUE),
            class = "numeric")
              )),
           column(6, wellPanel(
             sliderInput("N0", "log_10 of Nutrient influx at bottleneck",
                value = 4, min = 1, max = 8, step = 1),
             sliderInput("D", "log_10 of Dilution fraction at bottleneck",
                value = -1, min = -5, max = 0, step = 1)
           ))
         )
      ),
      tabPanel("Pharmacodynamics",
         fluidRow(
           column(12, wellPanel(
             matrixInput("pd",
                         value = pd_default,
                         rows = list(extend = FALSE, names = TRUE),
                         cols = list(extend = FALSE, names = TRUE),
                         class = "numeric")
             ))))
      # Add other tabs here
    ),
    plotOutput("simulation_plot")
  )),
               actionButton("run_simulation", "Run Simulation"),
               actionButton("reset_all", "Reset all parameters to defaults")
)

server <- function(input, output, session) {
  # Create a reactive expression for the simulation results
  simulation_result <- eventReactive(input$run_simulation, {
    simulate(
        rep = input$rep,
        pharmacokinetic = input$pharmacokinetic,
        stewardship = input$stewardship,
        time = input$time,
        freq = input$freq,
        dt = 10 ^ input$dt,
        N0 = round(10 ^ input$N0),
        D = 10 ^ input$D,
        HGT = 10 ^ input$HGT,
        m1 = 10 ^ input$drugs["A1", "mutation_rate"],
        m2 = 10 ^ input$drugs["A2", "mutation_rate"],
        d1 = input$pharmacokinetic * 10 ^ input$drugs["A1", "elimination_rate"],
        d2 = input$pharmacokinetic * 10 ^ input$drugs["A2", "elimination_rate"],
        influx = setNames(input$drugs[, "influx"], c("A1", "A2")),
        init = setNames(input$growth[, c("init")], c("S", "R1", "R2", "R12")),
        psi = input$pd[, "psi"],
        phi1 = input$pd[,"phi1"],
        zeta1 = input$pd[, "zeta1"],
        kappa1 = input$pd[, "kappa1"],
        phi2 = input$pd[, "phi2"],
        zeta2 = input$pd[, "zeta2"],
        kappa2 = input$pd[, "kappa2"],
        theta = input$pd[, "theta"],
        mu = input$growth[, "mu"],
        k = input$growth[, "k"],
        alpha = input$growth[, "alpha"]
    )
  })

  output$simulation_plot <- renderPlot({
    # Check if the simulation_result has been executed
    if (!is.null(simulation_result())) {
      # Create a plot of the simulation results
      log_plot(summarise(simulation_result()), IQR = input$display)
    }
  })
  # add code to reset all parameters to defaults
    observeEvent(input$reset_all, {
        shinyjs::reset("everything")
        # update the matrices to use default values
        updateMatrixInput(session, inputId = "drugs", value = drugs_default)
        updateMatrixInput(session, inputId = "growth", value = growth_default)
        updateMatrixInput(session, inputId = "pd", value = pd_default)
    })
}

# Run the application
shinyApp(ui = ui, server = server)

# To deploy the app use deployApp()
