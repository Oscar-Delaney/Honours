library(shiny)
library(shinyMatrix)
library(shinyjs)
library(rsconnect)

# Load the stochastic.R file
source("stochastic.R")

# Explanatory text
basic_text <- "This is a model describing the evolution of antibiotic resistance
in a bacterial population. It uses the Adaptive Tau package to stochastify
an ODE model. For details see Delaney, Engelstaedter, and Letten (forthcoming).
Contact Oscar Delaney on o.delaney@uq.net.au with any errors or suggestions."

drugs_text <- "A1 and A2 are arbitrary antibiotics. The mutation rate
is the proportion of genome replications that result in
resistance to that drug. The elimination rate is the rate at which the drug
degenerates in the body, in units of the MIC. This is only used if
the pharmacokinetic model is selected. The influx is the concentration of
each drug added at bottleneck event."

growth_text <- "The rows represent four bacterial strains: Susceptible, 
Antibiotic_1 resistant, Antibiotic_2 resistant, and double resistant. 
init is the starting population size of each straing. Mu is the growth
rate with unlimited nutrients. K is the nutrient concentration that produces 
half-maximal growth rate. alpha is the amount of nutrient used per new
bacterial cell."

pd_text <- "Psi is the growth rate in the absence of antibiotics.
Phi_i is the maximum reduction in growth rate caused by antibiotic i.
Zeta_i is the minimum inhibitory concentration of antibiotic i.
Kappa_i is the shape parameter, with smaller values representing a
shallower reduction in growth rate as antibiotic concentration increases, and 
larger values representing a steeper reduction in growth rate around the MIC.
Theta is the drug-drug interaction parameter, with values between -1 and 1."

# Default values
drugs_default <- matrix(c(
    1e-9, 1e-9, # mutation rates
    log(2) / 3.5, log(2) / 3.5, # drug elimination rates
    7, 7 # drug influx concentrations, MIC units
), nrow = 2, ncol = 3, dimnames = list(
    c("A1", "A2"),
    c("Mutation rate", "Elimination rate", "Influx")
))

growth_default <- matrix(
    c(
        1e10, 0, 0, 0, # init: initial populations
        0.8 * c(1, 0.9, 0.9, 0.81), # mu: growth rates
        1e2, 1e2, 1e2, 1e2, # k: half-maximal growth rates
        1, 1, 1, 1 # alpha: resources used per unit growth
    ),
    nrow = 4, ncol = 4,
    dimnames = list(
        c("S", "R1", "R2", "R12"),
        c("Init", "Mu", "K", "Alpha")
    )
)

pd_default <- matrix(
    c(
        0.3, 0.3, 0.3, 0.3, # psi
        0.6, 0.6, 0.6, 0.6, # phi1
        1, 28, 1, 28, # zeta1
        1, 1, 1, 1, # kappa1
        0.6, 0.6, 0.6, 0.6, # phi2
        1, 1, 28, 28, # zeta2
        1, 1, 1, 1, # kappa2
        0, 0, 0, 0 # theta
    ),
    nrow = 4, ncol = 8,
    dimnames = list(
        c("S", "R1", "R2", "R12"),
        c("Psi", "Phi1", "Zeta1", "Kappa1", "Phi2", "Zeta2", "Kappa2", "Theta")
    )
)

# UI
ui <- fluidPage(
    useShinyjs(),
    titlePanel("Stochastic Simulation of Antibiotic Resistance"),
    p(
    actionButton("run_simulation", "Run Simulation"),
    actionButton("reset_all", "Reset all parameters to defaults")
    ),
    div(
        id = "everything",
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "Basic",
                    fluidRow(
                        column(6, wellPanel(
                            p(basic_text),
                            numericInput("rep", "Number of Runs", value = 10, min = 1, step = 1),
                            numericInput("time", "Simulation Time (hours)", value = 100, step = 1),
                            selectInput("stewardship", "Antibiotic Stewardship Strategy",
                                choices = c(
                                    "Cycling" = "cycl",
                                    "Combination" = "comb",
                                    "Drug 1 Only" = "1_only",
                                    "Drug 2 Only" = "2_only"
                                )
                            ),
                            checkboxInput("pharmacokinetic", "Pharmacokinetic Model", FALSE),
                            checkboxInput("deterministic", "Deterministic Model", FALSE),
                            # sliderInput("dt", "log_10 of time between data points (hours)",
                            #     value = -1, min = -3, max = 0, step = 0.1
                            # ),
                            numericInput("HGT", "recombination rate",
                                value = 0, min = 0, max = 1, step = 1e-4
                            )
                        ))
                    )
                ),
                tabPanel(
                    "Drugs",
                    fluidRow(
                        column(6, wellPanel(
                            p(drugs_text),
                            matrixInput("drugs",
                                value = drugs_default,
                                rows = list(extend = FALSE, names = TRUE),
                                cols = list(extend = FALSE, names = TRUE),
                                class = "numeric"
                            )
                        ))
                    )
                ),
                tabPanel(
                    "Growth",
                    fluidRow(
                        column(6, wellPanel(
                            p(growth_text),
                            matrixInput("growth",
                                value = growth_default,
                                rows = list(extend = FALSE, names = TRUE),
                                cols = list(extend = FALSE, names = TRUE),
                                class = "numeric"
                            )
                        )),
                    )
                ),
                tabPanel(
                    "Bottlenecks",
                    fluidRow(
                        column(6, wellPanel(
                            sliderInput("freq", "Bottleneck Frequency (hours)",
                                value = 10, min = 1, max = 30, step = 1
                            ),
                            sliderInput("N0", "log_10 of Nutrient influx at bottleneck",
                                value = 15, min = 2, max = 20, step = 0.1
                            ),
                            sliderInput("D", "log_10 of Dilution fraction at bottleneck",
                                value = -1, min = -5, max = 0, step = 0.1
                            )
                        ))
                    )
                ),
                tabPanel(
                    "Pharmacodynamics",
                    fluidRow(
                        column(12, wellPanel(
                            p(pd_text),
                            matrixInput("pd",
                                value = pd_default,
                                rows = list(extend = FALSE, names = TRUE),
                                cols = list(extend = FALSE, names = TRUE),
                                class = "numeric"
                            )
                        ))
                    )
                ),
                tabPanel(
                    "Graph",
                    radioButtons("display", "Uncertainty Display",
                        choices = c(
                            "Median + 25th and 75th percentiles" = "median",
                            "Mean + 95% confidence interval" = "mean",
                            "Plot all the runs" = "all"
                        ),
                        selected = "median"
                    ),
                    plotOutput("simulation_plot", height = "800", width = "100%")
                )
                # Add other tabs here
            )
        )
    )
)

server <- function(input, output, session) {
    # Create a reactive expression for the simulation results
    simulation_result <- eventReactive(input$run_simulation, {
        simulate(
            rep = input$rep,
            pharmacokinetic = input$pharmacokinetic,
            deterministic = input$deterministic,
            stewardship = input$stewardship,
            time = input$time,
            freq = input$freq,
            # dt = 10^input$dt,
            N0 = round(10^input$N0),
            D = 10^input$D,
            HGT = input$HGT,
            m1 = input$drugs["A1", "Mutation rate"],
            m2 = input$drugs["A2", "Mutation rate"],
            d1 = input$pharmacokinetic * input$drugs["A1", "Elimination rate"],
            d2 = input$pharmacokinetic * input$drugs["A2", "Elimination rate"],
            influx = setNames(input$drugs[, "Influx"], c("A1", "A2")),
            init = setNames(input$growth[, c("Init")], c("S", "R1", "R2", "R12")),
            psi = input$pd[, "Psi"],
            phi1 = input$pd[, "Phi1"],
            zeta1 = input$pd[, "Zeta1"],
            kappa1 = input$pd[, "Kappa1"],
            phi2 = input$pd[, "Phi2"],
            zeta2 = input$pd[, "Zeta2"],
            kappa2 = input$pd[, "Kappa2"],
            theta = input$pd[, "Theta"],
            mu = input$growth[, "Mu"],
            k = input$growth[, "K"],
            alpha = input$growth[, "Alpha"]
        )
    })

    output$simulation_plot <- renderPlot({
        # Check if the simulation_result has been executed
        if (!is.null(simulation_result())) {
            # Create a plot of the simulation results
            log_plot(simulation_result(), type = input$display)
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