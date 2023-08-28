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

drugs_text <- "A and B are arbitrary antibiotics. The mutation rate
is the proportion of genome replications that result in resistance to 
that drug. The elimination rate is the rate at which the drug degenerates
in the body, in units of hours^-1. The influx is the concentration of
each drug added at bottleneck events, in units of zeta (see pharmacodynamics).
The recombination rate determines the rate at which the N_S + N_AB <--> N_A + N_B
reversible reaction occurs, with reasonable values being <1e-13 or so."

growth_text <- "The rows represent four bacterial strains: Susceptible, 
A-resistant, B-resistant, and double resistant. 
Init is the starting population size of each strain. Mu is the growth
rate with unlimited resources. K is the resource concentration that produces 
half the maximal growth rate. Alpha is the amount of resources used per new
bacterial cell."

events_text <- "The two possible time-discontinuous features of the model
are dosing with antibiotics, and bottlenecks where the population size is
reduced proportionally and new resources are added. Supply determines the rate
at which resources are continuously dded to the system. Drugs can be set to
disappear when the next dose arrives, which is less realistic but sometimes
convenient, or to persist through dosing."

pd_text <- "Bcidal_i is the maximum death rate caused by drug i. Bstatic_i is
the maximum proportional reduction in growth rate caused by drug i. Zeta_i is
the concentration of antibiotic i resulting in half the maximal death rate.
Kappa_i is the shape parameter, with smaller values representing a steep
initial increase in death rate, and shallow gradient about zeta_i.
Delta is the intrinsic death rate independent of drugs."

# Default values
drugs_default <- matrix(
    c(
        1e-9, 1e-9, # mutation rates
        rep(round(log(2) / 3.5, 3), 2), # drug elimination rates
        7, 7 # drug influx concentrations, MIC units
    ),
    nrow = 2, ncol = 3,
    dimnames = list(
        c("A", "B"),
        c("Mutation rate", "Elimination rate", "Influx")
    )
)

growth_default <- matrix(
    c(
        "1e+9", 0, 0, 0, # init: initial populations
        0.88 * c(1, 0.9, 0.9, 0.81), # mu: growth rates
        c(1e14, 1e14, 1e14, 1e14), # k: resources at half-maximal growth rate
        1, 1, 1, 1 # alpha: resources used per unit growth
    ),
    nrow = 4, ncol = 4,
    dimnames = list(
        c("N_S", "N_A", "N_B", "N_AB"),
        c("Init", "Mu", "K", "Alpha")
    )
)

pd_default <- matrix(
    c(
        1, 1, 1, 1, # Bcidal1
        0, 0, 0, 0, # Bstatic1
        1, 28, 1, 28, # zeta_A
        1, 1, 1, 1, # kappa_A
        1, 1, 1, 1, # Bcidal2
        0, 0, 0, 0, # Bstatic2
        1, 1, 28, 28, # zeta_B
        1, 1, 1, 1, # kappa_B
        0, 0, 0, 0 # delta
    ),
    nrow = 4, ncol = 9,
    dimnames = list(
        c("N_S", "N_A", "N_B", "N_AB"),
        c("Bcidal1", "Bstatic1", "Zeta1", "Kappa1", "Bcidal2", "Bstatic2", "Zeta2", "Kappa2", "Delta")
    )
)

basic_content <- wellPanel(
    p(basic_text),
    numericInput("rep", "Number of Runs", value = 1, min = 1, step = 1),
    numericInput("time", "Simulation Time (hours)", value = 100, step = 1),
    numericInput("dt", "Granularity (hours)", value = 0.1, step = 0.1),
    checkboxInput("deterministic", "Deterministic Model", FALSE),
    checkboxInput("cycl", "Cycle between drugs", TRUE),
    numericInput("seed", "Random Seed", value = NULL)
)

drugs_content <- wellPanel(
    p(drugs_text),
    matrixInput("drugs", value = drugs_default, class = "numeric"),
    numericInput("bactericidal", "Drug mode (bactericidal = 1 to
     bacteriostatic = 0)", value = 1, min = 0, max = 1, step = 1),
    numericInput("HGT", "recombination rate", value = 0, min = 0, step = 1e-15)
)

growth_content <- wellPanel(
    p(growth_text),
    matrixInput("growth", value = growth_default, class = "numeric")
)

events_content <- wellPanel(
    p(events_text),
    numericInput("tau", "Bottleneck frequency (hours)", value = 10, step = 1),
    numericInput("R0", "Bottleneck Resource pulse", value = "1e15", step = 1e9),
    numericInput("supply", "Resource supply rate (hours^-1)", value = 0, step = 1e9),
    numericInput("D", "Bottleneck Dilution fraction", value = 1e-1, step = 0.1),
    numericInput("dose_gap", "Gap between doses (hours)", value = 10, step = 1),
    numericInput("dose_rep", "Drug cycling period (doses)", value = 1, step = 1),
    checkboxInput("keep_old_drugs", "Old drugs persist through dosing", TRUE),
)

pd_content <- wellPanel(
    p(pd_text),
    matrixInput("pd", value = pd_default, class = "numeric")
)

graph_content <- tagList(
    radioButtons("display", "Uncertainty Display",
        choices = c(
            "Median + 25th and 75th percentiles" = "median",
            "Mean + 95% confidence interval" = "mean",
            "Plot all the runs" = "all"
        ),
        selected = "all"
    ),
    plotOutput("plot", height = "800", width = "100%")
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
                tabPanel("Basic", fluidRow(column(6, basic_content))),
                tabPanel("Drugs", fluidRow(column(6, drugs_content))),
                tabPanel("Growth", fluidRow(column(6, growth_content))),
                tabPanel("Events", fluidRow(column(6, events_content))),
                tabPanel("Pharmacodynamics", fluidRow(column(12, pd_content))),
                tabPanel("Graph", graph_content)
            )
        )
    )
)

server <- function(input, output, session) {
    # Create a reactive expression for the simulation results
    simulation_result <- eventReactive(input$run_simulation, {
        simulate(
            rep = input$rep,
            deterministic = input$deterministic,
            seed = input$seed,
            cycl = input$cycl,
            dose_rep = input$dose_rep,
            dose_gap = input$dose_gap,
            keep_old_drugs = input$keep_old_drugs,
            time = input$time,
            tau = input$tau,
            dt = input$dt,
            R0 = input$R0,
            D = input$D,
            HGT = input$HGT,
            m_A = input$drugs["A", "Mutation rate"],
            m_B = input$drugs["B", "Mutation rate"],
            d_A = input$drugs["A", "Elimination rate"],
            d_B = input$drugs["B", "Elimination rate"],
            influx = setNames(input$drugs[, "Influx"], c("C_A", "C_B")),
            init = input$growth[, c("Init")],
            bcidal_A = input$pd[, "Bcidal1"],
            bstatic_A = input$pd[, "Bstatic1"],
            zeta_A = input$pd[, "Zeta1"],
            kappa_A = input$pd[, "Kappa1"],
            bcidal_B = input$pd[, "Bcidal2"],
            bstatic_B = input$pd[, "Bstatic2"],
            zeta_B = input$pd[, "Zeta2"],
            kappa_B = input$pd[, "Kappa2"],
            delta = input$pd[, "Delta"],
            mu = input$growth[, "Mu"],
            k = input$growth[, "K"],
            alpha = input$growth[, "Alpha"]
        )
    })
    output$plot <- renderPlot({
        # Check if the simulation_result has been executed
        if (!is.null(simulation_result())) {
            # Create a plot of the simulation results
            log_plot(simulation_result()[[1]], type = input$display)
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