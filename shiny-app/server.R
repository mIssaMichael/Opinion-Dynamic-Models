#
# Author: Michael Issa
# Date: 11/30/2023
#

# server.R

library(shiny)
library(ggplot2)
library(patchwork)

# Set global parameters
R <- 15 # Number of rounds of the simulation
N <- 100 # Number of agents

# Function to simulate Hegselmann-Krause dynamics for one round
simulateRound <- function(currentOpinions, ConfidenceInterval, interactionTerm, RightBias) {
  ConfidenceIntervalLeft <- ConfidenceInterval * exp(-RightBias)
  ConfidenceIntervalRight <- ConfidenceInterval * exp(RightBias)
  newOpinions <- numeric(length(currentOpinions))
  
  for (agent in seq_along(currentOpinions)) {
    currentOpinion <- currentOpinions[agent]
    peersLeft <- which(currentOpinions > currentOpinion & currentOpinions <= currentOpinion + ConfidenceIntervalLeft)
    peersRight <- which(currentOpinions < currentOpinion & currentOpinions >= currentOpinion - ConfidenceIntervalRight)
    
    neighborhoodOpinions <- c(currentOpinion, currentOpinions[peersLeft], currentOpinions[peersRight])
    newOpinion <- interactionTerm * mean(neighborhoodOpinions)
    
    newOpinion <- newOpinion
    newOpinion <- pmax(-1, pmin(1, newOpinion)) # Ensure opinions are in the [-1, 1] interval
    newOpinions[agent] <- newOpinion
  }
  
  return(newOpinions)
}

# Function to run simulations and create plots
runSimulationAndPlot <- function(title, numRounds, ConfidenceInterval, interactionTerm, RightBias, numAgents) {
  # Initialize matrix to store opinions over rounds
  opinionsHistory <- matrix(data = NA, nrow = numRounds, ncol = numAgents + 1)
  opinionsHistory[1, ] <- seq(from = 0, to = 1, by = 1 / numAgents) # Adjust initial opinions
  
  # Run simulation for numRounds rounds
  for (round in 2:numRounds) {
    currentOpinions <- opinionsHistory[round - 1, ]
    newOpinions <- simulateRound(currentOpinions, ConfidenceInterval, interactionTerm, RightBias)
    
    # Ensure the number of agents remains constant
    opinionsHistory[round, 1:(numAgents + 1)] <- newOpinions[1:(numAgents + 1)]
  }
  
  # Reshape data for plotting
  plotData <- data.frame(
    Round = rep(1:numRounds, each = (numAgents + 1)),
    Opinion = as.vector(t(opinionsHistory)),
    Agent = as.factor(rep(1:(numAgents + 1), times = numRounds))
  )
  
  # Create the plot
  plot <- ggplot(plotData, aes(x = Round, y = Opinion, group = Agent, colour = Agent)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = rainbow(numAgents + 1)) +
    ggtitle(title) +
    labs(x = "Iterations", y = "Opinions") +
    theme_minimal() +
    theme(legend.position = "none") # Turn off legend
  
  return(plot)
}


# Define the server
function(input, output) {
  output$opinion_plot <- renderPlot({
    numRounds <- input$num_rounds
    ConfidenceInterval <- input$confidence_interval
    interactionTerm <- input$interaction_term
    RightBias <- input$right_bias
    numAgents <- input$num_agents
    
    title <- paste("Opinion Dynamics Simulation with",
                   "Rounds =", numRounds,
                   "CI =", ConfidenceInterval,
                   "Interaction Term =", interactionTerm,
                   "Right Bias =", RightBias,
                   "Agents =", numAgents)
    plot <- runSimulationAndPlot(title, numRounds, ConfidenceInterval, interactionTerm, RightBias, numAgents)
    print(plot)
  })
}
