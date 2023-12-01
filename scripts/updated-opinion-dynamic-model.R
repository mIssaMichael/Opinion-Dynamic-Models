#
# Author: Michael Issa
# Date: 11/30/2023
#


library(ggplot2)
library(patchwork)

# Set global parameters
R <- 15 # Number of rounds of the simulation
N <- 100 # Number of agents
interactionTerm <- 1 # alpha weight term for how much individuals give weight to their neighbors

# Function to simulate Hegselmann-Krause dynamics for one round
simulateRound <- function(currentOpinions, ConfidenceInterval) {
  RightBias <- 0.0 # Bias towards people
  ConfidenceIntervalLeft <- ConfidenceInterval * exp(-RightBias)
  ConfidenceIntervalRight <- ConfidenceInterval * exp(RightBias)
  newOpinions <- numeric(N + 1)
  
  for (agent in 1:(N + 1)) {
    currentOpinion <- currentOpinions[agent]
    peersLeft <- which(currentOpinions > currentOpinion & currentOpinions <= currentOpinion + ConfidenceIntervalLeft)
    peersRight <- which(currentOpinions < currentOpinion & currentOpinions >= currentOpinion - ConfidenceIntervalRight)
    
    neighborhoodOpinions <- c(currentOpinion, currentOpinions[peersLeft], currentOpinions[peersRight])
    newOpinion <- interactionTerm * mean(neighborhoodOpinions)
    
    newOpinion <- newOpinion
    newOpinion <- pmax(0, pmin(1, newOpinion)) # Ensure opinions are in the [0, 1] interval
    newOpinions[agent] <- newOpinion
  }
  
  return(newOpinions)
}

# Function to run simulations and create plots
runSimulationAndPlot <- function(title, ConfidenceInterval) {
  # Initialize matrix to store opinions over rounds
  opinionsHistory <- matrix(data = NA, nrow = R, ncol = N + 1)
  opinionsHistory[1, ] <- seq(from = 0, to = 1, by = 1 / N) # Adjust initial opinions
  
  # Run simulation for R rounds
  for (round in 2:R) {
    opinionsHistory[round, ] <- simulateRound(opinionsHistory[round - 1, ], ConfidenceInterval)
  }
  
  # Reshape data for plotting
  plotData <- data.frame(
    Round = rep(1:R, each = (N + 1)),
    Opinion = as.vector(t(opinionsHistory)),
    Agent = as.factor(rep(1:(N + 1), times = R))
  )
  
  # Create the plot
  plot <- ggplot(plotData, aes(x = Round, y = Opinion, group = Agent, colour = Agent)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = rainbow(N + 1)) +
    ggtitle(title) +
    labs(x = "Iterations", y = "Opinions") +
    theme_minimal() +
    theme(legend.position = "none") # Turn off legend
  
  return(plot)
}

# Run simulation and create plots
plot1 <- runSimulationAndPlot("Opinion Dynamics Simulation with CI=.05", 0.05)
plot2 <- runSimulationAndPlot("Opinion Dynamics Simulation with CI=.15", 0.15)
plot3 <- runSimulationAndPlot("Opinion Dynamics Simulation with CI=.25", 0.25)

# Combine the plots using patchwork
combined_plot <- plot1 + plot2 + plot3

# Display the plot
print(combined_plot)
