#
# Author: Michael Issa
# Date: 11/30/2023
#

# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Set global parameters
R <- 15 # Number of rounds of the simulation
N <- 100 # Number of agents
ConfidenceInterval <- 0.1
RightBias <- 0.0
ConfidenceIntervalLeft <- ConfidenceInterval * exp(-RightBias)
ConfidenceIntervalRight <- ConfidenceInterval * exp(RightBias)
sigma <- 0.05 # Noise SD
alpha_linear <- 0.1 # Linear weight term
alpha_exp <- 0.05 # Exponential weight term
alpha_nonlinear <- 0.1 # Non-linear weight term
beta <- 0.5 #

# Function to simulate Hegselmann-Krause dynamics for one round
simulateRound <- function(currentOpinions, interactionTermFn, useBeta = TRUE) {
  # Initialize a vector for new opinions
  newOpinions <- numeric(N + 1)
  
  # Loop through each agent
  for (agent in 1:(N + 1)) {
    currentOpinion <- currentOpinions[agent]
    
    # Identify neighbors to the left and right within confidence intervals
    peersLeft <- which(currentOpinions > currentOpinion & currentOpinions <= currentOpinion + ConfidenceIntervalLeft)
    peersRight <- which(currentOpinions < currentOpinion & currentOpinions >= currentOpinion - ConfidenceIntervalRight)
    
    # Combine current agent's opinion with neighbors' opinions
    neighborhoodOpinions <- c(currentOpinion, currentOpinions[peersLeft], currentOpinions[peersRight])
    
    # Calculate interaction term using specified function
    interactionTerm <- interactionTermFn(neighborhoodOpinions, useBeta)
    
    # Update agent's opinion based on the interaction term and noise
    newOpinion <- mean(neighborhoodOpinions) + interactionTerm + rnorm(1, mean = 0, sd = sigma)
    newOpinion <- pmax(-1, pmin(1, newOpinion))    
    newOpinions[agent] <- newOpinion
  }
  
  return(newOpinions)
}

# Interaction term functions
interaction_linear <- function(neighborhoodOpinions, useBeta) {
  # Calculate linear interaction term based on neighbors' opinions and beta
  if (useBeta) {
    return(alpha_linear * sum(beta * neighborhoodOpinions) / length(neighborhoodOpinions))
  } else {
    return(alpha_linear * sum(neighborhoodOpinions) / length(neighborhoodOpinions))
  }
}

interaction_noInteraction <- function(neighborhoodOpinions, useBeta) {
  # Calculate interaction term without considering specific interaction function
  if (useBeta) {
    return(mean(c(beta * neighborhoodOpinions, neighborhoodOpinions)))
  } else {
    return(mean(neighborhoodOpinions))
  }
}

interaction_exponential <- function(neighborhoodOpinions, useBeta) {
  # Calculate exponential interaction term based on neighbors' opinions and beta
  if (useBeta) {
    return(alpha_exp * sum(beta * exp(neighborhoodOpinions)) / length(neighborhoodOpinions))
  } else {
    return(alpha_exp * sum(exp(neighborhoodOpinions)) / length(neighborhoodOpinions))
  }
}

interaction_nonlinear <- function(neighborhoodOpinions, useBeta) {
  # Calculate nonlinear interaction term based on neighbors' opinions and beta
  if (useBeta) {
    return(alpha_nonlinear * sum(beta * sqrt(neighborhoodOpinions)) / length(neighborhoodOpinions))
  } else {
    return(alpha_nonlinear * sum(sqrt(neighborhoodOpinions)) / length(neighborhoodOpinions))
  }
}

# Function to run simulations and create plots
runSimulationAndPlot <- function(interactionTermFn, title, useBeta = TRUE) {
  # Initialize a matrix to store opinions over rounds
  opinionsHistory <- matrix(data = NA, nrow = R, ncol = N + 1)
  opinionsHistory[1, ] <- seq(from = -1, to = 1, by = 2 / N)
  
  # Loop through rounds and simulate opinions
  for (round in 2:R) {
    opinionsHistory[round, ] <- simulateRound(opinionsHistory[round - 1, ], interactionTermFn, useBeta)
  }
  
  # Create a plot for the simulation results
  plot <- ggplot(data.frame(Round = rep(1:round, each = (N + 1)), 
                            Opinion = as.vector(t(opinionsHistory[1:round, ])), 
                            Agent = as.factor(rep(1:(N + 1), times = round))), 
                 aes(x = Round, y = Opinion, group = Agent, colour = Agent)) +
    geom_line(size = 0.5) +
    scale_color_manual(values = rainbow(N + 1)) +
    ggtitle(title) +
    labs(x = "Iterations", y = "Opinions") +
    theme_minimal()
  
  return(plot)
}

# Run simulations and create plots for different interaction terms and beta settings
plot_linear_withBeta <- runSimulationAndPlot(interactionTermFn = interaction_linear, title = "Linear Interaction Term with Beta")
plot_noInteraction_withBeta <- runSimulationAndPlot(interactionTermFn = interaction_noInteraction, title = "No Interaction Term with Beta")
plot_exponential_withBeta <- runSimulationAndPlot(interactionTermFn = interaction_exponential, title = "Exponential Interaction Term with Beta")
plot_nonlinear_withBeta <- runSimulationAndPlot(interactionTermFn = interaction_nonlinear, title = "Nonlinear Interaction Term with Beta")

plot_linear_noBeta <- runSimulationAndPlot(interactionTermFn = interaction_linear, title = "Linear Interaction Term without Beta", useBeta = FALSE)
plot_noInteraction_noBeta <- runSimulationAndPlot(interactionTermFn = interaction_noInteraction, title = "No Interaction Term without Beta", useBeta = FALSE)
plot_exponential_noBeta <- runSimulationAndPlot(interactionTermFn = interaction_exponential, title = "Exponential Interaction Term without Beta", useBeta = FALSE)
plot_nonlinear_noBeta <- runSimulationAndPlot(interactionTermFn = interaction_nonlinear, title = "Nonlinear Interaction Term without Beta", useBeta = FALSE)

# Combine the main plot with the interaction function plot
combined_plot <- grid.arrange(
  arrangeGrob(
    plot_noInteraction_withBeta, plot_noInteraction_noBeta,
    plot_linear_withBeta, plot_linear_noBeta,
    plot_exponential_withBeta, plot_exponential_noBeta,
    plot_nonlinear_withBeta, plot_nonlinear_noBeta,
    ncol = 4
  )
)

# Show the combined plot
print(combined_plot)

