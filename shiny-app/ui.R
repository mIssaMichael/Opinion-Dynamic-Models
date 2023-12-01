#
# Author: Michael Issa
# Date: 11/30/2023
#

# ui.R
library(shiny)

fluidPage(
  titlePanel("Hegselmann-Krause Opinion Dynamics Simulation"),
  
  withMathJax(),
  fluidRow(style = "background-color:#F2F2F2; margin-top: 30px; margin-bottom: 30px; padding: 10px",
           column(
             width = 6,
             # Introduction text:
             p(
               tags$b("Description."),
              "Agent-based computer simulations aid in studying various social phenomenon. In this case, the emergence of consensus, polarization and fragementation is the topic of study.
              A model used to study those phenomena is the Hegselmann-Krause Opinion Dynamics model, which uses fairly simple principles to aid in understanding."
             ),
             p("First, pick a number of agents \\(N\\) who are given the task of determining the real value of some parameter lying in the open unit interval \\((0,1)\\). 
               The initial opinions of the agents are uniformly distrubted over the interval. Agents are then put through various rounds of opinion formation modeled as a time step \\(t\\) where they 
               update their opinion for the the next round \\(x_i(t+1)\\) according to a weighted average of all their neighboring agents who don't differ in their opinion any more than \\(\\epsilon_i\\) (their confidence interval).
               The original HK dynamics are given by this equation:"),
             p("\\( \\displaystyle x_i(t+1)=|I(i,x(t))|^{-1} \\sum_{j \\in I(i,x(t))} x_j(t)\\),", align = "center"),
             p("where \\(I(i,x) = \\{1 \\leq j \\leq N: |x_i - x_j| \\leq \\epsilon_i \\} \\).")
           ),
           column(
             width = 6,
             # Introduction continued:
             p("However, the implementation below uses a sligly different formulation of the dynamics that makes explicit the weighting factor \\(\\alpha_i\\) (interaction term below). The weighting factor informs us how much 
               information about the real value of the parameter the agents are getting by way of their peers rather than some other source of information. The dyanmics become: "),
             p("\\( \\displaystyle x_i(t+1)=\\alpha_i|I(i,x(t))|^{-1} \\sum_{j \\in I(i,x(t))} x_j(t) + (1- \\alpha_i)\\tau \\),", align = "center"),
             p("where \\(\\tau\\) is the true value of the parameter, which agents are attempting to find by talking to each other, doing experiments, or exploiting some other information channel. In the original dynamics, the value of
               \\(\\alpha_i = 1\\) because agents are only receiving information from each other. The dyanmics outlined here include more of the real situation of the agents."),
             p("Overall, this model exemplifies two key features we like in computational models: simple dynamics which give rise to complex behavior and room for expanding and including missing details"),
             br(),
             br(),
             # Refernces:
             tags$b("Citations"),
             p("Douven, I., & Riegler, A. (2010). Extending the Hegselmann-Krause Model I. Log. J. IGPL, 18, 323-335."),
             p("Hegselmann, R., & Krause, U. (2002). Opinion dynamics and bounded confidence: models, analysis and simulation. J. Artif. Soc. Soc. Simul., 5.")
           )),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("confidence_interval", "Confidence Interval", min = 0, max = 1.0, value = 0.05, step = 0.01),
      sliderInput("interaction_term", "Interaction Term", min = 0, max = 1, value = 1, step = 0.1),
      sliderInput("right_bias", "Right Bias", min = -1, max = 1, value = 0, step = 0.1),
      sliderInput("num_rounds", "Number of Rounds", min = 1, max = 100, value = 15, step = 1),
      sliderInput("num_agents", "Number of Agents", min = 1, max = 100, value = 100, step = 1)
    ),
    mainPanel(
      plotOutput("opinion_plot")
    )
  )
)
