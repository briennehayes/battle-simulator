# CS 558 Spring '20 Homework 3
# Brienne Hayes

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# define ui layout
ui <- fluidPage(
  
  titlePanel("Battle Simulation using Lanchester's Law for Aimed Fire"),
  
  sidebarLayout(
    sidebarPanel(
      
      # separate "about" dialogue
      actionButton("about", "About"),
      
      # Simulation value inputs
      numericInput("orcsInit",
                   "Initial army population for Orcs:",
                   value = 1000,
                   min = 0,
                   step = 100),
      sliderInput("orcsCoeff",
                  "Lethality coefficient for Orcs:",
                  min = 0, max = 1,
                  value = 0.8),
      numericInput("elvesInit",
                   "Initial army population for Elves:",
                   value = 800,
                   min = 0,
                   step = 100),
      sliderInput("elvesCoeff",
                  "Lethality coefficient for Elves:",
                  min = 0, max = 1,
                  value = 0.9),
      numericInput("endTime",
                   "Simulation end time:",
                   value = 2,
                   min = 0.1,
                   step = 0.5),
      numericInput("stepSize",
                   "Time step size:",
                   value = 0.1,
                   min = 0.01,
                   step = 0.1),
      checkboxInput("showReinforcements", "Allow Reinforcements"),
      
      # additional inputs for reinforcements
      conditionalPanel(
        condition = "input.showReinforcements == true",
        h4("Orc Reinforcements"),
        radioButtons("numReinforcementsOrcs",
                     "Number of reinforcement events:",
                     c("0" = 0,
                       "1" = 1,
                       "2" = 2,
                       "3" = 3), inline = TRUE),
        sliderInput("whenReinforcementOrcs",
                    "% of initial army population to call for reinforcements:",
                    value = 0.5,
                    min = 0.1, max = 0.8,
                    step = 0.01
        ),
        sliderInput("reinforcementSizeOrcs",
                    "Reinforcement size as % of initial army population:",
                    value = 0.25,
                    min = 0.1, max = 0.5,
                    step = 0.01),
        uiOutput("reinforcementCoeffsOrcs"),
        h4("Elf Reinforcements"),
        radioButtons("numReinforcementsElves",
                     "Number of reinforcement events:",
                     c("0" = 0,
                       "1" = 1,
                       "2" = 2,
                       "3" = 3), inline = TRUE),
        sliderInput("whenReinforcementElves",
                    "% of initial army population to call for reinforcements:",
                    value = 0.5,
                    min = 0.1, max = 0.8,
                    step = 0.01
        ),
        sliderInput("reinforcementSizeElves",
                    "Reinforcement size as % of initial army population:",
                    value = 0.25,
                    min = 0.1, max = 0.5,
                    step = 0.01),
        uiOutput("reinforcementCoeffsElves")
      )
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("simPlot")),
                  tabPanel("Table", dataTableOutput("battleHistory")))
    )
  )
)

# define server logic, perform calculations and build visualizations
server <- function(input, output) {
  
  observeEvent(input$about, {
    showModal(modalDialog(
      title = "About",
      p(h4("The Lanchester Simulation"),
        "This app simulates a battle between an army of orcs and an army of elves.",
        "The simulation follows Lanchester's Law for Aimed Fire, in which the size of each army decreases proportionally to the size of the enemy army.",
        "Each army's lethality coefficient represents the probability that one unit hits one enemy unit at each time step.",
        "Solutions to the Lanchester system are approximated using Euler's method of numerical integration.",
        "The default settings reflect the orc's numerical superiority and the battle prowess of the elves, but these values can be freely reconfigured."),
      p(h4("Reinforcements"),
        "Checking the Allow Reinforcements box enables an expanded interface to control battle reinforcements.",
        "Each army can have up to three waves of reinforcements. On each side, the three waves will consist of a fixed number of units",
        "and will occur when the current army population drops below a threshold. Both of these values are specified as proportions",
        "of the initial army size. Each of the three waves can, however, have different lethality coefficients.",
        "The army's new lethality coefficient following a reinforcement is the soldier-wise average of the new and remaining soldiers' coefficients."),
      p(h4("Known Issues"),
        "An error message may display for a fraction of a second when switching between numbers of reinforcement waves. This behavior",
        "should in no way impact the simulation output."),
      p(h4("Accreditation"),
        "This application was made by Brienne Hayes for Raymond Madachy's Computer Simulation class in Spring 2020.",
        "It is written in Shiny and dependent on tidyverse packages.",
        "Questions, issues, and suggestions for additional features should be sent to bhayes@sdsu.edu.")
      
    ))
  })
  
  # Approximates numerical integration for one time step
  euler <- function(x, dx, dt){return(x + (dx * dt))}
  
  # Defines derivatives under the Lanchester model
  deriv <- function(coeff, val){return(-coeff * val)}
  
  # Uses Euler's method to simulate a battle under the Lanchester model
  simulateBattle <- function(x0, y0, alpha, beta, endTime, stepSize, xBackup, yBackup){
    t_vals <- seq(from = 0, to = endTime, by = stepSize) 
    dt_vals <- diff(t_vals) # vector where all values are dt
    x_vals <- c(x0, rep(0, times = length(dt_vals) - 1))
    y_vals <- c(y0, rep(0, times = length(dt_vals) - 1))
    i <- 1
    
    # 
    # backupLog <- rep("", times = length(dt_vals)) # this will allow the data table to display when a backup happens
    
    # unpack backup variables
    if(!is.na(xBackup) && !is.na(yBackup)){
      backups <- TRUE # flag to check whether to perform backup computations
      
      xNumBackups <- xBackup[[1]]
      xBackupCounter <- 0
      xBackupThreshold <- xBackup[[2]] * x0 # to get threshold, multiply given proportion by original pop size
      xBackupSize <- xBackup[[3]] * x0 # to get backup size, multiply given proportion by original pop size
      
      # backup coefficients (may be NULL)
      xBackupCoeffs <- list(xBackup[[4]], xBackup[[5]], xBackup[[6]])
      
      yNumBackups <- yBackup[[1]]
      yBackupCounter <- 0
      yBackupThreshold <- yBackup[[2]] * y0 # to get threshold, multiply given proportion by original pop size
      yBackupSize <- yBackup[[3]] * y0 # to get backup size, multiply given proportion by original pop size
      
      # backup coefficients (may be NULL)
      yBackupCoeffs <- list(yBackup[[4]], yBackup[[5]], yBackup[[6]])
    }
    else{
      backups <- FALSE
    }
    
    # these will be overwritten in backup events
    xCoeff <- alpha
    yCoeff <- beta
    
    for(dt in dt_vals){
      i <- i + 1
      
      # calculate population change at time step
      dx <- deriv(yCoeff, y_vals[i - 1])
      dy <- deriv(xCoeff, x_vals[i - 1])
      
      # apply Euler's method at time step
      x_next <- euler(x_vals[i - 1], dx, dt)
      y_next <- euler(y_vals[i - 1], dy, dt)
      
      # population values should not go below zero
      if(x_next < 0){x_next <- 0}
      if(y_next < 0){y_next <- 0}
      
      # if population survives, consider reinforcement events
      # TODO this can probably be generalized so the code doesn't repeat
      if(backups){
        
        if (x_next < xBackupThreshold && xBackupCounter < xNumBackups){
          xBackupCounter <- xBackupCounter + 1 # increment backup counter
          
          xNewPop <- x_next + xBackupSize
          
          # compute new coefficient
          xCoeff <- ((x_next * xCoeff) + (xBackupSize * xBackupCoeffs[[xBackupCounter]])) / xNewPop
          
          x_vals[i] <- xNewPop
        }
        else{
          x_vals[i] <- x_next
        }
        
        if (y_next < yBackupThreshold && yBackupCounter < yNumBackups){
          yBackupCounter <- yBackupCounter + 1 # increment backup counter
          
          yNewPop <- y_next + yBackupSize
          
          # compute new coefficient
          yCoeff <- ((y_next * yCoeff) + (yBackupSize * yBackupCoeffs[[yBackupCounter]])) / yNewPop
          
          y_vals[i] <- yNewPop
        }
        else{
          y_vals[i] <- y_next
        }
      }
      else{
        x_vals[i] <- x_next
        y_vals[i] <- y_next
      }
    }
    
    data = data.frame(t_vals, x_vals, y_vals) %>%
      gather(key = "var", val = "pop", x_vals:y_vals)
    
    return(data)
  }
  
  results <- reactive({
    req(input$orcsInit, input$elvesInit) # ensures inputs are set
    
    req(input$endTime > input$stepSize & input$stepSize > 0) # ensures reasonable time values are given
    
    # packing and unpacking variables like this is probably very inefficient, but I'll find a better solution later
    if (input$showReinforcements){
      req(input$numReinforcementsOrcs)
      req(input$numReinforcementsElves)
      
      orcsBackup <- list(input$numReinforcementsOrcs, input$whenReinforcementOrcs,
                         input$reinforcementSizeOrcs, input$reinforcementCoeff1Orcs,
                         input$reinforcementCoeff2Orcs, input$reinforcementCoeff3Orcs)
      elvesBackup <- list(input$numReinforcementsElves, input$whenReinforcementElves,
                          input$reinforcementSizeElves, input$reinforcementCoeff1Elves,
                          input$reinforcementCoeff2Elves, input$reinforcementCoeff3Elves)
    }
    else{
      orcsBackup <- NA
      elvesBackup <- NA
    }
    
    
    simulateBattle(input$orcsInit, input$elvesInit, input$orcsCoeff, input$elvesCoeff,
                   input$endTime, input$stepSize, orcsBackup, elvesBackup)
  })
  
  output$simPlot <- renderPlot({
    ggplot(results(), aes(x = t_vals, y = pop)) +
      geom_line(aes(color = var), size = 2) +
      scale_colour_discrete(name  ="Army",
                            breaks=c("x_vals", "y_vals"),
                            labels=c("Orcs", "Elves")) +
      labs(title = "Simulation Results", x = "Time", y = "Army Population")
  })
  
  output$battleHistory <- renderDataTable({
    rename(spread(results(), key = "var", val = "pop"), Time = t_vals, Orcs = x_vals, Elves = y_vals)
  })
  
  # number of sliders to show for Orc reinforcements
  numInputsOrcs <- reactive({
    input$numReinforcementsOrcs
  })
  
  output$reinforcementCoeffsOrcs <- renderUI({
    req(input$showReinforcements) # ensure reinforcement panel is active
    numReinforcements <- numInputsOrcs() # get how many reinforcements are called for
    
    # generate one slider for each reinforcement wave
    output <- tagList()
    if (numReinforcements > 0){
      for(i in seq_along(1:numReinforcements)){
        output[[i]] <- tagList()
        output[[i]][[1]] <- sliderInput(paste0("reinforcementCoeff", i, "Orcs"),
                                        paste("Reinforcement", i, "lethality coefficient:"),
                                        min = 0, max = 1,
                                        value = 0.8)
      }
    }
    output
  })
  
  # number of sliders to show for Elf reinforcements
  numInputsElves <- reactive({
    input$numReinforcementsElves
  })
  
  output$reinforcementCoeffsElves <- renderUI({
    req(input$showReinforcements) # ensure reinforcement panel is active
    numReinforcements <- numInputsElves() # get how many reinforcements are called for
    
    # generate one slider for each reinforcement wave
    output <- tagList()
    if (numReinforcements > 0){
      for(i in seq_along(1:numReinforcements)){
        output[[i]] <- tagList()
        output[[i]][[1]] <- sliderInput(paste0("reinforcementCoeff", i, "Elves"),
                                        paste("Reinforcement", i, "lethality coefficient:"),
                                        min = 0, max = 1,
                                        value = 0.8)
      }
    }
    output
  })
  
}

# run the application 
shinyApp(ui = ui, server = server)



