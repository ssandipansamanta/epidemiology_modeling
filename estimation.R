path = 'C:/Users/samantas/Documents/MyLearning/Epidemiology_modeling/project'
file_name = 'input_data_MH'

input_path = paste0(path, '/Input/')
output_path = paste0(path, '/Output/')

input <- read.csv(paste0(input_path,file_name,'.csv'),header = TRUE)  
N = 121924973


#' @Required_Libraries:
#' 

library(deSolve)
library(rgenoud)
library(parallel)

 
#' @Model_Formulation:
#' 

SEIRD_model <- function(day,state, params){
  with(as.list(c(state, params)), {
    
    beta = params[1]; delta = params[2]; alpha = params[3]; gamma=params[4]; rho=params[5]
    
    dS = -beta * S * I / N
    dE = beta * S * I / N - delta * E
    dI = delta * E - (1 - alpha) * gamma * I - alpha * rho * I
    dR = (1 - alpha) * gamma * I
    dD = alpha * rho * I

    return(list(c(dS, dE, dI, dR, dD)))
  })
}


#' @Objective_Function:
#' 

obj_function <- function(ga_param){
  de_output   <- deSolve::ode(y=current_state, times=time_duration, func = SEIRD_model, parms = ga_param, method = ode_method)
  overall_error <- sum((input$infected - de_output[,4])^2 + (input$recovered - de_output[,5])^2 + (input$deceased - de_output[,6])^2)
  # overall_error = sum(abs(input$infected - de_output[,4]) + abs(input$recovered - de_output[,5]) + abs(input$deceased - de_output[,6]))
  
  return(overall_error)
}

#' @Initial_parameters:
#' 


initial_param  <- 
  c(
    beta  = 0.2,     # infected person infects 1 other person per day
    delta = 1/5,     # incubation period of 5 days
    alpha = 0.2,     # death Rate
    gamma = 1/4,    # infection last 10 days
    rho   = 1/5     # 15 days from infection to death 
  ) 

initial_state  <- c(S = N - 14, E = 0, I = 14, R = 0, D = 0)
times <- seq(1, nrow(input), by = 1)
current_state = initial_state; 
time_duration = times
ode_method = "lsode"


#' @Genetic_Algorithm_parameters:
#' 

pop_size <- 100
max_generations <- 200
PM1=0.3; PM2=0.2; PM3=0.2; PM4=0.2; PM5=0.8; PM6=0.8; PM7=0.3; PM8=0.8; PM9=0.8

lb <- c(0.1,0.02,0.01,0.1,0.15); ub <- c(0.5,0.5,0.25,0.5,0.5)
limit <- cbind(lb,ub)


#' @Optimization
#' 

Cluster_Export_DataSets <- c("input",'SEIRD_model',"ode_method", "current_state", "time_duration","N")
cl <- parallel::makeCluster(detectCores(), type="SOCK")
parallel::clusterExport(cl, Cluster_Export_DataSets)

opt_output <-rgenoud::genoud(fn=obj_function, 
                              nvars=length(initial_param), 
                              max=FALSE,
                              pop.size=pop_size,
                              max.generations=max_generations,
                              wait.generations=2,
                              hard.generation.limit=TRUE,
                              starting.values=initial_param,
                              MemoryMatrix=TRUE,
                              Domains=limit, default.domains=NULL, solution.tolerance=1E-2,
                              gr=NULL, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE,
                              BFGS=TRUE, data.type.int=FALSE, hessian=FALSE,
                              unif.seed=123, int.seed=456,
                              print.level=1,
                              P1=PM1, P2=PM2, P3=PM3, P4=PM4, P5=PM5, P6=PM6, P7=PM7, P8=PM8, P9=PM9,
                              P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp=NULL,
                              #control=list(abstol = 10000, reltol =10000),
                              optim.method="L-BFGS-B",#"BFGS",L-BFGS-B
                              transform=FALSE, debug=FALSE, cluster=cl, balance=TRUE)

stopCluster(cl)

#' @Final_Parameters:
#' 

final_params <- opt_output$par

