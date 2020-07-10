rm(list=ls(all=TRUE))
gc()

path = 'Path of root Directory'
file_name = 'Input file name'
initial_param_file_name = 'initial params bounds'

N = 8000 #1353890423

#' @Required_Libraries:
#' 

library(deSolve)
library(rgenoud)
library(parallel)
library(plotly)
library(data.table)
library(nsga2R)

#' @Inputs:
#' 

input_path = paste0(path, '/Input/')
output_path = paste0(path, '/Output/')

input <- data.table::fread(paste0(input_path,file_name,'.csv'),header = TRUE)

input[,cum_infected:=cumsum(infected)]
input[,cum_recovered:=cumsum(recovered)]
input[,cum_deceased:=cumsum(deceased)]

#' @Model_Formulation:
#' 

SEIRD_model <- function(day,state, params){
  with(as.list(c(state, params)), {
    
    lamda=params[1];
    beta1 = params[2];beta2 = params[3]; beta3 = params[4];  
    alpha1 = params[5]; alpha2 = params[6]; 
    r1=params[8];r2=params[9]; r3=params[10];
    d1=params[11]; d2=params[12]; d3=params[13];
    gama1=params[7];b1=params[14];h1=params[15];l=params[16];mu=params[17]
    
    dS = lamda +l*Q-(beta1 * S * I / N)-(beta2 * S * D / N)-(beta3 * S * H / N)-(alpha1 * S)-(mu * S)
    dE = (beta1 * S * I / N)+(beta2 * S * D / N)+(beta3 * S * H / N)-(gama1 * E)-(mu * E)
    dI = (gama1 * E) - I * (alpha2 + r1 + h1 + d1 + mu)
    dQ = (alpha1 * S)+(alpha2 * I)-(r2 * Q)-(d2 * Q)-(l * Q)-(mu * Q)
    dH = (h1 * I)-(r3 * H)-(d3 * H)-(mu * H)
    dR = (r1 * I)+(r2 * Q)+(r3 * H)-(mu * R)
    dD = (d1 * I)+(d2 * Q)+(d3 * H)-(b1 * D)
    
    return(list(c(dS, dE, dI, dQ, dH, dR, dD)))
  })
}


#' @Objective_Function: The objective function has single objctive value i.e. ESS or EAS
#' 

obj_function <- function(ga_param){
  de_output   <- as.data.frame(deSolve::ode(y=current_state, times=time_duration, func = SEIRD_model, parms = ga_param, method = ode_method))
  overall_error <- sum((input$cum_infected - de_output$I)^2 + (input$cum_recovered - de_output$R)^2 + (input$cum_deceased - de_output$D)^2)
  # overall_error = sum(abs(input$infected - de_output$I) + abs(input$recovered - de_output$R) + abs(input$deceased - de_output$D))
  return(1/overall_error)
}


moea_obj_function <- function(ga_param){
  de_output   <- as.data.frame(deSolve::ode(y=current_state, times=time_duration, func = SEIRD_model, parms = ga_param, method = ode_method))
  infected_error <- 1/sum(abs(input$cum_infected - de_output$I)) #1/sum((input$cum_infected - de_output$I)^2)
  recovered_error <- 1/sum(abs(input$cum_recovered - de_output$R))
  deceased_error <- 1/sum(abs(input$cum_deceased - de_output$D))
  return(c(-infected_error, -recovered_error, -deceased_error))
}


#' @Initial_parameters:
#' 


initial_param_input <- read.csv(paste0(input_path,initial_param_file_name,'.csv'),header = TRUE)  
initial_param <- initial_param_input$initial_param
initial_state  <- c(S = N - (input$infected[1]+input$recovered[1]+input$deceased[1]), 
                    E = 100, 
                    I = input$infected[1], 
                    Q = 10, 
                    H = (input$recovered[1] - input$deceased[1]), 
                    R = input$recovered[1], 
                    D = input$deceased[1])
times <- seq(1, nrow(input), by = 1)
current_state = initial_state; 
time_duration = times
ode_method = "lsoda"

# "lsode", "lsodes", "lsodar", "vode", "daspk",
# "euler", "rk4", "ode23", "ode45", "radau", 
# "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration"
# obj_function(ga_param = initial_param)

#' @Genetic_Algorithm_parameters:
#' 

pop_size <- 50*length(initial_param)
max_generations <- 50
PM1=0.9; PM2=0.3; PM3=0.5; PM4=0.8; PM5=0.8; PM6=0.8; PM7=0.3; PM8=0.8; PM9=0.5

lb <- initial_param_input$lower_limit
ub <- initial_param_input$upper_limit
limit <- cbind(lb,ub)


#' @GA_parameters:
#' 
 
# Cluster_Export_DataSets <- c("input",'SEIRD_model',"ode_method", "current_state", "time_duration","N")
# cl <- parallel::makeCluster(detectCores(), type="SOCK")
# parallel::clusterExport(cl, Cluster_Export_DataSets)
# 
# opt_output <-rgenoud::genoud(fn=obj_function,
#                              nvars=length(initial_param),
#                              max=TRUE,
#                              pop.size=pop_size,
#                              max.generations=max_generations,
#                              wait.generations=2,
#                              hard.generation.limit=TRUE,
#                              starting.values=initial_param,
#                              MemoryMatrix=TRUE,
#                              Domains=limit, default.domains=NULL, solution.tolerance=1E-10,
#                              gr=NULL, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE,
#                              BFGS=TRUE, data.type.int=FALSE, hessian=FALSE,
#                              unif.seed=123, int.seed=456,
#                              print.level=1,
#                              P1=PM1, P2=PM2, P3=PM3, P4=PM4, P5=PM5, P6=PM6, P7=PM7, P8=PM8, P9=PM9,
#                              P9mix=NULL, BFGSburnin=0, BFGSfn=NULL, BFGShelp=NULL,
#                              #control=list(abstol = 10000, reltol =10000),
#                              optim.method="L-BFGS-B",#"BFGS",L-BFGS-B
#                              transform=FALSE, debug=FALSE, cluster=cl, balance=TRUE)
# 
# stopCluster(cl)


#' @Final_Parameters:
#'

# final_params <- opt_output$par


#' @MOEA_parameters:
#' 

set.seed(123)
opt_output <- nsga2R::nsga2R(fn = moea_obj_function, 
                     varNo = length(initial_param), 
                     objDim = 3, 
                     generations = max_generations, 
                     popSize = pop_size,
                     lowerBounds = lb, 
                     upperBounds = ub,
                     tourSize = 2,
                     cprob = 0.8, XoverDistIdx = 5,
                     mprob = 0.3, MuDistIdx  = 10)

plot(opt_output$objectives)

pareto_rank <- as.data.table(opt_output$paretoFrontRank)
all_parameters <- as.data.table(cbind(opt_output$parameters,opt_output$paretoFrontRank))
all_objectives <- as.data.table(cbind(opt_output$objectives,opt_output$paretoFrontRank))


table(opt_output$paretoFrontRank)
all_possible_params <- all_parameters[all_parameters$V18==1,]


#' @Simulation:
#' 


# final_params <- initial_param #c( beta = final_params[1], delta = final_params[2], alpha = final_params[3], gamma=final_params[4], rho=final_params[5])
# days <- seq(0, nrow(input)-1, by = 1)
all_plot_output <- NULL
for(iter in (1:nrow(all_possible_params))){
  # iter <- 2
  final_params <- t(all_possible_params[iter,1:length(initial_param)])[,1]
  days <- seq(0, 300, by = 1)
  estimated_output <- as.data.table(deSolve::ode(y = initial_state, times = days, func = SEIRD_model, parms = final_params, method = ode_method))
  
  plot_output <- merge.data.table(estimated_output, input, by.x = "time", by.y = "day", all.x = TRUE)
  plot_output$generation <- iter
  all_plot_output <- as.data.table(rbind(all_plot_output,plot_output))
  
}



days <- seq(0, 300, by = 1)
estimated_output <- as.data.table(deSolve::ode(y = initial_state, times = days, func = SEIRD_model, parms = final_params, method = ode_method))
plot_output <- merge.data.table(estimated_output, input, by.x = "time", by.y = "day", all.x = TRUE)



fig <- plotly::plot_ly(plot_output, x = ~days)
# fig <- fig %>% add_trace(y = ~S, name = 'Susceptibles',mode = 'lines',line = list(color = 'rgba(65, 131, 215, 1)', width = 4))
# fig <- fig %>% add_trace(y = ~E, name = 'Exposed', mode = 'lines',line = list(color = 'rgba(247, 202, 24, 1)', width = 4))
fig <- fig %>% add_trace(y = ~I, name = 'Infecteds', mode = 'lines',line = list(color = 'rgba(30, 130, 76, 1)', width = 4))
# fig <- fig %>% add_trace(y = ~Q, name = 'Quarantined', mode = 'lines',line = list(color = 'rgba(242, 38, 19, 1)', width = 4))
# fig <- fig %>% add_trace(y = ~H, name = 'Hospitalized', mode = 'lines',line = list(color = 'rgba(165, 55, 253, 1)', width = 4))

fig <- fig %>% add_trace(y = ~R, name = 'Recovereds', mode = 'lines',line = list(color = 'rgba(118, 93, 105, 1)', width = 4))
fig <- fig %>% add_trace(y = ~D, name = 'Deaths', mode = 'lines',line = list(color = 'rgba(241, 130, 141,1)', width = 4))
fig <- fig %>% add_trace(y = ~cum_infected, name = 'actual infected', mode = 'markers',
                         marker = list(size = 10, color = 'rgba(30, 130, 76, 1)', line = list(color = 'rgba(152, 0, 0, .8)', width = 2)))
fig <- fig %>% add_trace(y = ~cum_recovered, name = 'actual recovered', mode = 'markers',
                         marker = list(size = 10, color = 'rgba(118, 93, 105, 1)', line = list(color = 'rgba(152, 0, 0, .8)', width = 2)))
fig <- fig %>% add_trace(y = ~cum_deceased, name = 'actual deaths', mode = 'markers',
                         marker = list(size = 10, color = 'rgba(241, 130, 141,1)', line = list(color = 'rgba(152, 0, 0, .8)', width = 2)))

fig

# round(head(plot_output[,c("E","I","cum_infected","Q","R","cum_recovered","H","D","cum_deceased")],n=20),0)
# round(head(plot_output[,c("I","cum_infected","R","cum_recovered","D","cum_deceased")],n=30),0)
# cbind(initial_param_input[,1],round(initial_param_input[,2:5],3))
# fwrite(x=plot_output,file=paste0(output_path, 'Fitted_Output_with_',max_generations,'_gen.csv'),row.names = FALSE)
# fwrite(x=initial_param_input,file=paste0(output_path, 'Params_estimates_with_',max_generations,'_gen.csv'),row.names = FALSE)
fwrite(x=all_plot_output,file=paste0(output_path, 'all_Output_with_500',max_generations,'_gen.csv'),row.names = FALSE)
fwrite(x=all_possible_params,file=paste0(output_path, 'all_possible_params_with_500',max_generations,'_gen.csv'),row.names = FALSE)
fwrite(x=all_objectives,file=paste0(output_path, 'all_objective_value_with_500',max_generations,'_gen.csv'),row.names = FALSE)
