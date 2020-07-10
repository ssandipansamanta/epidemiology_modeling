rm(list=ls(all=TRUE))
gc()

path = 'C:/Users/samantas/Documents/MyLearning/Epidemiology_modeling/project'
output_path = paste0(path, '/Output/')


seird_generator <- odin::odin({
  
  #' @Transition_Equations: Core equations for transitions between compartments:
  
  update(S) <- S + nBirth - nSE - nSQ + nQS - nSND
  update(E) <- E + nSE - nEI - nEND
  update(I) <- I + nEI - nIQ - nIR - nIH - nID - nIND
  update(Q) <- Q + nSQ + nIQ - nQR - nQS - nQD - nQND
  update(H) <- H + nIH - nHR - nHD - nHND
  update(R) <- R + nIR + nQR + nHR - nRND
  update(D) <- D + nID + nQD + nHD - nDBu 
  
  #' @Transition_Probabilities: Individual probabilities of transition:
  
  pSE  <- 1 - exp(-(beta1 * I + beta2 * H + beta3 * D) / N)
  pSQ  <- 1 - exp(-alpha1)
  pEI  <- 1 - exp(-gamma1)
  pIQ  <- 1 - exp(-alpha2)
  pIR  <- 1 - exp(-r1)
  pIH  <- 1 - exp(-h1)
  pID  <- 1 - exp(-d1)
  pQR  <- 1 - exp(-r2)
  pQS  <- 1 - exp(-l1)
  pQD  <- 1 - exp(-d2)
  pHR  <- 1 - exp(-r3)
  pHD  <- 1 - exp(-d3)
  pmu  <- 1 - exp(-mu)
  pbu  <- 1 - exp(-b1)
  
  #' @Drawing_Samples: Draws from binomial/Poisson distributions for numbers changing between compartments:
  
  nSE <- rbinom(S, pSE)
  nSQ <- rbinom(S, pSQ)
  nEI <- rbinom(E, pEI)
  nIQ <- rbinom(I, pIQ)
  nIR <- rbinom(I, pIR)
  nIH <- rbinom(I, pIH)
  nID <- rbinom(I, pID)
  nQR <- rbinom(Q, pQR)
  nQS <- rbinom(Q, pQS)
  nQD <- rbinom(Q, pQD)
  nHR <- rbinom(H, pHR)
  nHD <- rbinom(H, pHD)
  
  nSND <- rbinom(S, pmu)
  nEND <- rbinom(E, pmu)
  nIND <- rbinom(I, pmu)
  nQND <- rbinom(Q, pmu)
  nHND <- rbinom(H, pmu)
  nRND <- rbinom(R, pmu)
  nDBu <- rbinom(D, pbu)
  
  nBirth  <- rpois(lambda)
  
  N <- 1375980000
  
  #' @Initial_States
  initial(S) <- S_initial_value
  initial(E) <- E_initial_value
  initial(I) <- I_initial_value
  initial(Q) <- Q_initial_value
  initial(H) <- H_initial_value
  initial(R) <- R_initial_value
  initial(D) <- D_initial_value
  
  #' @parameters_Initialvalues - default in parentheses:
  S_initial_value <- user(1353890331)
  E_initial_value <- user(100)
  I_initial_value <- user(81)
  Q_initial_value <- user(10)
  H_initial_value <- user(7)
  R_initial_value <- user(9)
  D_initial_value <- user(2)
  
  
  lambda  <- user(1.11754218134554)
  beta1   <- user(0.551673049254618)
  beta2   <- user(0.0233876939989058)
  beta3   <- user(0.625131799387952)
  alpha1  <- user(4.4154696721785E-06)
  alpha2  <- user(0.0115667407044624)
  gamma1  <- user(0.2257890639763671)
  r1      <- user(0.0228005792053569)
  r2      <- user(8.26090265620624E-06)
  r3      <- user(0.371321586930202)
  d1      <- user(0.00960087708935922)
  d2      <- user(0.00175163503548737)
  d3      <- user(0.537274111027604)
  b1      <- user(0.510538433439926)
  h1      <- user(0.00308863244054197)
  l1      <- user(0.457429079552529)
  mu      <- user(0.0235627671386538)

  # lambda <-  user(1.12716195470781);	beta1 <-  user(0.599412888056549);	beta2 <-  user(0.671223526202388);	beta3 <-  user(0.0257791033900429);	alpha1 <-  user(1.67800924080207E-06);	alpha2 <-  user(0.0306467536266609);	gamma1 <-  user(0.205403628049027);	r1 <-  user(0.0228005792053569);	r2 <-  user(0.000399488982911961);	r3 <-  user(0.288692026793548);	d1 <-  user(0.00948111776918047);	d2 <-  user(0.000973999169501765);	d3 <-  user(0.289149552404002);	b1 <-  user(0.352604613157015);	h1 <-  user(0.00254891022923103);	l1 <-  user(0.455842525657497);	mu <-  user(0.0235849858531086);
  
}, verbose = FALSE)

## x: instance of odin model
## t: time steps
## n: number of replicates
run_model <- function(x, t = 0:100, n = 1, ...) {
  res <- x$run(t, replicate = n, ...)
  res <- x$transform_variables(res)
  res <- cbind.data.frame(t = res[[1]], res[-1])
  attr(res, "n_compartments") <- length(x$names) - 1
  attr(res, "n_replicates") <- n
  attr(res, "compartments") <- x$names[-1]
  class(res) <- c("pretty_odin", class(res))
  res
}

no_data_points <- 300
no_valid_simulation <- 0
seed <- 1
all_output <- data.frame(t=0:no_data_points)
while (no_valid_simulation <= 100){
  set.seed(seed)
  seird <- seird_generator()
  model_output <- run_model(seird, t = 0:no_data_points, n = 100)
  model_output_with_full_data <- model_output[, colSums(is.na(model_output)) == 0]
  if(!is.null(ncol(model_output_with_full_data))){
    edited_col_names <- stringr::str_replace(names(model_output_with_full_data)[!names(model_output_with_full_data) %in% 't'], "[.]", "_")
    final_col_names <- c("t",paste0(edited_col_names,"_",seed))
    colnames(model_output_with_full_data) <- final_col_names
    data.table::fwrite(x=model_output_with_full_data,file=paste0(output_path, '/simulation_result/Stochastic_output_with_seed_',seed,'.csv'),row.names = FALSE)
    all_output <- merge(x=all_output,y=model_output_with_full_data,by="t")
    no_valid_simulation <- (ncol(all_output) - 1)/7
    cat("\nValid Simuation:",no_valid_simulation,"; # Iteration:", seed-1)
    
  }
  seed = seed + 1
  gc()
}

data.table::fwrite(x=all_output,file=paste0(output_path, '/simulation_result/all_Stochastic_output.csv'),row.names = FALSE)
all_infection_cols <- stringr::str_subset(names(all_output),"I_")
all_recovered_cols <- stringr::str_subset(names(all_output),"R_")
all_deaths_cols <- stringr::str_subset(names(all_output),"D_")
all_hospitalization_cols <- stringr::str_subset(names(all_output),"H_")

chart_output <- all_output[,c("t", all_infection_cols, all_recovered_cols, all_deaths_cols, all_hospitalization_cols)]
data.table::fwrite(x=chart_output,file=paste0(output_path, '/simulation_result/Stochastic_chart.csv'),row.names = FALSE)


chart_output <- data.table::fread(input = paste0(output_path, '/simulation_result/Stochastic_chart.csv'),header = TRUE)
transpose_chart_output <- data.table::melt.data.table(data.table::setDT(chart_output),  id.vars = 't', variable.name = 'series')
summary_chart_output <- transpose_chart_output[,.(max_value = max(value)), keyby = series]

order_output <- summary_chart_output[order(max_value,decreasing=TRUE),]
subset_chart_data <- order_output[1:80,'series']


goodHosp <- c(stringr::str_replace_all(subset_chart_data$series,"[I]","I"),
              stringr::str_replace_all(subset_chart_data$series,"[I]","R"),
              stringr::str_replace_all(subset_chart_data$series,"[I]","H"),
              stringr::str_replace_all(subset_chart_data$series,"[I]","D"))
goodHosp_match <- paste0(goodHosp, collapse = "|")

library(data.table)
plot_output <- transpose_chart_output[series %like% goodHosp_match]
data.table::fwrite(x=plot_output,file=paste0(output_path, '/simulation_result/Stochastic_chart_v01.csv'),row.names = FALSE)

library(ggplot2)
ggplot(plot_output, aes(t,value)) + geom_line(aes(colour = series),show.legend = FALSE)
