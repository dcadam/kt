#### BRANCHING PROCCESS SIMULATION

library(epi.branch.sim)
library(epitrix)
library(tidyverse)
library(epicontacts)
library(foreach)
library(doParallel)
source('src/branching-params.R')
source('src/util.R')


## Vary R and k inputs the epidemics. 
k <- c(0.1, 0.5, 1, 10, 100)
R <- 2

Rk_iteration <- crossing(R, k)
iter <- 300 # Each combination R & k repeats 300 times. 

### SIMULATION

ncpus <- parallel::detectCores()

cluster <- makeCluster(ncpus)
registerDoParallel(cluster)

set.seed(123)
out <- vector(mode = "list", length = 5)
for (i in 1:length(out)) {
  
  sims_out <- foreach(s = 1:iter, .packages=c('dplyr', 'epi.branch.sim')) %dopar% {
    
    R0 <- Rk_iteration$R[[i]]
    disp <- Rk_iteration$k[[i]]
    sec_infect_params <- list(type = "Hellewell", disp = disp)
    
    # Create initial simulation objects
    sim_params <- initialize_sim_params(
      R0, infect_dur, do_variable_trace, p_trace,
      p_trace_app, p_trace_app_comp, p_symp, dt,
      incub_params, generation_int_params,
      iso_delay_params, sec_infect_params,
      import_params, pd_params
    )
    
    repeat {
      # initialize outbreak at time zero
      sim_status <- initialize_sim_status(start_time, n_initial)
      state_df <- create_state_df(n_initial, sim_params, sim_status, initialize = TRUE)
      record_df <- create_record_df(state_df, sim_status, initialize = TRUE)
      # exit if more than one case
      # Continue outbreak until t max
      
      for (t in (1 + sim_status$t):tmax) { # t takes values 2,3,4,5,...  60 (tmax)
        out <- step_simulation(sim_status, state_df, record_df, sim_params)
        sim_status <- out$status # update sim_status
        state_df <- out$state # update state_df
        record_df <- out$record # update record_df
      }
      
      if (max(record_df$t_inf) >= pd_params$pd_change_t & nrow(record_df) > 3000) break ##### no stocastic extiction, each iteration repeats until an epidemic is sustained for 60 days, at least 3000 cases
      
    }
    
    
    return(record_df)
    
    
  } 
  
  
  out[[i]] <- sims_out
  
}


### CLEAN FOR ANALYSIS
set.seed(123)
out <- vector(mode = "list", length = 5)
for (i in 1:nrow(Rk_iteration)) {
  
  ## determine epidemic ID numbers
  epi_out_tbl <- list_rbind(sims_out[[i]]) |> 
    as_tibble() |> 
    (function(x) {x |> mutate(n_row = 1:nrow(x))}) ()
  
  epi_index <- epi_out_tbl |> 
    dplyr::select(case_id, n_row) |> 
    filter(case_id == 1)
  
  epi_case_count <- epi_index |> 
    dplyr::select(n_row_2 = n_row) |> 
    slice(2:nrow(epi_index)) |> 
    bind_rows(tibble(n_row_2 = nrow(epi_out_tbl))) |> 
    bind_cols(epi_index) |> 
    mutate(dif = n_row_2 - n_row)
  
  print("Assign epidemic ID")
  epi_id <- c()
  for (j in 1:nrow(epi_case_count)) {
    
    ind <- rep(j, times = epi_case_count$dif[j])
    
    epi_id <- c(epi_id, ind)
    
  }
  
  ## Add epidemic index to and split into list
  ## Calculate individual onset date from models infection time and incubation period. report date is uniform draw between 0 and 9 days
  epi_out_list <- epi_out_tbl |> 
    mutate(epidemic_id = c(epi_id, iter)) |> 
    group_split(epidemic_id) |> 
    map(function(x) x |> 
          mutate(t_symp = round(t_inf + d_incub)) |> 
          mutate(t_report = round(t_symp + runif(n = nrow(x), min = 0, max = 9))))
  
  
  
  ## Determine observed degree (infected) when control is applied on day t = 60
  ## First get true contact pairs
  epi_out_list_pairs <- map(epi_out_list, function(x) x |> 
                              slice(2:nrow(x)) |> 
                              transmute(
                                from = as.numeric(source),
                                to = as.numeric(case_id)))
  
  
  ## Calculate degree (observed infected)
  print("Calculating observed degree")
  epi_out_list_pairs_degree <- map2(epi_out_list, epi_out_list_pairs, function(x,y) make_epicontacts(x, y, directed = TRUE) |> 
                                      get_degree(type = "out"), .progress = TRUE)
  
  ## now append to final dataset
  epi_out_list_obs_n <- map2(epi_out_list, epi_out_list_pairs_degree, function(x, y) x |> 
                               left_join(tibble(obs_n = y, case_id = as.numeric(names(y))), by = "case_id"), .progress = TRUE)
  
  
  out[[i]] <- list_rbind(epi_out_list_obs_n) |>
    dplyr::select(epidemic_id, case_id, source, d_incub, t_inf, t_symp, t_report, n_sec_infects, obs_n) |> 
    mutate(controlled = n_sec_infects != obs_n,
           R = Rk_iteration$R[[i]],
           k = Rk_iteration$k[[i]])
  
  
  
  
  
}

saveRDS(out, file = "data/simulations/rds/sim-output-raw.rds")


