# FIXED PARAMS

##incubation period - approx Weibull dist mean = 4.82 and sd = 2.45 from  https://doi.org/10.1093/infdis/jiab424
w_shape = (sqrt(log(4 + (2.45^2 / 4.82^2))))^2
w_scale = 4.82 / gamma(1 + (1/w_shape))

## generation time. gamma distribution with mean = 5.71 and sd = 1.72 from  https://doi.org/10.1093/infdis/jiab424
g_shape <- gamma_mucv2shapescale(mu = 5.71, cv = 1.72 / 5.71)$shape
g_scale <- gamma_mucv2shapescale(mu = 5.71, cv = 1.72 / 5.71)$scale

infect_dur <- 21 # in days
p_symp <- 0 # percentage of cases that will be symptomatic

incub_params <- list(dist = "weibull", shape = w_shape, scale = w_scale) # distribution returns a value in days
generation_int_params <- list(dist = "gamma", shape = g_shape, rate = 1/g_scale) # distribution returns a value in days

# Contact tracing parameters - assumes uncontrolled epidemic growth. 
do_variable_trace <- FALSE # use a constant tracing value (see next line)
p_trace <- 0 # all cases have a 0% chance of being traced
p_trace_app <- 0 # percentage of population using a contact tracing app
p_trace_app_comp <- 0 # percentage of population complying with contact tracing app isolation instructions



# Set delays to isolation
iso_delay_params <- list(
  dist = "uniform", # the following are min/max ranges of a uniform distribution
  traced_min = 21, # traced secondary cases isolated 1-2 days after index case isolated
  traced_max = 21,
  untraced_min = 21, # untraced secondary cases isolated 3-4 days after *symptom onset*
  untraced_max = 21,
  untraced_pd_min = 21, # untraced secondary cases which are distancing isolated earlier
  untraced_pd_max = 21
)

# Physical (or social) distancing parameters. We don't apply control to the data generating process. as the sec cases are drawn from e.g. v = R0 * 0.8 where each individual distancing reduces their number of contacts as a proportion e.g. 80% of typical contacts R0 = 3; v = 2.4, but the dispersion is fixed thus the generation process is defining the offspring now. 
pd_params <- list(
  pd_pop_frac = 1, # 100% of population is effected vs unaffected by 
  pd_contact_rate1 = 1.0, # initial fraction of contacts for distancing group
  pd_contact_rate2 = 0.25, # contact rate change at t = 60 is Rt = 0.5 (2 * 0.25 = 0.5)
  pd_change_t = 60 # simulation day number where the contact rate changes
)


# The potential number of secondary infections is drawn from a negative binomial distribution with a mean value of R0 * contact_rate

# Define simulation control parameters (including imports)
start_time <- 0 # start at time 0 (recommended)
dt <- 1 # timestep, in days
tmax <- 1000 # total number of days in simulation
n_initial <- 1 # number of cases on day 0

## No importations during the period. Only care about local transmissions
import_params <- list(
  type = "constant",
  rate = 0, # no new import every day
  iso_lengths = 0, # number of days that imported cases might self-isolate upon arrival
  iso_p.group = 1 # probability that an imported case will do the above self-isolation
)
