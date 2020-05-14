library(epimod)

# model_generation(net_fname = "Net/Pertussis.PNPRO",
#                  functions_fname = "Cpp/transitions.cpp")
# 
sensitivity_analysis(n_config = 2^14,
                     parameters_fname = "./input/Functions_list.csv",
                     functions_fname = "./Rfunction/Functions.R",
                     solver_fname = "./Net/Pertussis.solver",
                     f_time = 365*21,
                     s_time = 365,
                     timeout = "1d",
                     parallel_processors=16,
                     reference_data ="./input/reference_data.csv",
                     distance_measure_fname="./Rfunction/msqd.R",
                     target_value_fname="./Rfunction/Select.R")
# 
# Optimization to get someting close to the mean
# model_calibration(parameters_fname = "./input/Functions_list.csv",
#                   functions_fname = "./Rfunction/Functions.R",
#                   solver_fname = "./Net/Pertussis.solver",
#                   f_time = 365*21,
#                   s_time = 365,
#                   reference_data = "./input/reference_data.csv",
#                   distance_measure_fname = "./Rfunction/msqd.R",
#                   # Vectors to control the optimization
#                   ini_v = c(0.0025,0.0031,0.0023,
#                             0.8615141,0.0507066,
#                             0.16041882,0.175,0.0181,0.02621356,0.04956218,
#                             0.99994,0.9987036,0.10108811,0.1520784,0.2098173),
#                   ub_v = c(0.0025,0.01,0.0025,
#                            1, 1,
#                            1, 1, 1, 1, 1,
#                            1, 1, 1, 1, 1),
#                   lb_v = c(0,0.0025,0,
#                            1e-7, 1e-7,
#                            1e-7, 1e-7, 1e-7, 1e-7, 1e-7,
#                            1e-7, 1e-7, 1e-7, 1e-7, 1e-7),
#                   ini_vector_mod = TRUE
# )

# results_deterministic_model_calibration
optim <- c(0.002474758,0.002537443,0.002458887,
           0.931635,7.669116e-06,
           0.9963428,2.133152e-07,1e-07,1e-07,1e-07,
           0.9999986,0.005559236,1e-07,1e-07,1e-07)

# Optimization to get someting close to the mean
model_calibration(parameters_fname = "./input/Functions_list.csv",
                  functions_fname = "./Rfunction/Functions.R",
                  solver_fname = "./Net/Pertussis.solver",
                  f_time = 365*21,
                  s_time = 365,
                  n_run = 2^8,
                  reference_data = "./input/reference_data.csv",
                  distance_measure_fname = "./Rfunction/aic.R",
                  # Vectors to control the optimization
                  ini_v = optim[1:3],
                  ub_v = c(rep(0.0026,2),0.0025),
                  lb_v = c(rep(0.0025,2),0.0024),
                  ini_vector_mod = TRUE,
                  solver_type = "TAUG",
                  parallel_processors = 16
)

rnk <- read.csv(file="~/results_stochastic_model_calibration/Pertussis-calibration_optim-trace.csv", sep = "")
rnk <- rnk[order(rnk$distance),]
optim <- c(rnk[1,3:5],
           0.931635,7.669116e-06,
           0.9963428,2.133152e-07,1e-07,1e-07,1e-07,
           0.9999986,0.005559236,1e-07,1e-07,1e-07)

model_analysis(solver_fname = "./Net/Pertussis.solver",
               f_time = 365*43,
               s_time = 365,
               n_config = 1,
               n_run = 2^10,
               parallel_processors = 16,
               solver_type = "TAUG",
               parameters_fname = "./input/Functions_list.csv",
               functions_fname = "./Rfunction/Functions.R",
               ini_v = optim ,
               ini_vector_mod = TRUE)

