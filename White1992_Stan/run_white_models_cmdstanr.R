# Run from the same folder as these Stan files and the JSON.
# Requires CmdStan + R packages cmdstanr, posterior, and jsonlite.
# No microscopy-stage sample size is invented.

library(cmdstanr)
library(posterior)
library(jsonlite)

setwd("D:\\Data-Work\\Staging\\code\\Stan\\White1992_Stan_First_Fit_Package")

input <- jsonlite::fromJSON("white_stage_M8.json", simplifyVector = TRUE)
ref <- read.csv("White1992_12Patients_Mathematica_patient_fits.csv")
dir.create("stan_outputs", showWarnings = FALSE)

clamp <- function(x, low, high) pmin(high, pmax(low, x))

run_white <- function(model = c("C", "B"), mult = 8, seed = 20261015) {
  model <- match.arg(model)
  stanfile <- if (model == "C") {
    "white_stage_model3_gradual26.stan"
  } else {
    "white_stage_model2_S8.stan"
  }
  model_path <- cmdstan_model(stanfile)
  data <- input
  data$M <- as.numeric(mult)

  ref_i <- subset(ref, M == mult & Model == model)
  ref_i <- ref_i[order(ref_i$Patient), ]
  stopifnot(nrow(ref_i) == data$N,
            identical(as.integer(ref_i$Patient), seq_len(data$N)))

  initials <- function() {
    pars <- list(
      mu = clamp((ref_i$Mu_hours + rnorm(data$N, 0, 0.7)) %% 48, 0.01, 47.99),
      sigma = clamp(ref_i$Sigma_hours * exp(rnorm(data$N, 0, 0.05)), 1.01, 17.99),
      stage_error_sd = 0.04
    )
    if (model == "C") {
      ks <- ref_i$k_per_hour[1]
      pars$k_s <- clamp(ks * exp(rnorm(1, 0, 0.05)), 0.051, 1.99)
    }
    pars
  }

  fit <- model_path$sample(
    data = data, chains = 4, parallel_chains = 4,
    iter_warmup = 3000, iter_sampling = 3000,
    seed = seed, init = initials,
    adapt_delta = 0.97, max_treedepth = 12,
    output_dir = "stan_outputs"
  )
  filename <- sprintf("stan_outputs/%s_M%s_posterior_summary.csv", model, mult)
  variables <- c("mu", "sigma", "stage_error_sd", "q_circulating",
                 "total_density_N0", "sequestered_density")
  if (model == "C") variables <- c(variables, "k_s")
  summary <- fit$summary(variables = variables)
  write.csv(summary, filename, row.names = FALSE)
  cat("\n", model, " M=", mult, "\n", sep="")
  print(fit$diagnostic_summary())
  print(summary[summary$variable %in% c("k_s", "stage_error_sd"), ])
  invisible(fit)
}

# START WITH M=8. Run BOTH models before sensitivity scenarios.
fit_C_M8 <- run_white("C", 8)
fit_B_M8 <- run_white("B", 8)

# Later, after checking divergences, R-hat, ESS, posterior predictive fit:
# fits_C <- lapply(c(3,5,8,10,15), function(m) run_white("C",m))
# fits_B <- lapply(c(3,5,8,10,15), function(m) run_white("B",m))

library(bayesplot)

mcmc_trace(
  fit_B_M8$draws(
    variables = c(
      "mu[1]", "mu[4]",
      "sigma[1]", "sigma[4]"
    )
  )
)
