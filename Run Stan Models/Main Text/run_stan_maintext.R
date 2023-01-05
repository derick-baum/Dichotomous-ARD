#  "Uses and Limitations of Dichotomous Aggregate Relational Data"
#  R Script for running Stan Models
#
#  Author: Derick S.Baum <derick_baum@g.harvard.edu>

################################ Main Text ################################

#  Load packages ---------------------------------------------

library(rstan)
library(cmdstanr)
library(NSUM)

#  Subpopulation proportions ---------------------------------------------

data("McCarty")
pop_size <- 2.8e+08
mccarty_props <- McCarty$known / pop_size

#  Maltiel et al. random degree model (McCarty Data) ---------------------------------------------

## Count data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyCount_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyCount_dta.rds"
    )
  )

# set maximum ARD item responses as lower bounds for network size:
maximum_xiks <-
  apply(mccartyCount_dta, 1, max)

mccartyCount_MaltielRDM_StanData <-
  list(
    N = nrow(mccartyCount_dta),
    K = ncol(mccartyCount_dta),
    y = mccartyCount_dta,
    m = mccarty_props,
    L = maximum_xiks
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCount_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_count.stan"
  )

### Run Stan model ---------------------------------------------

mccartyCount_MaltielRDM_StanFit <-
  mccartyCount_MaltielRDM_StanMod$sample(
    data = mccartyCount_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyCount_MaltielRDM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyCount_MaltielRDM_StanFitObj <-
  read_stan_csv(mccartyCount_MaltielRDM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyCount_MaltielRDM_Extract <-
  rstan::extract(mccartyCount_MaltielRDM_StanFitObj)

## Dichotomous data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyDichotomous_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyDichotomous_dta.rds"
    )
  )

mccartyDichotomous_MaltielRDM_StanData <-
  list(
    N = nrow(mccartyDichotomous_dta),
    K = ncol(mccartyDichotomous_dta),
    y = mccartyDichotomous_dta,
    m = mccarty_props
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomous_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_dichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyDichotomous_MaltielRDM_StanFit <-
  mccartyDichotomous_MaltielRDM_StanMod$sample(
    data = mccartyDichotomous_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyDichotomous_MaltielRDM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyDichotomous_MaltielRDM_StanFitObj <-
  read_stan_csv(mccartyDichotomous_MaltielRDM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyDichotomous_MaltielRDM_Extract <-
  rstan::extract(mccartyDichotomous_MaltielRDM_StanFitObj)


#  Zheng et al. group prevalence model (McCarty Data) ---------------------------------------------

## Count data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyCount_ZhengGP_StanData <-
  list(
    N = nrow(mccartyCount_dta),
    K = ncol(mccartyCount_dta),
    y = mccartyCount_dta
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCount_ZhengGP_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. GP/ZhengGP_count.stan"
  )

### Run Stan model ---------------------------------------------

mccartyCount_ZhengGP_StanFit <-
  mccartyCount_ZhengGP_StanMod$sample(
    data = mccartyCount_ZhengGP_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyCount_ZhengGP_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyCount_ZhengGP_StanFitObj <-
  read_stan_csv(mccartyCount_ZhengGP_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyCount_ZhengGP_Extract <-
  rstan::extract(mccartyCount_ZhengGP_StanFitObj)

## Dichotomous data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyDichotomous_ZhengGP_StanData <-
  list(
    N = nrow(mccartyDichotomous_dta),
    K = ncol(mccartyDichotomous_dta),
    y = mccartyDichotomous_dta
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomous_ZhengGP_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. GP/ZhengGP_dichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyDichotomous_ZhengGP_StanFit <-
  mccartyDichotomous_ZhengGP_StanMod$sample(
    data = mccartyDichotomous_ZhengGP_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyDichotomous_ZhengGP_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyDichotomous_ZhengGP_StanFitObj <-
  read_stan_csv(mccartyDichotomous_ZhengGP_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyDichotomous_ZhengGP_Extract <-
  rstan::extract(mccartyDichotomous_ZhengGP_StanFitObj)

#  Maltiel et al. barrier effects model (McCarty Data) ---------------------------------------------

## Count data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyCount_MaltielBEM_StanData <-
  list(
    N = nrow(mccartyCount_dta),
    K = ncol(mccartyCount_dta),
    y = mccartyCount_dta,
    m = mccarty_props,
    L = maximum_xiks
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCount_MaltielBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. BEM/MaltielBEM_count.stan"
  )

### Run Stan model ---------------------------------------------

mccartyCount_MaltielBEM_StanFit <-
  mccartyCount_MaltielBEM_StanMod$sample(
    data = mccartyCount_MaltielBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

# If the command above takes too long to run,
# try running the model using the stan() function instead:

N <- nrow(mccartyCount_dta)
K <- ncol(mccartyCount_dta)
y <- mccartyCount_dta
m <- mccarty_props
L <- maximum_xiks

mccartyCount_MaltielBEM_StanFit_alt <-
  stan(
    file = "MaltielBEM_count.stan",
    data = c("N", "K", "y", "m", "L"),
    iter = 2000,
    chains = 4
  )

### Diagnostics ---------------------------------------------

mccartyCount_MaltielBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyCount_MaltielBEM_StanFitObj <-
  read_stan_csv(mccartyCount_MaltielBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyCount_MaltielBEM_Extract <-
  rstan::extract(mccartyCount_MaltielBEM_StanFitObj)

## Dichotomous data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyDichotomous_MaltielBEM_StanData <-
  list(
    N = nrow(mccartyDichotomous_dta),
    K = ncol(mccartyDichotomous_dta),
    y = mccartyDichotomous_dta,
    m = mccarty_props
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomous_MaltielBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. BEM/MaltielBEM_dichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyDichotomous_MaltielBEM_StanFit <-
  mccartyDichotomous_MaltielBEM_StanMod$sample(
    data = mccartyDichotomous_MaltielBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyDichotomous_MaltielBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyDichotomous_MaltielBEM_StanFitObj <-
  read_stan_csv(mccartyDichotomous_MaltielBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyDichotomous_MaltielBEM_Extract <-
  rstan::extract(mccartyDichotomous_MaltielBEM_StanFitObj)

## 0/1/2+ data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyTrichotomous_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyTrichotomous_dta.rds"
    )
  )

mccartyTrichotomous_MaltielBEM_StanData <-
  list(
    N = nrow(mccartyTrichotomous_dta),
    K = ncol(mccartyTrichotomous_dta),
    y = mccartyTrichotomous_dta,
    m = mccarty_props
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyTrichotomous_MaltielBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. BEM/MaltielBEM_trichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyTrichotomous_MaltielBEM_StanFit <-
  mccartyTrichotomous_MaltielBEM_StanMod$sample(
    data = mccartyTrichotomous_MaltielBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyTrichotomous_MaltielBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyTrichotomous_MaltielBEM_StanFitObj <-
  read_stan_csv(mccartyTrichotomous_MaltielBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyTrichotomous_MaltielBEM_Extract <-
  rstan::extract(mccartyTrichotomous_MaltielBEM_StanFitObj)

#  Zheng et al. barrier effects model (McCarty Data) ---------------------------------------------

## Count data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyCount_ZhengBEM_StanData <-
  list(
    N = nrow(mccartyCount_dta),
    K = ncol(mccartyCount_dta),
    y = mccartyCount_dta
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCount_ZhengBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. BEM/ZhengBEM_count.stan"
  )

### Run Stan model ---------------------------------------------

mccartyCount_ZhengBEM_StanFit <-
  mccartyCount_ZhengBEM_StanMod$sample(
    data = mccartyCount_ZhengBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyCount_ZhengBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyCount_ZhengBEM_StanFitObj <-
  read_stan_csv(mccartyCount_ZhengBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyCount_ZhengBEM_Extract <-
  rstan::extract(mccartyCount_ZhengBEM_StanFitObj)

## Dichotomous data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyDichotomous_ZhengBEM_StanData <-
  list(
    N = nrow(mccartyDichotomous_dta),
    K = ncol(mccartyDichotomous_dta),
    y = mccartyDichotomous_dta
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomous_ZhengBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. BEM/ZhengBEM_dichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyDichotomous_ZhengBEM_StanFit <-
  mccartyDichotomous_ZhengBEM_StanMod$sample(
    data = mccartyDichotomous_ZhengBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyDichotomous_ZhengBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyDichotomous_ZhengBEM_StanFitObj <-
  read_stan_csv(mccartyDichotomous_ZhengBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyDichotomous_ZhengBEM_Extract <-
  rstan::extract(mccartyDichotomous_ZhengBEM_StanFitObj)

## 0/1/2+ data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyTrichotomous_ZhengBEM_StanData <-
  list(
    N = nrow(mccartyTrichotomous_dta),
    K = ncol(mccartyTrichotomous_dta),
    y = mccartyTrichotomous_dta
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyTrichotomous_ZhengBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. BEM/ZhengBEM_trichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyTrichotomous_ZhengBEM_StanFit <-
  mccartyTrichotomous_ZhengBEM_StanMod$sample(
    data = mccartyTrichotomous_ZhengBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

### Diagnostics ---------------------------------------------

mccartyTrichotomous_ZhengBEM_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyTrichotomous_ZhengBEM_StanFitObj <-
  read_stan_csv(mccartyTrichotomous_ZhengBEM_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyTrichotomous_ZhengBEM_Extract <-
  rstan::extract(mccartyTrichotomous_ZhengBEM_StanFitObj)

#  Covariate-based model (McCarty Data) ---------------------------------------------

## Count data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyCountCov_dta <-
  readRDS("Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyCountCov_dta.rds")

# Load packages for data preparation:

library(dplyr)
library(fastDummies)

mccartyCountCov_dta <- mccartyCountCov_dta %>%
  mutate(across(
    isex,
    ~ case_when(. == 1 ~ "Male",
                . == 0 ~ "Female",
                TRUE ~ NA_character_)
  )) %>%
  # center age variable
  mutate(age = as.vector(scale(age, center = TRUE, scale = FALSE)))

dummies_sex <-
  fastDummies::dummy_cols(mccartyCountCov_dta[, 1])

dummies_sex_renamed <- dummies_sex %>%
  dplyr::select(!where(is.factor) & !where(is.character)) %>%
  rename(female = isex_Female,
         male = isex_Male)

y_count <- as.matrix(mccartyCountCov_dta[, 3:14])
male <- dummies_sex_renamed$male

# center sex dummy:
male_center <- as.vector(scale(male, center = TRUE, scale = FALSE))
age <- mccartyCountCov_dta$age

mccartyCount_CovMod_StanData <-
  list(
    N = nrow(y_count),
    K = ncol(y_count),
    y = y_count,
    age = age,
    male = male_center
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCount_CovMod_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Covariate Model/CovMod_count.stan"
  )

### Run Stan model ---------------------------------------------

mccartyCount_CovMod_StanFit <-
  mccartyCount_CovMod_StanMod$sample(
    data = mccartyCount_CovMod_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    init = 0
  )

### Diagnostics ---------------------------------------------

mccartyCount_CovMod_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyCount_CovMod_StanFitObj <-
  read_stan_csv(mccartyCount_CovMod_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyCount_CovMod_Extract <-
  rstan::extract(mccartyCount_CovMod_StanFitObj)

## Dichotomous data ---------------------------------------------

### Prepare Stan data ---------------------------------------------

mccartyDichotomousCov_dta <-
  readRDS(
    "Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyDichotomousCov_dta.rds"
  )

y_dichotomous <- as.matrix(mccartyDichotomousCov_dta[, 3:14])

mccartyDichotomous_CovMod_StanData <-
  list(
    N = nrow(y_dichotomous),
    K = ncol(y_dichotomous),
    y = y_dichotomous,
    age = age,
    male = male_center
  )

### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomous_CovMod_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Covariate Model/CovMod_dichotomous.stan"
  )

### Run Stan model ---------------------------------------------

mccartyDichotomous_CovMod_StanFit <-
  mccartyDichotomous_CovMod_StanMod$sample(
    data = mccartyDichotomous_CovMod_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    init = 0
  )

### Diagnostics ---------------------------------------------

mccartyDichotomous_CovMod_StanFit$cmdstan_diagnose()

### Create stanfit object ---------------------------------------------

mccartyDichotomous_CovMod_StanFitObj <-
  read_stan_csv(mccartyDichotomous_CovMod_StanFit$output_files())

### Extract samples ---------------------------------------------

mccartyDichotomous_CovMod_Extract <-
  rstan::extract(mccartyDichotomous_CovMod_StanFitObj)
