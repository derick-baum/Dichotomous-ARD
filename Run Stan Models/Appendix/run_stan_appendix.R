#  "Uses and Limitations of Dichotomous Aggregate Relational Data"
#  R Script for running Stan Models
#
#  Author: Derick S.Baum <derick_baum@g.harvard.edu>

################################ Appendix ################################

#  Load packages ---------------------------------------------

library(rstan)
library(cmdstanr)
library(NSUM)

#  Appendix A: Results with the Spanish Data ---------------------------------------------

##  Subpopulation proportions ---------------------------------------------

spanish_props <-
  readRDS("Dichotomous-ARD-Replication-Materials/Data/Appendix/spanishProps.rds")
spanish_props <- head(spanish_props$prop,-1)

##  Maltiel et al. random degree model (Spanish Data) ---------------------------------------------

### Count data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishCount_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Appendix/spanishCount_dta.rds"
    )
  )

# set maximum ARD item responses as lower bounds for network size:
maximum_xiks <-
  apply(spanishCount_dta, 1, max)

spanishCount_MaltielRDM_StanData <-
  list(
    N = nrow(spanishCount_dta),
    K = ncol(spanishCount_dta),
    y = spanishCount_dta,
    m = spanish_props,
    L = maximum_xiks
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishCount_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_count.stan"
  )

#### Run Stan model ---------------------------------------------

spanishCount_MaltielRDM_StanFit <-
  spanishCount_MaltielRDM_StanMod$sample(
    data = spanishCount_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishCount_MaltielRDM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishCount_MaltielRDM_StanFitObj <-
  read_stan_csv(spanishCount_MaltielRDM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishCount_MaltielRDM_Extract <-
  rstan::extract(spanishCount_MaltielRDM_StanFitObj)

### Dichotomous data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishDichotomous_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Appendix/spanishDichotomous_dta.rds"
    )
  )

spanishDichotomous_MaltielRDM_StanData <-
  list(
    N = nrow(spanishDichotomous_dta),
    K = ncol(spanishDichotomous_dta),
    y = spanishDichotomous_dta,
    m = spanish_props
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishDichotomous_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_dichotomous.stan"
  )

#### Run Stan model ---------------------------------------------

spanishDichotomous_MaltielRDM_StanFit <-
  spanishDichotomous_MaltielRDM_StanMod$sample(
    data = spanishDichotomous_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishDichotomous_MaltielRDM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishDichotomous_MaltielRDM_StanFitObj <-
  read_stan_csv(spanishDichotomous_MaltielRDM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishDichotomous_MaltielRDM_Extract <-
  rstan::extract(spanishDichotomous_MaltielRDM_StanFitObj)

##  Zheng et al. group prevalence model (McCarty Data) ---------------------------------------------

### Count data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishCount_ZhengGP_StanData <-
  list(
    N = nrow(spanishCount_dta),
    K = ncol(spanishCount_dta),
    y = spanishCount_dta
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishCount_ZhengGP_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. GP/ZhengGP_count.stan"
  )

#### Run Stan model ---------------------------------------------

spanishCount_ZhengGP_StanFit <-
  spanishCount_ZhengGP_StanMod$sample(
    data = spanishCount_ZhengGP_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishCount_ZhengGP_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishCount_ZhengGP_StanFitObj <-
  read_stan_csv(spanishCount_ZhengGP_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishCount_ZhengGP_Extract <-
  rstan::extract(spanishCount_ZhengGP_StanFitObj)

### Dichotomous data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishDichotomous_ZhengGP_StanData <-
  list(
    N = nrow(spanishDichotomous_dta),
    K = ncol(spanishDichotomous_dta),
    y = spanishDichotomous_dta
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishDichotomous_ZhengGP_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. GP/ZhengGP_dichotomous.stan"
  )

#### Run Stan model ---------------------------------------------

spanishDichotomous_ZhengGP_StanFit <-
  spanishDichotomous_ZhengGP_StanMod$sample(
    data = spanishDichotomous_ZhengGP_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishDichotomous_ZhengGP_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishDichotomous_ZhengGP_StanFitObj <-
  read_stan_csv(spanishDichotomous_ZhengGP_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishDichotomous_ZhengGP_Extract <-
  rstan::extract(spanishDichotomous_ZhengGP_StanFitObj)

##  Maltiel et al. barrier effects model (McCarty Data) ---------------------------------------------

### Count data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishCount_MaltielBEM_StanData <-
  list(
    N = nrow(spanishCount_dta),
    K = ncol(spanishCount_dta),
    y = spanishCount_dta,
    m = spanish_props,
    L = maximum_xiks
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishCount_MaltielBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. BEM/MaltielBEM_count.stan"
  )

#### Run Stan model ---------------------------------------------

spanishCount_MaltielBEM_StanFit <-
  spanishCount_MaltielBEM_StanMod$sample(
    data = spanishCount_MaltielBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishCount_MaltielBEM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishCount_MaltielBEM_StanFitObj <-
  read_stan_csv(spanishCount_MaltielBEM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishCount_MaltielBEM_Extract <-
  rstan::extract(spanishCount_MaltielBEM_StanFitObj)

### Dichotomous data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishDichotomous_MaltielBEM_StanData <-
  list(
    N = nrow(spanishDichotomous_dta),
    K = ncol(spanishDichotomous_dta),
    y = spanishDichotomous_dta,
    m = spanish_props
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishDichotomous_MaltielBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. BEM/MaltielBEM_dichotomous.stan"
  )

#### Run Stan model ---------------------------------------------

spanishDichotomous_MaltielBEM_StanFit <-
  spanishDichotomous_MaltielBEM_StanMod$sample(
    data = spanishDichotomous_MaltielBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishDichotomous_MaltielBEM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishDichotomous_MaltielBEM_StanFitObj <-
  read_stan_csv(spanishDichotomous_MaltielBEM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishDichotomous_MaltielBEM_Extract <-
  rstan::extract(spanishDichotomous_MaltielBEM_StanFitObj)

##  Zheng et al. barrier effects model (spanish Data) ---------------------------------------------

### Count data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishCount_ZhengBEM_StanData <-
  list(
    N = nrow(spanishCount_dta),
    K = ncol(spanishCount_dta),
    y = spanishCount_dta
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishCount_ZhengBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. BEM/ZhengBEM_count.stan"
  )

#### Run Stan model ---------------------------------------------

spanishCount_ZhengBEM_StanFit <-
  spanishCount_ZhengBEM_StanMod$sample(
    data = spanishCount_ZhengBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishCount_ZhengBEM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishCount_ZhengBEM_StanFitObj <-
  read_stan_csv(spanishCount_ZhengBEM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishCount_ZhengBEM_Extract <-
  rstan::extract(spanishCount_ZhengBEM_StanFitObj)

### Dichotomous data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

spanishDichotomous_ZhengBEM_StanData <-
  list(
    N = nrow(spanishDichotomous_dta),
    K = ncol(spanishDichotomous_dta),
    y = spanishDichotomous_dta
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
spanishDichotomous_ZhengBEM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Zheng et al. BEM/ZhengBEM_dichotomous.stan"
  )

#### Run Stan model ---------------------------------------------

spanishDichotomous_ZhengBEM_StanFit <-
  spanishDichotomous_ZhengBEM_StanMod$sample(
    data = spanishDichotomous_ZhengBEM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

spanishDichotomous_ZhengBEM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

spanishDichotomous_ZhengBEM_StanFitObj <-
  read_stan_csv(spanishDichotomous_ZhengBEM_StanFit$output_files())

#### Extract samples ---------------------------------------------

spanishDichotomous_ZhengBEM_Extract <-
  rstan::extract(spanishDichotomous_ZhengBEM_StanFitObj)

#  Appendix B: Number of Subpopulations ---------------------------------------------

##  Subpopulation proportions ---------------------------------------------

data("McCarty")
pop_size <- 2.8e+08
mccarty_props <- McCarty$known / pop_size
mccarty_propsReduced <- mccarty_props[1:12]

##  Maltiel et al. random degree model (Reduced McCarty Data) ---------------------------------------------

### Count data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

mccartyCountReduced_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Appendix/mccartyCountReduced_dta.rds"
    )
  )

# set maximum ARD item responses as lower bounds for network size:
maximum_xiks <-
  apply(mccartyCountReduced_dta, 1, max)

mccartyCountReduced_MaltielRDM_StanData <-
  list(
    N = nrow(mccartyCountReduced_dta),
    K = ncol(mccartyCountReduced_dta),
    y = mccartyCountReduced_dta,
    m = mccarty_propsReduced,
    L = maximum_xiks
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyCountReduced_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_count.stan"
  )

#### Run Stan model ---------------------------------------------

mccartyCountReduced_MaltielRDM_StanFit <-
  mccartyCountReduced_MaltielRDM_StanMod$sample(
    data = mccartyCountReduced_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

mccartyCountReduced_MaltielRDM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

mccartyCountReduced_MaltielRDM_StanFitObj <-
  read_stan_csv(mccartyCountReduced_MaltielRDM_StanFit$output_files())

#### Extract samples ---------------------------------------------

mccartyCountReduced_MaltielRDM_Extract <-
  rstan::extract(mccartyCountReduced_MaltielRDM_StanFitObj)

### Dichotomous data ---------------------------------------------

#### Prepare Stan data ---------------------------------------------

mccartyDichotomousReduced_dta <-
  as.matrix(
    readRDS(
      "Dichotomous-ARD-Replication-Materials/Data/Appendix/mccartyDichotomousReduced_dta.rds"
    )
  )

mccartyDichotomousReduced_MaltielRDM_StanData <-
  list(
    N = nrow(mccartyDichotomousReduced_dta),
    K = ncol(mccartyDichotomousReduced_dta),
    y = mccartyDichotomousReduced_dta,
    m = mccarty_propsReduced
  )

#### Stan setup ---------------------------------------------

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
mccartyDichotomousReduced_MaltielRDM_StanMod <-
  cmdstan_model(
    "Dichotomous-ARD-Replication-Materials/Stan Files/Maltiel et al. RDM/MaltielRDM_dichotomous.stan"
  )

#### Run Stan model ---------------------------------------------

mccartyDichotomousReduced_MaltielRDM_StanFit <-
  mccartyDichotomousReduced_MaltielRDM_StanMod$sample(
    data = mccartyDichotomousReduced_MaltielRDM_StanData,
    seed = 123,
    chains = 4,
    refresh = 200,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000
  )

#### Diagnostics ---------------------------------------------

mccartyDichotomousReduced_MaltielRDM_StanFit$cmdstan_diagnose()

#### Create stanfit object ---------------------------------------------

mccartyDichotomousReduced_MaltielRDM_StanFitObj <-
  read_stan_csv(mccartyDichotomousReduced_MaltielRDM_StanFit$output_files())

#### Extract samples ---------------------------------------------

mccartyDichotomousReduced_MaltielRDM_Extract <-
  rstan::extract(mccartyDichotomousReduced_MaltielRDM_StanFitObj)
