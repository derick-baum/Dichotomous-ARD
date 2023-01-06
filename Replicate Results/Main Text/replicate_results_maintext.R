#  "Uses and Limitations of Dichotomous Aggregate Relational Data"
#  R Script for replicating results in the article
#
#  Author: Derick S.Baum <derick_baum@g.harvard.edu>

################################ Main Text ################################

#  Load packages ---------------------------------------------

library(tidyverse)
library(NSUM)
library(ggrepel)
library(naniar)
library(scales)

#  Figure 1 ---------------------------------------------

## Extract degree samples  ---------------------------------------------

# Count:

mccartyCount_MaltielRDM_ExtractDegrees <-
  mccartyCount_MaltielRDM_Extract$d

# Dichotomous:

mccartyDichotomous_MaltielRDM_ExtractDegrees <-
  mccartyDichotomous_MaltielRDM_Extract$d

## Compute posterior means  ---------------------------------------------

# Count:

mccartyCount_MaltielRDM_DegreeEsts <-
  colMeans(mccartyCount_MaltielRDM_ExtractDegrees)

# Dichotomous:

mccartyDichotomous_MaltielRDM_DegreeEsts <-
  colMeans(mccartyDichotomous_MaltielRDM_ExtractDegrees)

## Create plot  ---------------------------------------------

# Joint count and dichotomous estimates into a sigle data frame:

mccarty_MaltielRDM_DegreeEsts_df <-
  data.frame(networkdegree_Count = mccartyCount_MaltielRDM_DegreeEsts,
             networkdegree_Dichotomous = mccartyDichotomous_MaltielRDM_DegreeEsts)

# Lengthen the data:

mccarty_MaltielRDM_DegreeEsts_df_long <-
  mccarty_MaltielRDM_DegreeEsts_df %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "data_type"),
    names_sep = "_"
  )

mccarty_MaltielRDM_DegreeEsts_df_long$data_type <-
  factor(mccarty_MaltielRDM_DegreeEsts_df_long$data_type,
         levels = c("Count", "Dichotomous"))

# Create kernel density plot:

mccarty_MaltielRDM_DegreeEsts_KDP <-
  ggplot(
    mccarty_MaltielRDM_DegreeEsts_df_long,
    aes(networkdegree, color = data_type, linetype = data_type)
  ) +
  geom_density() +
  coord_cartesian(xlim = c(0, 1000)) +
  labs(color = "",
       linetype = "",
       x = "Posterior Mean Network Size") +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  scale_linetype_manual(values = c("Count" = "solid", "Dichotomous" = "dashed")) +
  theme_classic() +
  theme(legend.text = element_text(size = 10))

#  Table 1 ---------------------------------------------

summary(mccartyCount_MaltielRDM_DegreeEsts)
summary(mccartyDichotomous_MaltielRDM_DegreeEsts)

#  Figure 2 ---------------------------------------------

mccarty_MaltielRDM_DegreeEsts_Scatterdf <-
  data.frame(Count = mccartyCount_MaltielRDM_DegreeEsts,
             Dichotomous = mccartyDichotomous_MaltielRDM_DegreeEsts)

mccarty_MaltielRDM_DegreeEsts_Scatter <-
  ggplot(mccarty_MaltielRDM_DegreeEsts_Scatterdf,
         aes(x = Count, y = Dichotomous)) +
  geom_jitter() +
  theme_classic()

#  Table 2 ---------------------------------------------

# Count:

mccartyCount_MaltielRDM_QuantileSamples <-
  apply(
    mccartyCount_MaltielRDM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyCount_MaltielRDM_QuantileEsts <-
  apply(mccartyCount_MaltielRDM_QuantileSamples, 1, median)
mccartyCount_MaltielRDM_QuantileIQRs <-
  apply(mccartyCount_MaltielRDM_QuantileSamples, 1, IQR)

# Dichotomous:

mccartyDichotomous_MaltielRDM_QuantileSamples <-
  apply(
    mccartyDichotomous_MaltielRDM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyDichotomous_MaltielRDM_QuantileEsts <-
  apply(mccartyDichotomous_MaltielRDM_QuantileSamples, 1, median)
mccartyDichotomous_MaltielRDM_QuantileIQRs <-
  apply(mccartyDichotomous_MaltielRDM_QuantileSamples, 1, IQR)

mccartyCount_MaltielRDM_QuantileEsts
mccartyDichotomous_MaltielRDM_QuantileEsts
mccartyCount_MaltielRDM_QuantileIQRs
mccartyDichotomous_MaltielRDM_QuantileIQRs

#  Figure 3 ---------------------------------------------

## McCarty et al. data ---------------------------------------------

### Renormalize betas for all subpopulations ---------------------------------------------

# Setup:

data("McCarty")
pop_size <- 2.8e+08
mccarty_props <- McCarty$known / pop_size

mccarty_sumprops <- sum(mccarty_props)
mccarty_indsubpops <- 1:29

#### Count ---------------------------------------------

# Extract betas:

mccartyCount_ZhengGP_ExtractBetas <-
  mccartyCount_ZhengGP_Extract$beta

# Compute renormalization constant:

mccartyCount_ZhengGP_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyCount_ZhengGP_ExtractBetas)) {
  mccartyCount_ZhengGP_Constant[i] = log(sum(exp(mccartyCount_ZhengGP_ExtractBetas[i, mccarty_indsubpops])) / mccarty_sumprops)
}

# Renormalize betas:

mccartyCount_ZhengGP_NormalizedBetas <-
  matrix(NA, nrow = 4000, ncol = 29)

for (i in 1:nrow(mccartyCount_ZhengGP_ExtractBetas)) {
  mccartyCount_ZhengGP_NormalizedBetas[i,] = mccartyCount_ZhengGP_ExtractBetas[i,] - mccartyCount_ZhengGP_Constant[i]
}

#### Dichotomous ---------------------------------------------

# Extract betas:

mccartyDichotomous_ZhengGP_ExtractBetas <-
  mccartyDichotomous_ZhengGP_Extract$beta

# Compute renormalization constant:

mccartyDichotomous_ZhengGP_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyDichotomous_ZhengGP_ExtractBetas)) {
  mccartyDichotomous_ZhengGP_Constant[i] = log(sum(exp(
    mccartyDichotomous_ZhengGP_ExtractBetas[i, mccarty_indsubpops]
  )) / mccarty_sumprops)
}

# Renormalize betas:

mccartyDichotomous_ZhengGP_NormalizedBetas <-
  matrix(NA, nrow = 4000, ncol = 29)

for (i in 1:nrow(mccartyDichotomous_ZhengGP_ExtractBetas)) {
  mccartyDichotomous_ZhengGP_NormalizedBetas[i,] = mccartyDichotomous_ZhengGP_ExtractBetas[i,] - mccartyDichotomous_ZhengGP_Constant[i]
}

### Compute posterior means  ---------------------------------------------

# Count:

mccartyCount_ZhengGP_NormalizedExpBetas <-
  exp(mccartyCount_ZhengGP_NormalizedBetas)
mccartyCount_ZhengGP_NormalizedExpBetaEsts <-
  colMeans(mccartyCount_ZhengGP_NormalizedExpBetas)

# Dichotomous:

mccartyDichotomous_ZhengGP_NormalizedExpBetas <-
  exp(mccartyDichotomous_ZhengGP_NormalizedBetas)
mccartyDichotomous_ZhengGP_NormalizedExpBetaEsts <-
  colMeans(mccartyDichotomous_ZhengGP_NormalizedExpBetas)

## Lubbers et al. data ---------------------------------------------

### Renormalize betas for all subpopulations ---------------------------------------------

# Setup:

spanish_subpops <-
  readRDS("Dichotomous-ARD-Replication-Materials/Data/Appendix/spanishProps.rds")
spanish_props <- head(spanish_subpops$prop, -1)

spanish_sumprops <- sum(spanish_props)
spanish_indsubpops <- 1:14

#### Count ---------------------------------------------

# Extract betas:

spanishCount_ZhengGP_ExtractBetas <-
  spanishCount_ZhengGP_Extract$beta

# Compute renormalization constant:

spanishCount_ZhengGP_Constant <- c() # vector of constants

for (i in 1:nrow(spanishCount_ZhengGP_ExtractBetas)) {
  spanishCount_ZhengGP_Constant[i] = log(sum(exp(spanishCount_ZhengGP_ExtractBetas[i, spanish_indsubpops])) / spanish_sumprops)
}

# Renormalize betas:

spanishCount_ZhengGP_NormalizedBetas <-
  matrix(NA, nrow = 4000, ncol = 14)

for (i in 1:nrow(spanishCount_ZhengGP_ExtractBetas)) {
  spanishCount_ZhengGP_NormalizedBetas[i,] = spanishCount_ZhengGP_ExtractBetas[i,] - spanishCount_ZhengGP_Constant[i]
}

#### Dichotomous ---------------------------------------------

# Extract betas:

spanishDichotomous_ZhengGP_ExtractBetas <-
  spanishDichotomous_ZhengGP_Extract$beta

# Compute renormalization constant:

spanishDichotomous_ZhengGP_Constant <- c() # vector of constants

for (i in 1:nrow(spanishDichotomous_ZhengGP_ExtractBetas)) {
  spanishDichotomous_ZhengGP_Constant[i] = log(sum(exp(
    spanishDichotomous_ZhengGP_ExtractBetas[i, spanish_indsubpops]
  )) / spanish_sumprops)
}

# Renormalize betas:

spanishDichotomous_ZhengGP_NormalizedBetas <-
  matrix(NA, nrow = 4000, ncol = 14)

for (i in 1:nrow(spanishDichotomous_ZhengGP_ExtractBetas)) {
  spanishDichotomous_ZhengGP_NormalizedBetas[i,] = spanishDichotomous_ZhengGP_ExtractBetas[i,] - spanishDichotomous_ZhengGP_Constant[i]
}

### Compute posterior means  ---------------------------------------------

# Count:

spanishCount_ZhengGP_NormalizedExpBetas <-
  exp(spanishCount_ZhengGP_NormalizedBetas)
spanishCount_ZhengGP_NormalizedExpBetaEsts <-
  colMeans(spanishCount_ZhengGP_NormalizedExpBetas)

# Dichotomous:

spanishDichotomous_ZhengGP_NormalizedExpBetas <-
  exp(spanishDichotomous_ZhengGP_NormalizedBetas)
spanishDichotomous_ZhengGP_NormalizedExpBetaEsts <-
  colMeans(spanishDichotomous_ZhengGP_NormalizedExpBetas)

## Prepare data frame ---------------------------------------------

# McCarty et al. data:

mccarty_columnnames <- c("group", "Count",
                         "Dichotomous")

mccarty_subpopnames <- names(McCarty$known)
mccarty_subpopnames[26] <- "Male in prison"
mccarty_subpopnames[20] <- "Jaycees"

mccarty_ZhengGP_ExpBetaEsts_df <-
  data.frame(
    mccarty_subpopnames,
    mccartyCount_ZhengGP_NormalizedExpBetaEsts,
    mccartyDichotomous_ZhengGP_NormalizedExpBetaEsts
  )

colnames(mccarty_ZhengGP_ExpBetaEsts_df) <- mccarty_columnnames
mccarty_ZhengGP_ExpBetaEsts_df <-
  cbind(mccarty_ZhengGP_ExpBetaEsts_df, proportions = mccarty_props)

mccarty_ZhengGP_ExpBetaEsts_df_long <-
  pivot_longer(
    mccarty_ZhengGP_ExpBetaEsts_df,
    cols = 2:3,
    names_to = "data_type",
    values_to = "prevalences"
  )

mccarty_ZhengGP_ExpBetaEsts_df_long$Data <- "McCarty et al. Data"

mccarty_subsetindices <-
  c(1:2,
    5:6,
    11:12,
    17:18,
    21:22,
    25:26,
    27:28,
    31:32,
    41:42,
    45:46,
    47:48,
    51:52)

mccartySubset_ZhengGP_ExpBetaEsts_df_long <-
  mccarty_ZhengGP_ExpBetaEsts_df_long[mccarty_subsetindices, ]

# Lubbers et al. data:

spanish_columnnames <- c("group", "Count",
                         "Dichotomous")

spanish_ZhengGP_ExpBetaEsts_df <-
  data.frame(
    spanish_subpops$name[1:14],
    spanishCount_ZhengGP_NormalizedExpBetaEsts,
    spanishDichotomous_ZhengGP_NormalizedExpBetaEsts
  )

colnames(spanish_ZhengGP_ExpBetaEsts_df) <- spanish_columnnames
spanish_ZhengGP_ExpBetaEsts_df <-
  cbind(spanish_ZhengGP_ExpBetaEsts_df, proportions = spanish_props)

spanish_ZhengGP_ExpBetaEsts_df_long <-
  pivot_longer(
    spanish_ZhengGP_ExpBetaEsts_df,
    cols = 2:3,
    names_to = "data_type",
    values_to = "prevalences"
  )

spanish_ZhengGP_ExpBetaEsts_df_long$Data <- "Lubbers et al. Data"

spanish_subsetindices <- c(1:2, 5:6, 7:8, 11:12, 17:22, 25:26)

spanishSubset_ZhengGP_ExpBetaEsts_df_long <-
  spanish_ZhengGP_ExpBetaEsts_df_long[spanish_subsetindices, ]

# Join both data frames:

Joint_ZhengGP_ExpBetaEsts_df_long <-
  rbind(mccarty_ZhengGP_ExpBetaEsts_df_long,
        spanish_ZhengGP_ExpBetaEsts_df_long)

Joint_ZhengGP_ExpBetaEsts_df_long$Data <-
  factor(
    Joint_ZhengGP_ExpBetaEsts_df_long$Data,
    levels = c("McCarty et al. Data", "Lubbers et al. Data")
  )

SubsetJoint_ZhengGP_ExpBetaEsts_df_long <-
  rbind(
    mccartySubset_ZhengGP_ExpBetaEsts_df_long,
    spanishSubset_ZhengGP_ExpBetaEsts_df_long
  )

SubsetJoint_ZhengGP_ExpBetaEsts_df_long$Data <-
  factor(
    SubsetJoint_ZhengGP_ExpBetaEsts_df_long$Data,
    levels = c("McCarty et al. Data", "Lubbers et al. Data")
  )

## Create plot ---------------------------------------------

uniqueInitials <- c("C", "D")
initialShapes <- unname(sapply(uniqueInitials, utf8ToInt))

ZhengGP_ComparePrevalenceProportion <-
  ggplot(
    Joint_ZhengGP_ExpBetaEsts_df_long,
    aes(
      x = proportions,
      y = prevalences,
      color = data_type,
      shape = data_type
    )
  ) +
  geom_jitter(size = 3) +
  xlab(expression(Fraction ~ "in" ~ Population ~ "(" * italic(m[k]) * ")")) +
  ylab(expression(
    Posterior ~ Mean ~ Group ~ Prevalence ~ Estimate ~ "(" * italic(b[k]) * ")"
  )) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  expand_limits(x = 0, y = 0) +
  geom_label_repel(
    data = SubsetJoint_ZhengGP_ExpBetaEsts_df_long,
    aes(label = group),
    min.segment.length = 0,
    seed = 1234,
    box.padding = 0.3,
    size = 2.5,
    show.legend = FALSE
  ) +
  labs(color = "Group Prevalence",
       shape = "Group Prevalence") +
  scale_shape_manual(values = initialShapes) +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  facet_wrap(~ Data,
             scales = "free",
             nrow = 2) +
  theme_classic() +
  theme(legend.text = element_text(size = 10))

#  Table 3 ---------------------------------------------

## McCarty et al. data ---------------------------------------------

### Renormalize alphas for all subpopulations ---------------------------------------------

#### Count ---------------------------------------------

# Extract alphas:

mccartyCount_ZhengGP_ExtractAlphas <-
  mccartyCount_ZhengGP_Extract$alpha

# Renormalize alphas:

mccartyCount_ZhengGP_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyCount_ZhengGP_ExtractAlphas)) {
  mccartyCount_ZhengGP_NormalizedAlphas[i,] = mccartyCount_ZhengGP_ExtractAlphas[i,] + mccartyCount_ZhengGP_Constant[i]
}

#### Dichotomous ---------------------------------------------

# Extract alphas:

mccartyDichotomous_ZhengGP_ExtractAlphas <-
  mccartyDichotomous_ZhengGP_Extract$alpha

# Renormalize alphas:

mccartyDichotomous_ZhengGP_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyDichotomous_ZhengGP_ExtractAlphas)) {
  mccartyDichotomous_ZhengGP_NormalizedAlphas[i,] = mccartyDichotomous_ZhengGP_ExtractAlphas[i,] + mccartyDichotomous_ZhengGP_Constant[i]
}

### Compute posterior means  ---------------------------------------------

# Count:

mccartyCount_ZhengGP_NormalizedExpAlphas <-
  exp(mccartyCount_ZhengGP_NormalizedAlphas)
mccartyCount_ZhengGP_NormalizedDegreeEsts <-
  colMeans(mccartyCount_ZhengGP_NormalizedExpAlphas)

# Dichotomous:

mccartyDichotomous_ZhengGP_NormalizedExpAlphas <-
  exp(mccartyDichotomous_ZhengGP_NormalizedAlphas)
mccartyDichotomous_ZhengGP_NormalizedDegreeEsts <-
  colMeans(mccartyDichotomous_ZhengGP_NormalizedExpAlphas)

## Lubbers et al. data ---------------------------------------------

### Renormalize alphas for all subpopulations ---------------------------------------------

#### Count ---------------------------------------------

# Extract alphas:

spanishCount_ZhengGP_ExtractAlphas <-
  spanishCount_ZhengGP_Extract$alpha

# Renormalize alphas:

spanishCount_ZhengGP_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 2371)

for (i in 1:nrow(spanishCount_ZhengGP_ExtractAlphas)) {
  spanishCount_ZhengGP_NormalizedAlphas[i,] = spanishCount_ZhengGP_ExtractAlphas[i,] + spanishCount_ZhengGP_Constant[i]
}

#### Dichotomous ---------------------------------------------

# Extract alphas:

spanishDichotomous_ZhengGP_ExtractAlphas <-
  spanishDichotomous_ZhengGP_Extract$alpha

# Renormalize alphas:

spanishDichotomous_ZhengGP_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 2371)

for (i in 1:nrow(spanishDichotomous_ZhengGP_ExtractAlphas)) {
  spanishDichotomous_ZhengGP_NormalizedAlphas[i,] = spanishDichotomous_ZhengGP_ExtractAlphas[i,] + spanishDichotomous_ZhengGP_Constant[i]
}

### Compute posterior means  ---------------------------------------------

# Count:

spanishCount_ZhengGP_NormalizedExpAlphas <-
  exp(spanishCount_ZhengGP_NormalizedAlphas)
spanishCount_ZhengGP_NormalizedDegreeEsts <-
  colMeans(spanishCount_ZhengGP_NormalizedExpAlphas)

# Dichotomous:

spanishDichotomous_ZhengGP_NormalizedExpAlphas <-
  exp(spanishDichotomous_ZhengGP_NormalizedAlphas)
spanishDichotomous_ZhengGP_NormalizedDegreeEsts <-
  colMeans(spanishDichotomous_ZhengGP_NormalizedExpAlphas)

## Generate table results ---------------------------------------------

summary(mccartyCount_ZhengGP_NormalizedDegreeEsts)
summary(mccartyDichotomous_ZhengGP_NormalizedDegreeEsts)
summary(spanishCount_ZhengGP_NormalizedDegreeEsts)
summary(spanishDichotomous_ZhengGP_NormalizedDegreeEsts)

#  Table 4 ---------------------------------------------

## McCarty et al. data ---------------------------------------------

# Count:

mccartyCount_ZhengGP_QuantileSamples <-
  apply(
    mccartyCount_ZhengGP_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyCount_ZhengGP_QuantileEsts <-
  apply(mccartyCount_ZhengGP_QuantileSamples, 1, median)
mccartyCount_ZhengGP_QuantileIQRs <-
  apply(mccartyCount_ZhengGP_QuantileSamples, 1, IQR)

# Dichotomous:

mccartyDichotomous_ZhengGP_QuantileSamples <-
  apply(
    mccartyDichotomous_ZhengGP_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyDichotomous_ZhengGP_QuantileEsts <-
  apply(mccartyDichotomous_ZhengGP_QuantileSamples, 1, median)
mccartyDichotomous_ZhengGP_QuantileIQRs <-
  apply(mccartyDichotomous_ZhengGP_QuantileSamples, 1, IQR)

## Lubbers et al. data ---------------------------------------------

# Count:

spanishCount_ZhengGP_QuantileSamples <-
  apply(
    spanishCount_ZhengGP_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishCount_ZhengGP_QuantileEsts <-
  apply(spanishCount_ZhengGP_QuantileSamples, 1, median)
spanishCount_ZhengGP_QuantileIQRs <-
  apply(spanishCount_ZhengGP_QuantileSamples, 1, IQR)

# Dichotomous:

spanishDichotomous_ZhengGP_QuantileSamples <-
  apply(
    spanishDichotomous_ZhengGP_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishDichotomous_ZhengGP_QuantileEsts <-
  apply(spanishDichotomous_ZhengGP_QuantileSamples, 1, median)
spanishDichotomous_ZhengGP_QuantileIQRs <-
  apply(spanishDichotomous_ZhengGP_QuantileSamples, 1, IQR)

## Generate table results ---------------------------------------------

# McCarty et al. Data:

mccartyCount_ZhengGP_QuantileEsts
mccartyDichotomous_ZhengGP_QuantileEsts
mccartyCount_ZhengGP_QuantileIQRs
mccartyDichotomous_ZhengGP_QuantileIQRs

# Lubbers et al. Data:

spanishCount_ZhengGP_QuantileEsts
spanishDichotomous_ZhengGP_QuantileEsts
spanishCount_ZhengGP_QuantileIQRs
spanishDichotomous_ZhengGP_QuantileIQRs

#  Table 5 ---------------------------------------------

## Renormalize alphas for 3 rarest name subpopulations ---------------------------------------------

# Setup:

mccarty_sumrarenames <-
  sum(mccarty_props[order(mccarty_props[1:12])[1:3]])
mccarty_indrarenames <- order(mccarty_props[1:12])[1:3]

### Count ---------------------------------------------

# Compute renormalization constant using the 3 rarest name subpopulations:

mccartyCount_ZhengGP_Constant_RareNames <- c() # vector of constants

for (i in 1:nrow(mccartyCount_ZhengGP_ExtractBetas)) {
  mccartyCount_ZhengGP_Constant_RareNames[i] = log(sum(exp(mccartyCount_ZhengGP_ExtractBetas[i, mccarty_indrarenames])) / mccarty_sumrarenames)
}

# Renormalize alphas:

mccartyCount_ZhengGP_NormalizedAlphas_RareNames <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyCount_ZhengGP_ExtractAlphas)) {
  mccartyCount_ZhengGP_NormalizedAlphas_RareNames[i,] = mccartyCount_ZhengGP_ExtractAlphas[i,] + mccartyCount_ZhengGP_Constant_RareNames[i]
}

### Dichotomous ---------------------------------------------

# Compute renormalization constant using the 3 rarest name subpopulations:

mccartyDichotomous_ZhengGP_Constant_RareNames <-
  c() # vector of constants

for (i in 1:nrow(mccartyDichotomous_ZhengGP_ExtractBetas)) {
  mccartyDichotomous_ZhengGP_Constant_RareNames[i] = log(sum(exp(
    mccartyDichotomous_ZhengGP_ExtractBetas[i, mccarty_indrarenames]
  )) / mccarty_sumrarenames)
}

# Renormalize alphas:

mccartyDichotomous_ZhengGP_NormalizedAlphas_RareNames <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyDichotomous_ZhengGP_ExtractAlphas)) {
  mccartyDichotomous_ZhengGP_NormalizedAlphas_RareNames[i,] = mccartyDichotomous_ZhengGP_ExtractAlphas[i,] + mccartyDichotomous_ZhengGP_Constant_RareNames[i]
}

## Compute posterior means  ---------------------------------------------

# Count:

mccartyCount_ZhengGP_NormalizedExpAlphas_RareNames <-
  exp(mccartyCount_ZhengGP_NormalizedAlphas_RareNames)
mccartyCount_ZhengGP_NormalizedDegreeEsts_RareNames <-
  colMeans(mccartyCount_ZhengGP_NormalizedExpAlphas_RareNames)

# Dichotomous:

mccartyDichotomous_ZhengGP_NormalizedExpAlphas_RareNames <-
  exp(mccartyDichotomous_ZhengGP_NormalizedAlphas_RareNames)
mccartyDichotomous_ZhengGP_NormalizedDegreeEsts_RareNames <-
  colMeans(mccartyDichotomous_ZhengGP_NormalizedExpAlphas_RareNames)

## Generate table results ---------------------------------------------

summary(mccartyCount_ZhengGP_NormalizedDegreeEsts_RareNames)
summary(mccartyDichotomous_ZhengGP_NormalizedDegreeEsts_RareNames)

#  Table 6 ---------------------------------------------

## Renormalize alphas for all subpopulations ---------------------------------------------

### Count ---------------------------------------------

# Extract betas:

mccartyCount_ZhengBEM_ExtractBetas <-
  mccartyCount_ZhengBEM_Extract$beta

# Extract alphas:

mccartyCount_ZhengBEM_ExtractAlphas <-
  mccartyCount_ZhengBEM_Extract$alpha

# Compute renormalization constant:

mccartyCount_ZhengBEM_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyCount_ZhengBEM_ExtractBetas)) {
  mccartyCount_ZhengBEM_Constant[i] = log(sum(exp(mccartyCount_ZhengBEM_ExtractBetas[i, mccarty_indsubpops])) / mccarty_sumprops)
}

# Renormalize alphas:

mccartyCount_ZhengBEM_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyCount_ZhengBEM_ExtractAlphas)) {
  mccartyCount_ZhengBEM_NormalizedAlphas[i,] = mccartyCount_ZhengBEM_ExtractAlphas[i,] + mccartyCount_ZhengBEM_Constant[i]
}

### Dichotomous ---------------------------------------------

# Extract betas:

mccartyDichotomous_ZhengBEM_ExtractBetas <-
  mccartyDichotomous_ZhengBEM_Extract$beta

# Extract alphas:

mccartyDichotomous_ZhengBEM_ExtractAlphas <-
  mccartyDichotomous_ZhengBEM_Extract$alpha

# Compute renormalization constant:

mccartyDichotomous_ZhengBEM_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyDichotomous_ZhengBEM_ExtractBetas)) {
  mccartyDichotomous_ZhengBEM_Constant[i] = log(sum(exp(
    mccartyDichotomous_ZhengBEM_ExtractBetas[i, mccarty_indsubpops]
  )) / mccarty_sumprops)
}

# Renormalize alphas:

mccartyDichotomous_ZhengBEM_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyDichotomous_ZhengBEM_ExtractAlphas)) {
  mccartyDichotomous_ZhengBEM_NormalizedAlphas[i,] = mccartyDichotomous_ZhengBEM_ExtractAlphas[i,] + mccartyDichotomous_ZhengBEM_Constant[i]
}

## Generate table results ---------------------------------------------

# Count:

mccartyCount_ZhengBEM_NormalizedExpAlphas <-
  exp(mccartyCount_ZhengBEM_NormalizedAlphas)

mccartyCount_ZhengBEM_QuantileSamples <-
  apply(
    mccartyCount_ZhengBEM_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyCount_ZhengBEM_QuantileEsts <-
  apply(mccartyCount_ZhengBEM_QuantileSamples, 1, median)
mccartyCount_ZhengBEM_QuantileIQRs <-
  apply(mccartyCount_ZhengBEM_QuantileSamples, 1, IQR)

# Dichotomous:

mccartyDichotomous_ZhengBEM_NormalizedExpAlphas <-
  exp(mccartyDichotomous_ZhengBEM_NormalizedAlphas)

mccartyDichotomous_ZhengBEM_QuantileSamples <-
  apply(
    mccartyDichotomous_ZhengBEM_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyDichotomous_ZhengBEM_QuantileEsts <-
  apply(mccartyDichotomous_ZhengBEM_QuantileSamples, 1, median)
mccartyDichotomous_ZhengBEM_QuantileIQRs <-
  apply(mccartyDichotomous_ZhengBEM_QuantileSamples, 1, IQR)

mccartyCount_ZhengBEM_QuantileEsts
mccartyDichotomous_ZhengBEM_QuantileEsts
mccartyCount_ZhengBEM_QuantileIQRs
mccartyDichotomous_ZhengBEM_QuantileIQRs

#  Table 7 ---------------------------------------------

## Compute posterior means from model with count data ---------------------------------------------

mccartyCount_ZhengBEM_NormalizedExpAlphas <-
  exp(mccartyCount_ZhengBEM_NormalizedAlphas)
mccartyCount_ZhengBEM_NormalizedDegreeEsts <-
  colMeans(mccartyCount_ZhengBEM_NormalizedExpAlphas)

## Prepare data ---------------------------------------------

# Load dichotomous data:

mccartyDichotomous_dta <-
  readRDS(
    "Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyDichotomous_dta.rds"
  )

# Identify cases with missing data:

mccartyDichotomous_dta <- mccartyDichotomous_dta %>%
  replace_with_na_all(condition = ~ .x == -99)

mccartyMissingInd <- c()
for (i in 1:nrow(mccartyDichotomous_dta)) {
  if (anyNA(mccartyDichotomous_dta[i, ])) {
    mccartyMissingInd[i] = TRUE
  } else {
    mccartyMissingInd[i] = FALSE
  }
}

# Remove missing cases:

mccartyDichotomousComplete_dta <- na.omit(mccartyDichotomous_dta)

# Extract Nicole column:

zik_NicoleDichotomous <- mccartyDichotomousComplete_dta$n12a

# Count data network size estimates for complete cases:

mccartyCount_ZhengBEM_NormalizedDegreeEsts_Complete <-
  mccartyCount_ZhengBEM_NormalizedDegreeEsts[!mccartyMissingInd]

## Likelihood function ---------------------------------------------

ZhengBEM_Likelihood_fn <- function(theta, degrees, zik) {
  NLL <- c()
  for (i in 1:length(zik)) {
    if (zik[i] == 0) {
      NLL[i] = (degrees[i] * theta[1] * log(theta[2])) / (theta[2] - 1)
    } else if (zik[i] == 1) {
      NLL[i] = -log(1 - (1 / theta[2]) ^ (degrees[i] * theta[1] / (theta[2] - 1)))
    }
  }
  sum(NLL)
}

## Run optimization algorithm with different initial values ---------------------------------------------

init_values <- list(c(0.1, 10), c(0.5, 2), c(0.9, 15))

ZhengBEM_OptimResults <-
  map(
    init_values,
    ~ optim(
      .x,
      ZhengBEM_Likelihood_fn,
      degrees = mccartyCount_ZhengBEM_NormalizedDegreeEsts_Complete,
      zik = zik_NicoleDichotomous,
      method = "Nelder-Mead",
      hessian = FALSE
    )
  )

## Generate table results ---------------------------------------------

bk_omegak_pairs <- map(ZhengBEM_OptimResults, ~ .x$par)
s_bkomegak <-
  map_dbl(bk_omegak_pairs, ~ map2_dbl(.x[[1]], .x[[2]], ~ .x * log(.y) / (.y - 1)))
NLL <- map_dbl(ZhengBEM_OptimResults, ~ .x$value)

bk_omegak_pairs
s_bkomegak
NLL

#  Table 8 ---------------------------------------------

## Extract degree samples ---------------------------------------------

# Count:

mccartyCount_MaltielBEM_ExtractDegrees <-
  mccartyCount_MaltielBEM_Extract$d

# Dichotomous:

mccartyDichotomous_MaltielBEM_ExtractDegrees <-
  mccartyDichotomous_MaltielBEM_Extract$d

## Generate table results ---------------------------------------------

# Count:

mccartyCount_MaltielBEM_QuantileSamples <-
  apply(
    mccartyCount_MaltielBEM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyCount_MaltielBEM_QuantileEsts <-
  apply(mccartyCount_MaltielBEM_QuantileSamples, 1, median)
mccartyCount_MaltielBEM_QuantileIQRs <-
  apply(mccartyCount_MaltielBEM_QuantileSamples, 1, IQR)

# Dichotomous:

mccartyDichotomous_MaltielBEM_QuantileSamples <-
  apply(
    mccartyDichotomous_MaltielBEM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyDichotomous_MaltielBEM_QuantileEsts <-
  apply(mccartyDichotomous_MaltielBEM_QuantileSamples, 1, median)
mccartyDichotomous_MaltielBEM_QuantileIQRs <-
  apply(mccartyDichotomous_MaltielBEM_QuantileSamples, 1, IQR)

mccartyCount_MaltielBEM_QuantileEsts
mccartyDichotomous_MaltielBEM_QuantileEsts
mccartyCount_MaltielBEM_QuantileIQRs
mccartyDichotomous_MaltielBEM_QuantileIQRs

#  Table 9 ---------------------------------------------

## Load count data ---------------------------------------------

mccartyCount_dta <-
  readRDS("Dichotomous-ARD-Replication-Materials/Data/Main Text/mccartyCount_dta.rds")

mccartyCount_dta <- mccartyCount_dta %>%
  replace_with_na_all(condition = ~ .x == -99)

## Average number of Michaels known ---------------------------------------------

yik_MichaelCount <- mccartyCount_dta$n1a
mean_yik_michaels <- mean(yik_MichaelCount, na.rm = TRUE)

## Proportion of people knowing Michaels ---------------------------------------------

zik_MichaelDichotomous <- mccartyDichotomous_dta$n1a
prop_zik_michaels <- mean(zik_MichaelDichotomous, na.rm = TRUE)

## Extract parameters ---------------------------------------------

### Fraction in the population ---------------------------------------------

mk_michael <- mccarty_props[1]

### rhok parameter ---------------------------------------------

# Count:

mccartyCount_MaltielBEM_ExtractRho <-
  mccartyCount_MaltielBEM_Extract$rho

rhok_michaelCount <- mccartyCount_MaltielBEM_ExtractRho[, 1]

# Dichotomous:

mccartyDichotomous_MaltielBEM_ExtractRho <-
  mccartyDichotomous_MaltielBEM_Extract$rho

rhok_michaelDichotomous <-
  mccartyDichotomous_MaltielBEM_ExtractRho[, 1]

## Simulating prior pik distributions using mk and rhok (Using Equation 6) ---------------------------------------------

set.seed(1234)

### Count ---------------------------------------------

piksims_michaelCount <- matrix(data = NA, nrow = 4000, ncol = 1365)
for (i in 1:nrow(piksims_michaelCount)) {
  gamma_k <- mk_michael * (1 / rhok_michaelCount[i] - 1)
  eta_k <- (1 - mk_michael) * (1 / rhok_michaelCount[i] - 1)
  piksims_michaelCount[i, ] <-
    rbeta(ncol(piksims_michaelCount), gamma_k, eta_k)
}

### Dichotomous ---------------------------------------------

piksims_michaelDichotomous <-
  matrix(data = NA, nrow = 4000, ncol = 1365)
for (i in 1:nrow(piksims_michaelDichotomous)) {
  gamma_k <- mk_michael * (1 / rhok_michaelDichotomous[i] - 1)
  eta_k <- (1 - mk_michael) * (1 / rhok_michaelDichotomous[i] - 1)
  piksims_michaelDichotomous[i, ] <-
    rbeta(ncol(piksims_michaelDichotomous), gamma_k, eta_k)
}

## Simulating di distributions using the simulated pik ---------------------------------------------

### Count (Equation 12) ---------------------------------------------

disims_michaelCount <- matrix(data = NA, nrow = 4000, ncol = 1365)
for (i in 1:nrow(disims_michaelCount)) {
  for (k in 1:ncol(disims_michaelCount)) {
    disims_michaelCount[i, k] <-
      mean_yik_michaels / piksims_michaelCount[i, k]
  }
}

### Dichotomous (Equation 13) ---------------------------------------------

disims_michaelDichotomous <-
  matrix(data = NA, nrow = 4000, ncol = 1365)
for (i in 1:nrow(disims_michaelDichotomous)) {
  for (k in 1:ncol(disims_michaelDichotomous)) {
    disims_michaelDichotomous[i, k] <-
      log(1 - prop_zik_michaels) / log(1 - piksims_michaelDichotomous[i, k])
  }
}

## Compute estimates and IQRs ---------------------------------------------

### Count ---------------------------------------------

#### pik ---------------------------------------------

# Quantiles:

piksims_michaelCount_QuantileSamples <- apply(piksims_michaelCount,
                                              1,
                                              quantile,
                                              probs = c(0, 0.25, 0.50, 0.75, 1))

piksims_michaelCount_QuantileEsts <-
  apply(piksims_michaelCount_QuantileSamples, 1, median)
piksims_michaelCount_QuantileIQRs <-
  apply(piksims_michaelCount_QuantileSamples, 1, IQR)

# Mean:

piksims_michaelCount_MeanSamples <- apply(piksims_michaelCount,
                                          1,
                                          mean)

piksims_michaelCount_MeanEst <-
  median(piksims_michaelCount_MeanSamples)
piksims_michaelCount_MeanIQR <-
  IQR(piksims_michaelCount_MeanSamples)

#### di ---------------------------------------------

# Quantiles:

disims_michaelCount_QuantileSamples <- apply(disims_michaelCount,
                                             1,
                                             quantile,
                                             probs = c(0, 0.25, 0.50, 0.75, 1))

disims_michaelCount_QuantileEsts <-
  apply(disims_michaelCount_QuantileSamples, 1, median)
disims_michaelCount_QuantileIQRs <-
  apply(disims_michaelCount_QuantileSamples, 1, IQR)

# Mean:

disims_michaelCount_MeanSamples <- apply(disims_michaelCount,
                                         1,
                                         mean)

disims_michaelCount_MeanEst <-
  median(disims_michaelCount_MeanSamples)
disims_michaelCount_MeanIQR <- IQR(disims_michaelCount_MeanSamples)

### Dichotomous ---------------------------------------------

#### pik ---------------------------------------------

# Quantiles:

piksims_michaelDichotomous_QuantileSamples <-
  apply(piksims_michaelDichotomous,
        1,
        quantile,
        probs = c(0, 0.25, 0.50, 0.75, 1))

piksims_michaelDichotomous_QuantileEsts <-
  apply(piksims_michaelDichotomous_QuantileSamples, 1, median)
piksims_michaelDichotomous_QuantileIQRs <-
  apply(piksims_michaelDichotomous_QuantileSamples, 1, IQR)

# Mean:

piksims_michaelDichotomous_MeanSamples <-
  apply(piksims_michaelDichotomous,
        1,
        mean)

piksims_michaelDichotomous_MeanEst <-
  median(piksims_michaelDichotomous_MeanSamples)
piksims_michaelDichotomous_MeanIQR <-
  IQR(piksims_michaelDichotomous_MeanSamples)

#### di ---------------------------------------------

# Quantiles:

disims_michaelDichotomous_QuantileSamples <-
  apply(disims_michaelDichotomous,
        1,
        quantile,
        probs = c(0, 0.25, 0.50, 0.75, 1))

disims_michaelDichotomous_QuantileEsts <-
  apply(disims_michaelDichotomous_QuantileSamples, 1, median)
disims_michaelDichotomous_QuantileIQRs <-
  apply(disims_michaelDichotomous_QuantileSamples, 1, IQR)

# Mean:

disims_michaelDichotomous_MeanSamples <-
  apply(disims_michaelDichotomous,
        1,
        mean)

disims_michaelDichotomous_MeanEst <-
  median(disims_michaelDichotomous_MeanSamples)
disims_michaelDichotomous_MeanIQR <-
  IQR(disims_michaelDichotomous_MeanSamples)

## Generate table results ---------------------------------------------

### Panel A ---------------------------------------------

#### Count ---------------------------------------------

# Estimates:

pik_michaelCount_Ests <-
  c(
    piksims_michaelCount_QuantileEsts[1:3],
    piksims_michaelDichotomous_MeanEst,
    piksims_michaelCount_QuantileEsts[4:5]
  )

# IQRs:

pik_michaelCount_IQRs <-
  c(
    piksims_michaelCount_QuantileIQRs[1:3],
    piksims_michaelDichotomous_MeanIQR,
    piksims_michaelCount_QuantileIQRs[4:5]
  )

#### Dichotomous ---------------------------------------------

# Estimates:

pik_michaelDichotomous_Ests <-
  c(
    piksims_michaelDichotomous_QuantileEsts[1:3],
    piksims_michaelDichotomous_MeanEst,
    piksims_michaelDichotomous_QuantileEsts[4:5]
  )

# IQRs:

pik_michaelDichotomous_IQRs <-
  c(
    piksims_michaelDichotomous_QuantileIQRs[1:3],
    piksims_michaelDichotomous_MeanIQR,
    piksims_michaelDichotomous_QuantileIQRs[4:5]
  )

### Panel B ---------------------------------------------

#### Count ---------------------------------------------

# Estimates:

di_michaelCount_Ests <-
  c(
    disims_michaelCount_QuantileEsts[1:3],
    disims_michaelCount_MeanEst,
    disims_michaelCount_QuantileEsts[4:5]
  )

# IQRs:

di_michaelCount_IQRs <-
  c(
    disims_michaelCount_QuantileIQRs[1:3],
    disims_michaelCount_MeanIQR,
    disims_michaelCount_QuantileIQRs[4:5]
  )

#### Dichotomous ---------------------------------------------

# Estimates:

di_michaelDichotomous_Ests <-
  c(
    disims_michaelDichotomous_QuantileEsts[1:3],
    disims_michaelDichotomous_MeanEst,
    disims_michaelDichotomous_QuantileEsts[4:5]
  )

# IQRs:

di_michaelDichotomous_IQRs <-
  c(
    disims_michaelDichotomous_QuantileIQRs[1:3],
    disims_michaelDichotomous_MeanIQR,
    disims_michaelDichotomous_QuantileIQRs[4:5]
  )

#  Table 10 ---------------------------------------------

## Setup ---------------------------------------------

mccarty_firstnames <- names(mccarty_props)[1:12]
mccarty_firstnames_reordered <- c(
  "Anthony",
  "Christopher",
  "David",
  "James",
  "Michael",
  "Robert",
  "Christina",
  "Jacqueline",
  "Jennifer",
  "Kimberly",
  "Nicole",
  "Stephanie"
)
mccarty_firstnames_column <-
  rep(mccarty_firstnames_reordered, each = 2)

## Sex coefficients ---------------------------------------------

### Extract coefficient samples  ---------------------------------------------

#### Count ---------------------------------------------

mccartyCount_CovMod_ExtractSexCoeffs <-
  mccartyCount_CovMod_Extract$gamma_sex

mccartyCount_CovMod_ExtractSexCoeffs_df <-
  data.frame(mccartyCount_CovMod_ExtractSexCoeffs)

colnames(mccartyCount_CovMod_ExtractSexCoeffs_df) <-
  mccarty_firstnames

# Reorganize columns:

mccartyCount_CovMod_ExtractSexCoeffs_df <-
  mccartyCount_CovMod_ExtractSexCoeffs_df %>%
  dplyr::select(mccarty_firstnames_reordered)

#### Dichomotous ---------------------------------------------

mccartyDichotomous_CovMod_ExtractSexCoeffs <-
  mccartyDichotomous_CovMod_Extract$gamma_sex

mccartyDichotomous_CovMod_ExtractSexCoeffs_df <-
  data.frame(mccartyDichotomous_CovMod_ExtractSexCoeffs)

colnames(mccartyDichotomous_CovMod_ExtractSexCoeffs_df) <-
  mccarty_firstnames

# Reorganize columns:

mccartyDichotomous_CovMod_ExtractSexCoeffs_df <-
  mccartyDichotomous_CovMod_ExtractSexCoeffs_df %>%
  dplyr::select(mccarty_firstnames_reordered)

### Compute posterior means and credible intervals ---------------------------------------------

#### Count ---------------------------------------------

mccartyCount_CovMod_EstsCIsSexCoeffs <- c()
for (i in 1:ncol(mccartyCount_CovMod_ExtractSexCoeffs_df)) {
  mccartyCount_CovMod_EstsCIsSexCoeffs[2 * i - 1] <-
    formatC(round(mean(
      mccartyCount_CovMod_ExtractSexCoeffs_df[, i]
    ), 4), format = "f")
  mccartyCount_CovMod_EstsCIsSexCoeffs[2 * i] <- paste0("[",
                                                        formatC(round(
                                                          quantile(mccartyCount_CovMod_ExtractSexCoeffs_df[, i], 0.025),
                                                          4
                                                        ), format = "f"),
                                                        ", ",
                                                        formatC(round(
                                                          quantile(mccartyCount_CovMod_ExtractSexCoeffs_df[, i], 0.975),
                                                          4
                                                        ), format = "f"),
                                                        "]")
}

#### Dichotomous ---------------------------------------------

mccartyDichotomous_CovMod_EstsCIsSexCoeffs <- c()
for (i in 1:ncol(mccartyDichotomous_CovMod_ExtractSexCoeffs_df)) {
  mccartyDichotomous_CovMod_EstsCIsSexCoeffs[2 * i - 1] <-
    formatC(round(mean(
      mccartyDichotomous_CovMod_ExtractSexCoeffs_df[, i]
    ), 4), format = "f")
  mccartyDichotomous_CovMod_EstsCIsSexCoeffs[2 * i] <- paste0("[",
                                                              formatC(round(
                                                                quantile(mccartyDichotomous_CovMod_ExtractSexCoeffs_df[, i], 0.025),
                                                                4
                                                              ), format = "f"),
                                                              ", ",
                                                              formatC(round(
                                                                quantile(mccartyDichotomous_CovMod_ExtractSexCoeffs_df[, i], 0.975),
                                                                4
                                                              ), format = "f"),
                                                              "]")
}

## Age coefficients ---------------------------------------------

### Extract coefficient samples  ---------------------------------------------

#### Count ---------------------------------------------

mccartyCount_CovMod_ExtractAgeCoeffs <-
  mccartyCount_CovMod_Extract$gamma_age

mccartyCount_CovMod_ExtractAgeCoeffs_df <-
  data.frame(mccartyCount_CovMod_ExtractAgeCoeffs)

colnames(mccartyCount_CovMod_ExtractAgeCoeffs_df) <-
  mccarty_firstnames

# Reorganize columns:

mccartyCount_CovMod_ExtractAgeCoeffs_df <-
  mccartyCount_CovMod_ExtractAgeCoeffs_df %>%
  dplyr::select(mccarty_firstnames_reordered)

#### Dichomotous ---------------------------------------------

mccartyDichotomous_CovMod_ExtractAgeCoeffs <-
  mccartyDichotomous_CovMod_Extract$gamma_age

mccartyDichotomous_CovMod_ExtractAgeCoeffs_df <-
  data.frame(mccartyDichotomous_CovMod_ExtractAgeCoeffs)

colnames(mccartyDichotomous_CovMod_ExtractAgeCoeffs_df) <-
  mccarty_firstnames

# Reorganize columns:

mccartyDichotomous_CovMod_ExtractAgeCoeffs_df <-
  mccartyDichotomous_CovMod_ExtractAgeCoeffs_df %>%
  dplyr::select(mccarty_firstnames_reordered)

### Compute posterior means and credible intervals ---------------------------------------------

#### Count ---------------------------------------------

mccartyCount_CovMod_EstsCIsAgeCoeffs <- c()
for (i in 1:ncol(mccartyCount_CovMod_ExtractAgeCoeffs_df)) {
  mccartyCount_CovMod_EstsCIsAgeCoeffs[2 * i - 1] <-
    formatC(round(mean(
      mccartyCount_CovMod_ExtractAgeCoeffs_df[, i]
    ), 4), format = "f")
  mccartyCount_CovMod_EstsCIsAgeCoeffs[2 * i] <- paste0("[",
                                                        formatC(round(
                                                          quantile(mccartyCount_CovMod_ExtractAgeCoeffs_df[, i], 0.025),
                                                          4
                                                        ), format = "f"),
                                                        ", ",
                                                        formatC(round(
                                                          quantile(mccartyCount_CovMod_ExtractAgeCoeffs_df[, i], 0.975),
                                                          4
                                                        ), format = "f"),
                                                        "]")
}

#### Dichotomous ---------------------------------------------

mccartyDichotomous_CovMod_EstsCIsAgeCoeffs <- c()
for (i in 1:ncol(mccartyDichotomous_CovMod_ExtractAgeCoeffs_df)) {
  mccartyDichotomous_CovMod_EstsCIsAgeCoeffs[2 * i - 1] <-
    formatC(round(mean(
      mccartyDichotomous_CovMod_ExtractAgeCoeffs_df[, i]
    ), 4), format = "f")
  mccartyDichotomous_CovMod_EstsCIsAgeCoeffs[2 * i] <- paste0("[",
                                                              formatC(round(
                                                                quantile(mccartyDichotomous_CovMod_ExtractAgeCoeffs_df[, i], 0.025),
                                                                4
                                                              ), format = "f"),
                                                              ", ",
                                                              formatC(round(
                                                                quantile(mccartyDichotomous_CovMod_ExtractAgeCoeffs_df[, i], 0.975),
                                                                4
                                                              ), format = "f"),
                                                              "]")
}

## Prepare data frame with results ---------------------------------------------

mccarty_CovMod_EstsCIsCoeffs <-
  data.frame(
    name = mccarty_firstnames_column,
    CountSexCoeffs = mccartyCount_CovMod_EstsCIsSexCoeffs,
    DichotomousSexCoeffs = mccartyDichotomous_CovMod_EstsCIsSexCoeffs,
    CountAgeCoeffs = mccartyCount_CovMod_EstsCIsAgeCoeffs,
    DichotomousAgeCoeffs = mccartyDichotomous_CovMod_EstsCIsAgeCoeffs
  )

#  Table 11 ---------------------------------------------

## Renormalize alphas for the 12 name subpopulations ---------------------------------------------

# Setup:

mccarty_sumfirstnames <-
  sum(mccarty_props[1:12])
mccarty_indfirstnames <- 1:12

### Covariate-based ---------------------------------------------

# Extract betas:

mccartyDichotomous_CovMod_ExtractBetas <-
  mccartyDichotomous_CovMod_Extract$beta

# Extract alphas:

mccartyDichotomous_CovMod_ExtractAlphas <-
  mccartyDichotomous_CovMod_Extract$alpha

# Compute renormalization constant:

mccartyDichotomous_CovMod_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyDichotomous_CovMod_ExtractBetas)) {
  mccartyDichotomous_CovMod_Constant[i] = log(sum(exp(
    mccartyDichotomous_CovMod_ExtractBetas[i, mccarty_indfirstnames]
  )) / mccarty_sumfirstnames)
}

# Renormalize alphas:

mccartyDichotomous_CovMod_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1340)

for (i in 1:nrow(mccartyDichotomous_CovMod_ExtractAlphas)) {
  mccartyDichotomous_CovMod_NormalizedAlphas[i,] = mccartyDichotomous_CovMod_ExtractAlphas[i,] + mccartyDichotomous_CovMod_Constant[i]
}

### Group prevalence ---------------------------------------------

# Extract betas:

mccartyDichotomousReduced_ZhengGP_ExtractBetas <-
  mccartyDichotomousReduced_ZhengGP_Extract$beta

# Extract alphas:

mccartyDichotomousReduced_ZhengGP_ExtractAlphas <-
  mccartyDichotomousReduced_ZhengGP_Extract$alpha

# Compute renormalization constant:

mccartyDichotomousReduced_ZhengGP_Constant <-
  c() # vector of constants

for (i in 1:nrow(mccartyDichotomousReduced_ZhengGP_ExtractBetas)) {
  mccartyDichotomousReduced_ZhengGP_Constant[i] = log(sum(exp(
    mccartyDichotomousReduced_ZhengGP_ExtractBetas[i, mccarty_indfirstnames]
  )) / mccarty_sumfirstnames)
}

# Renormalize alphas:

mccartyDichotomousReduced_ZhengGP_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1340)

for (i in 1:nrow(mccartyDichotomousReduced_ZhengGP_ExtractAlphas)) {
  mccartyDichotomousReduced_ZhengGP_NormalizedAlphas[i,] = mccartyDichotomousReduced_ZhengGP_ExtractAlphas[i,] + mccartyDichotomousReduced_ZhengGP_Constant[i]
}

## Compute posterior means  ---------------------------------------------

# Covariate-based:

mccartyDichotomous_CovMod_NormalizedExpAlphas <-
  exp(mccartyDichotomous_CovMod_NormalizedAlphas)
mccartyDichotomous_CovMod_NormalizedDegreeEsts <-
  colMeans(mccartyDichotomous_CovMod_NormalizedExpAlphas)

# Group prevalence:

mccartyDichotomousReduced_ZhengGP_NormalizedExpAlphas <-
  exp(mccartyDichotomousReduced_ZhengGP_NormalizedAlphas)
mccartyDichotomousReduced_ZhengGP_NormalizedDegreeEsts <-
  colMeans(mccartyDichotomousReduced_ZhengGP_NormalizedExpAlphas)

## Generate table results ---------------------------------------------

summary(mccartyDichotomous_CovMod_NormalizedDegreeEsts)
summary(mccartyDichotomousReduced_ZhengGP_NormalizedDegreeEsts)

#  Figure 4 ---------------------------------------------

## Create plot  ---------------------------------------------

# Joint covariate-based and group prevalence estimates into a sigle data frame:

mccarty_CovGP_DegreeEsts_df <-
  data.frame(networkdegree_Cov = mccartyDichotomous_CovMod_NormalizedDegreeEsts,
             networkdegree_GP = mccartyDichotomousReduced_ZhengGP_NormalizedDegreeEsts)

# Lengthen the data:

mccarty_CovGP_DegreeEsts_df_long <-
  mccarty_CovGP_DegreeEsts_df %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "model"),
    names_sep = "_"
  )

mccarty_CovGP_DegreeEsts_df_long$model <-
  factor(
    mccarty_CovGP_DegreeEsts_df_long$model,
    levels = c("Cov", "GP"),
    labels = c("Covariate-Based", "Group Prevalence")
  )

# Create kernel density plot:

mccarty_CovGP_DegreeEsts_KDP <-
  ggplot(mccarty_CovGP_DegreeEsts_df_long,
         aes(networkdegree, color = model, linetype = model)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 1100)) +
  labs(color = "",
       linetype = "",
       x = "Posterior Mean Network Size") +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  scale_linetype_manual(values = c("Covariate-Based" = "solid", "Group Prevalence" = "dashed")) +
  theme_classic() +
  theme(legend.text = element_text(size = 10))

#  Figure 5 ---------------------------------------------

## Compute posterior means for barrier effect parameters ---------------------------------------------

### rho ---------------------------------------------

# Count:

mccartyCount_MaltielBEM_RhoEsts <-
  colMeans(mccartyCount_MaltielBEM_ExtractRho)

# 0/1/2+

mccartyTrichotomous_MaltielBEM_ExtractRho <-
  mccartyTrichotomous_MaltielBEM_Extract$rho

mccartyTrichotomous_MaltielBEM_RhoEsts <-
  colMeans(mccartyTrichotomous_MaltielBEM_ExtractRho)

### omega ---------------------------------------------

# Count:

mccartyCount_ZhengBEM_ExtractInvomega <-
  mccartyCount_ZhengBEM_Extract$inv_omega

mccartyCount_ZhengBEM_ExtractOmega <-
  1 / mccartyCount_ZhengBEM_ExtractInvomega

mccartyCount_ZhengBEM_OmegaEsts <-
  colMeans(mccartyCount_ZhengBEM_ExtractOmega)

# 0/1/2+:

mccartyTrichotomous_ZhengBEM_ExtractInvomega <-
  mccartyTrichotomous_ZhengBEM_Extract$inv_omega

mccartyTrichotomous_ZhengBEM_ExtractOmega <-
  1 / mccartyTrichotomous_ZhengBEM_ExtractInvomega

mccartyTrichotomous_ZhengBEM_OmegaEsts <-
  colMeans(mccartyTrichotomous_ZhengBEM_ExtractOmega)

## Setup ---------------------------------------------

mccarty_subpopnames_adj <- names(McCarty$known)
mccarty_subpopnames_adj[26] <- "Male in prison"
mccarty_subpopnames_adj[20] <- "Jaycees*"

## Construct credible intervals for rho parameters ---------------------------------------------

# Count:

mccartyCount_MaltielBEM_RhoCIs <- matrix(nrow = 29, ncol = 2)
for (i in 1:ncol(mccartyCount_MaltielBEM_ExtractRho)) {
  mccartyCount_MaltielBEM_RhoCIs[i, 1] <-
    quantile(mccartyCount_MaltielBEM_ExtractRho[, i], 0.025)
  mccartyCount_MaltielBEM_RhoCIs[i, 2] <-
    quantile(mccartyCount_MaltielBEM_ExtractRho[, i], 0.975)
}

# 0/1/2+:

mccartyTrichotomous_MaltielBEM_RhoCIs <- matrix(nrow = 29, ncol = 2)
for (i in 1:ncol(mccartyTrichotomous_MaltielBEM_ExtractRho)) {
  mccartyTrichotomous_MaltielBEM_RhoCIs[i, 1] <-
    quantile(mccartyTrichotomous_MaltielBEM_ExtractRho[, i], 0.025)
  mccartyTrichotomous_MaltielBEM_RhoCIs[i, 2] <-
    quantile(mccartyTrichotomous_MaltielBEM_ExtractRho[, i], 0.975)
}

## Construct credible intervals for omega parameters ---------------------------------------------

# Count:

mccartyCount_ZhengBEM_OmegaCIs <- matrix(nrow = 29, ncol = 2)
for (i in 1:ncol(mccartyCount_ZhengBEM_ExtractOmega)) {
  mccartyCount_ZhengBEM_OmegaCIs[i, 1] <-
    quantile(mccartyCount_ZhengBEM_ExtractOmega[, i], 0.025)
  mccartyCount_ZhengBEM_OmegaCIs[i, 2] <-
    quantile(mccartyCount_ZhengBEM_ExtractOmega[, i], 0.975)
}

# 0/1/2+:

mccartyTrichotomous_ZhengBEM_OmegaCIs <- matrix(nrow = 29, ncol = 2)
for (i in 1:ncol(mccartyTrichotomous_ZhengBEM_ExtractOmega)) {
  mccartyTrichotomous_ZhengBEM_OmegaCIs[i, 1] <-
    quantile(mccartyTrichotomous_ZhengBEM_ExtractOmega[, i], 0.025)
  mccartyTrichotomous_ZhengBEM_OmegaCIs[i, 2] <-
    quantile(mccartyTrichotomous_ZhengBEM_ExtractOmega[, i], 0.975)
}

## Create a data frame with the information above ---------------------------------------------

mccarty_CountTrichotomous_Segregation_df <-
  data.frame(
    SubPop = mccarty_subpopnames_adj,
    Count_Rho_Estimate = mccartyCount_MaltielBEM_RhoEsts,
    Trichotomous_Rho_Estimate = mccartyTrichotomous_MaltielBEM_RhoEsts,
    Count_Omega_Estimate = mccartyCount_ZhengBEM_OmegaEsts,
    Trichotomous_Omega_Estimate = mccartyTrichotomous_ZhengBEM_OmegaEsts,
    Count_Rho_lb = mccartyCount_MaltielBEM_RhoCIs[, 1],
    Trichotomous_Rho_lb = mccartyTrichotomous_MaltielBEM_RhoCIs[, 1],
    Count_Rho_ub = mccartyCount_MaltielBEM_RhoCIs[, 2],
    Trichotomous_Rho_ub = mccartyTrichotomous_MaltielBEM_RhoCIs[, 2],
    Count_Omega_lb = mccartyCount_ZhengBEM_OmegaCIs[, 1],
    Trichotomous_Omega_lb = mccartyTrichotomous_ZhengBEM_OmegaCIs[, 1],
    Count_Omega_ub = mccartyCount_ZhengBEM_OmegaCIs[, 2],
    Trichotomous_Omega_ub = mccartyTrichotomous_ZhengBEM_OmegaCIs[, 2]
  )

# Lengthen data:

mccarty_CountTrichotomous_Segregation_df_long <-
  mccarty_CountTrichotomous_Segregation_df %>%
  pivot_longer(
    cols = -SubPop,
    names_to = c("DataType", "ParameterType", ".value"),
    names_pattern = "(.*)_(.*)_(.*)"
  )

# Handle the Jaycees value for the omega parameter (see Figure 5 footnote):

mccarty_CountTrichotomous_Segregation_df_long[79:80, 4:6] <- NA
mccarty_CountTrichotomous_Segregation_df_long[77:78, 1] <- "Jaycees"
mccarty_subpopnames_JayceesNew <-
  c(mccarty_subpopnames_adj[1:20],
    "Jaycees",
    mccarty_subpopnames_adj[21:29])

# Change factor variables:

mccarty_CountTrichotomous_Segregation_df_long$DataType <-
  factor(
    mccarty_CountTrichotomous_Segregation_df_long$DataType,
    levels = c("Count", "Trichotomous"),
    labels = c("Count", "0/1/2+")
  )

mccarty_CountTrichotomous_Segregation_df_long$SubPop <-
  factor(
    mccarty_CountTrichotomous_Segregation_df_long$SubPop,
    levels = rev(mccarty_subpopnames_JayceesNew)
  )

mccarty_CountTrichotomous_Segregation_df_long$ParameterType <-
  factor(
    mccarty_CountTrichotomous_Segregation_df_long$ParameterType,
    labels = c(
      "omega~(Zheng~et~al.~Barrier~Effects~Model)",
      "rho~(Maltiel~et~al.~Barrier~Effects~Model)"
    )
  )

## Create plot ---------------------------------------------

pd <- position_dodge(width = 0.5)
Segregation_CompareCountTrichotomous <-
  ggplot(
    mccarty_CountTrichotomous_Segregation_df_long,
    aes(x = Estimate, y = SubPop, color = DataType)
  ) +
  geom_point(
    aes(shape = DataType),
    size = 1.2,
    position = pd,
    na.rm = TRUE
  ) +
  labs(
    x = "",
    y = "Subpopulation",
    color = "",
    shape = ""
  ) +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  scale_x_continuous(breaks = scales::pretty_breaks(7)) +
  geom_errorbar(
    aes(xmin = lb , xmax = ub),
    width = 0,
    position = pd,
    linewidth = 0.2
  ) +
  facet_wrap( ~ ParameterType,
              labeller = label_parsed,
              scales = "free") +
  theme_classic(base_size = 10)

# Table 12 ---------------------------------------------

# We have already produced results for count and dichotomous data
# in Tables 6 and 8 above

## Zheng et al. ---------------------------------------------

### 0/1/2+ ---------------------------------------------

# Extract betas:

mccartyTrichotomous_ZhengBEM_ExtractBetas <-
  mccartyTrichotomous_ZhengBEM_Extract$beta

# Extract alphas:

mccartyTrichotomous_ZhengBEM_ExtractAlphas <-
  mccartyTrichotomous_ZhengBEM_Extract$alpha

# Compute renormalization constant:

mccartyTrichotomous_ZhengBEM_Constant <- c() # vector of constants

for (i in 1:nrow(mccartyTrichotomous_ZhengBEM_ExtractBetas)) {
  mccartyTrichotomous_ZhengBEM_Constant[i] = log(sum(exp(
    mccartyTrichotomous_ZhengBEM_ExtractBetas[i, mccarty_indsubpops]
  )) / mccarty_sumprops)
}

# Renormalize alphas:

mccartyTrichotomous_ZhengBEM_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 1365)

for (i in 1:nrow(mccartyTrichotomous_ZhengBEM_ExtractAlphas)) {
  mccartyTrichotomous_ZhengBEM_NormalizedAlphas[i,] = mccartyTrichotomous_ZhengBEM_ExtractAlphas[i,] + mccartyTrichotomous_ZhengBEM_Constant[i]
}

# Generate table results:

mccartyTrichotomous_ZhengBEM_NormalizedExpAlphas <-
  exp(mccartyTrichotomous_ZhengBEM_NormalizedAlphas)

mccartyTrichotomous_ZhengBEM_QuantileSamples <-
  apply(
    mccartyTrichotomous_ZhengBEM_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyTrichotomous_ZhengBEM_QuantileEsts <-
  apply(mccartyTrichotomous_ZhengBEM_QuantileSamples, 1, median)
mccartyTrichotomous_ZhengBEM_QuantileIQRs <-
  apply(mccartyTrichotomous_ZhengBEM_QuantileSamples, 1, IQR)

mccartyCount_ZhengBEM_QuantileEsts
mccartyDichotomous_ZhengBEM_QuantileEsts
mccartyTrichotomous_ZhengBEM_QuantileEsts
mccartyCount_ZhengBEM_QuantileIQRs
mccartyDichotomous_ZhengBEM_QuantileIQRs
mccartyTrichotomous_ZhengBEM_QuantileIQRs

## Maltiel et al. ---------------------------------------------

### 0/1/2+ ---------------------------------------------

# Extract degree samples:

mccartyTrichotomous_MaltielBEM_ExtractDegrees <-
  mccartyTrichotomous_MaltielBEM_Extract$d

# Generate table results:

mccartyTrichotomous_MaltielBEM_QuantileSamples <-
  apply(
    mccartyTrichotomous_MaltielBEM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

mccartyTrichotomous_MaltielBEM_QuantileEsts <-
  apply(mccartyTrichotomous_MaltielBEM_QuantileSamples, 1, median)
mccartyTrichotomous_MaltielBEM_QuantileIQRs <-
  apply(mccartyTrichotomous_MaltielBEM_QuantileSamples, 1, IQR)

mccartyCount_MaltielBEM_QuantileEsts
mccartyDichotomous_MaltielBEM_QuantileEsts
mccartyTrichotomous_MaltielBEM_QuantileEsts
mccartyCount_MaltielBEM_QuantileIQRs
mccartyDichotomous_MaltielBEM_QuantileIQRs
mccartyTrichotomous_MaltielBEM_QuantileIQRs
