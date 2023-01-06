#  "Uses and Limitations of Dichotomous Aggregate Relational Data"
#  R Script for replicating results in the article
#
#  Author: Derick S.Baum <derick_baum@g.harvard.edu>

################################ Appendix ################################

#  Load packages ---------------------------------------------

library(tidyverse)
library(NSUM)
library(ggrepel)
library(naniar)
library(scales)

#  Figure A1 ---------------------------------------------

## Extract degree samples  ---------------------------------------------

# Count:

spanishCount_MaltielRDM_ExtractDegrees <-
  spanishCount_MaltielRDM_Extract$d

# Dichotomous:

spanishDichotomous_MaltielRDM_ExtractDegrees <-
  spanishDichotomous_MaltielRDM_Extract$d

## Compute posterior means  ---------------------------------------------

# Count:

spanishCount_MaltielRDM_DegreeEsts <-
  colMeans(spanishCount_MaltielRDM_ExtractDegrees)

# Dichotomous:

spanishDichotomous_MaltielRDM_DegreeEsts <-
  colMeans(spanishDichotomous_MaltielRDM_ExtractDegrees)

## Create plot  ---------------------------------------------

# Joint count and dichotomous estimates into a sigle data frame:

spanish_MaltielRDM_DegreeEsts_df <-
  data.frame(networkdegree_Count = spanishCount_MaltielRDM_DegreeEsts,
             networkdegree_Dichotomous = spanishDichotomous_MaltielRDM_DegreeEsts)

# Lengthen the data:

spanish_MaltielRDM_DegreeEsts_df_long <-
  spanish_MaltielRDM_DegreeEsts_df %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "data_type"),
    names_sep = "_"
  )

spanish_MaltielRDM_DegreeEsts_df_long$data_type <-
  factor(spanish_MaltielRDM_DegreeEsts_df_long$data_type,
         levels = c("Count", "Dichotomous"))

# Create kernel density plot:

spanish_MaltielRDM_DegreeEsts_KDP <-
  ggplot(
    spanish_MaltielRDM_DegreeEsts_df_long,
    aes(networkdegree, color = data_type, linetype = data_type)
  ) +
  geom_density() +
  coord_cartesian(xlim = c(0, 3000)) +
  labs(color = "",
       linetype = "",
       x = "Posterior Mean Network Size") +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  scale_linetype_manual(values = c("Count" = "solid", "Dichotomous" = "dashed")) +
  theme_classic() +
  theme(legend.text = element_text(size = 10))

#  Table A1 ---------------------------------------------

summary(spanishCount_MaltielRDM_DegreeEsts)
summary(spanishDichotomous_MaltielRDM_DegreeEsts)

#  Table A2 ---------------------------------------------

# Count:

spanishCount_MaltielRDM_QuantileSamples <-
  apply(
    spanishCount_MaltielRDM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishCount_MaltielRDM_QuantileEsts <-
  apply(spanishCount_MaltielRDM_QuantileSamples, 1, median)
spanishCount_MaltielRDM_QuantileIQRs <-
  apply(spanishCount_MaltielRDM_QuantileSamples, 1, IQR)

# Dichotomous:

spanishDichotomous_MaltielRDM_QuantileSamples <-
  apply(
    spanishDichotomous_MaltielRDM_ExtractDegrees,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishDichotomous_MaltielRDM_QuantileEsts <-
  apply(spanishDichotomous_MaltielRDM_QuantileSamples, 1, median)
spanishDichotomous_MaltielRDM_QuantileIQRs <-
  apply(spanishDichotomous_MaltielRDM_QuantileSamples, 1, IQR)

spanishCount_MaltielRDM_QuantileEsts
spanishDichotomous_MaltielRDM_QuantileEsts
spanishCount_MaltielRDM_QuantileIQRs
spanishDichotomous_MaltielRDM_QuantileIQRs

#  Figure A3 ---------------------------------------------

## Renormalize alphas for all subpopulations ---------------------------------------------

# Setup:

spanish_subpops <-
  readRDS("Dichotomous-ARD-Replication-Materials/Data/Appendix/spanishProps.rds")
spanish_props <- head(spanish_subpops$prop, -1)

spanish_sumprops <- sum(spanish_props)
spanish_indsubpops <- 1:14

### Count ---------------------------------------------

# Extract betas:

spanishCount_ZhengBEM_ExtractBetas <-
  spanishCount_ZhengBEM_Extract$beta

# Extract alphas:

spanishCount_ZhengBEM_ExtractAlphas <-
  spanishCount_ZhengBEM_Extract$alpha

# Compute renormalization constant:

spanishCount_ZhengBEM_Constant <- c() # vector of constants

for (i in 1:nrow(spanishCount_ZhengBEM_ExtractBetas)) {
  spanishCount_ZhengBEM_Constant[i] = log(sum(exp(spanishCount_ZhengBEM_ExtractBetas[i, spanish_indsubpops])) / spanish_sumprops)
}

# Renormalize alphas:

spanishCount_ZhengBEM_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 2371)

for (i in 1:nrow(spanishCount_ZhengBEM_ExtractAlphas)) {
  spanishCount_ZhengBEM_NormalizedAlphas[i,] = spanishCount_ZhengBEM_ExtractAlphas[i,] + spanishCount_ZhengBEM_Constant[i]
}

### Dichotomous ---------------------------------------------

# Extract betas:

spanishDichotomous_ZhengBEM_ExtractBetas <-
  spanishDichotomous_ZhengBEM_Extract$beta

# Extract alphas:

spanishDichotomous_ZhengBEM_ExtractAlphas <-
  spanishDichotomous_ZhengBEM_Extract$alpha

# Compute renormalization constant:

spanishDichotomous_ZhengBEM_Constant <- c() # vector of constants

for (i in 1:nrow(spanishDichotomous_ZhengBEM_ExtractBetas)) {
  spanishDichotomous_ZhengBEM_Constant[i] = log(sum(exp(spanishDichotomous_ZhengBEM_ExtractBetas[i, spanish_indsubpops])) / spanish_sumprops)
}

# Renormalize alphas:

spanishDichotomous_ZhengBEM_NormalizedAlphas <-
  matrix(NA, nrow = 4000, ncol = 2371)

for (i in 1:nrow(spanishDichotomous_ZhengBEM_ExtractAlphas)) {
  spanishDichotomous_ZhengBEM_NormalizedAlphas[i,] = spanishDichotomous_ZhengBEM_ExtractAlphas[i,] + spanishDichotomous_ZhengBEM_Constant[i]
}

## Generate table results ---------------------------------------------

# Count:

spanishCount_ZhengBEM_NormalizedExpAlphas <-
  exp(spanishCount_ZhengBEM_NormalizedAlphas)

spanishCount_ZhengBEM_QuantileSamples <-
  apply(
    spanishCount_ZhengBEM_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishCount_ZhengBEM_QuantileEsts <-
  apply(spanishCount_ZhengBEM_QuantileSamples, 1, median)
spanishCount_ZhengBEM_QuantileIQRs <-
  apply(spanishCount_ZhengBEM_QuantileSamples, 1, IQR)

# Dichotomous:

spanishDichotomous_ZhengBEM_NormalizedExpAlphas <-
  exp(spanishDichotomous_ZhengBEM_NormalizedAlphas)

spanishDichotomous_ZhengBEM_QuantileSamples <-
  apply(
    spanishDichotomous_ZhengBEM_NormalizedExpAlphas,
    1,
    quantile,
    probs = c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
  )

spanishDichotomous_ZhengBEM_QuantileEsts <-
  apply(spanishDichotomous_ZhengBEM_QuantileSamples, 1, median)
spanishDichotomous_ZhengBEM_QuantileIQRs <-
  apply(spanishDichotomous_ZhengBEM_QuantileSamples, 1, IQR)

spanishCount_ZhengBEM_QuantileEsts
spanishDichotomous_ZhengBEM_QuantileEsts
spanishCount_ZhengBEM_QuantileIQRs
spanishDichotomous_ZhengBEM_QuantileIQRs

#  Figure B1 ---------------------------------------------

## Extract degree samples  ---------------------------------------------

# Count:

mccartyCountReduced_MaltielRDM_ExtractDegrees <-
  mccartyCountReduced_MaltielRDM_Extract$d

# Dichotomous:

mccartyDichotomousReduced_MaltielRDM_ExtractDegrees <-
  mccartyDichotomousReduced_MaltielRDM_Extract$d

## Compute posterior means  ---------------------------------------------

# Count:

mccartyCountReduced_MaltielRDM_DegreeEsts <-
  colMeans(mccartyCountReduced_MaltielRDM_ExtractDegrees)

# Dichotomous:

mccartyDichotomousReduced_MaltielRDM_DegreeEsts <-
  colMeans(mccartyDichotomousReduced_MaltielRDM_ExtractDegrees)

## Create plot  ---------------------------------------------

# Joint count and dichotomous estimates into a sigle data frame:

mccartyReduced_MaltielRDM_DegreeEsts_df <-
  data.frame(networkdegree_Count = mccartyCountReduced_MaltielRDM_DegreeEsts,
             networkdegree_Dichotomous = mccartyDichotomousReduced_MaltielRDM_DegreeEsts)

# Lengthen the data:

mccartyReduced_MaltielRDM_DegreeEsts_df_long <-
  mccartyReduced_MaltielRDM_DegreeEsts_df %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "data_type"),
    names_sep = "_"
  )

mccartyReduced_MaltielRDM_DegreeEsts_df_long$data_type <-
  factor(mccartyReduced_MaltielRDM_DegreeEsts_df_long$data_type,
         levels = c("Count", "Dichotomous"))

# Create kernel density plot:

mccartyReduced_MaltielRDM_DegreeEsts_KDP <-
  ggplot(
    mccartyReduced_MaltielRDM_DegreeEsts_df_long,
    aes(networkdegree, color = data_type, linetype = data_type)
  ) +
  geom_density() +
  coord_cartesian(xlim = c(0, 2500)) +
  labs(color = "",
       linetype = "",
       x = "Posterior Mean Network Size") +
  scale_color_manual(values = c("#00239CFF", "#E10600FF")) +
  scale_linetype_manual(values = c("Count" = "solid", "Dichotomous" = "dashed")) +
  theme_classic() +
  theme(legend.text = element_text(size = 10))

#  Table B1 ---------------------------------------------

summary(mccartyCountReduced_MaltielRDM_DegreeEsts)
summary(mccartyDichotomousReduced_MaltielRDM_DegreeEsts)
