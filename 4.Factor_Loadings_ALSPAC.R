library(lavaan)
library(dplyr)
library(purrr)
library(mice)
library(readr)

# =========================
# LOAD DATA
# =========================

## Epigenetic age sample
df_epi <- readRDS(
  "...rds"
)

## Brain age sample
df_brain <- readRDS(
  "...rds"
)

## Cognition sample
df_cognition <- readRDS(
  "...rds"
)

df_cognition <- complete(df_cognition, action = "all")
df_brain <- complete(df_brain, action = "all")
df_epi <- complete(df_epi, action = "all")

# =========================
# BUILT ENVIRONMENT VARIABLES
# =========================

built_env_vars <- c(
  "building_density_age12child",
  "connectivity_density_age12child",
  "facility_richness_age12child",
  "greenspace_age12child",
  "pop_density_age12child",
  "walkability_age12child"
)

# =========================
# CFA MODEL
# =========================

cfa_model <- '
  BuiltEnv =~ building_density_age12child +
              connectivity_density_age12child +
              facility_richness_age12child +
              greenspace_age12child +
              pop_density_age12child +
              walkability_age12child
'

# =========================
# FUNCTION TO EXTRACT LOADINGS
# =========================

extract_loadings <- function(df_list, sample_name) {
  
  map_dfr(seq_along(df_list), function(i) {
    
    df <- df_list[[i]]
    
    # standardize variables
    df[built_env_vars] <- scale(df[built_env_vars])
    
    # fit CFA
    fit <- cfa(
      cfa_model,
      data = df,
      estimator = "MLR"
    )
    
    # extract standardized loadings
    standardizedSolution(fit) %>%
      filter(op == "=~") %>%
      mutate(
        sample = sample_name,
        imputation = i
      ) %>%
      select(
        sample,
        imputation,
        indicator = rhs,
        loading = est.std,
        se,
        z,
        pvalue
      )
  })
}

# =========================
# RUN FOR EACH SAMPLE
# =========================

epi_loadings <- extract_loadings(df_epi, "Epigenetic")

brain_loadings <- extract_loadings(df_brain, "Brain")

cognition_loadings <- extract_loadings(df_cognition, "Cognition")

# combine all results
all_loadings <- bind_rows(
  epi_loadings,
  brain_loadings,
  cognition_loadings
)

# =========================
# VIEW RESULTS
# =========================

print(all_loadings)

# =========================
# AVERAGE LOADINGS ACROSS IMPUTATIONS
# =========================

mean_loadings <- all_loadings %>%
  group_by(sample, indicator) %>%
  summarise(
    mean_loading = mean(loading, na.rm = TRUE),
    sd_loading = sd(loading, na.rm = TRUE),
    .groups = "drop"
  )

print(mean_loadings)

write_csv(
  mean_loadings,
  "...csv"
)
