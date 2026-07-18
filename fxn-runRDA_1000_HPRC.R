# Function to run highly robust RDA and ANOVA analysis for Gulf data
library(tidyverse)
library(compositions)
library(vegan)

load(file = "/scratch/group/hu-lab/GoM-amplicon-analysis/GoM-r-outputs/RDA_Robjects_07012026.RData", verbose = TRUE)


set.seed(10234)

# Function to create dataframes from output
run_rda <- function(mod_metadata, mod_asv_long, 
                    DATASET, features = c("include", "exclude")){
  # Select categorical data
  # mod_metadata <- surf_table
  # mod_asv_long <- surf_asv
  if (features == "include") {
    categorical <- mod_metadata %>% 
      select(STN_NISKIN, WATER_MASS, OilRig_Count, DAY_NIGHT, ends_with("_THRESHOLD")) %>% 
      distinct()
  } else {
    categorical <- mod_metadata %>% 
      select(STN_NISKIN, OilRig_Count, DAY_NIGHT) %>% distinct()}
  #
  # Isolate numerical data as matrix
  metadata_for_rda_0 <- mod_metadata %>% 
    # Ratios for consideration
    mutate(SIL_NO3 = (SIL/NO3)) %>% 
    mutate(T_S = (Temperature/Salinity)) %>%
    mutate(DIC_TA = (DIC/TA)) %>% 
    # mutate(across(where(is.numeric), ~na_if(abs(.), Inf))) %>% 
    # select(STN_NISKIN, Depth, DIST_OUTFLOW, DIC, NH4, NO2, NO3, 
    #        Oxygen_CTD, PO4, SIL, Salinity, TA, Temperature, pH, SIL_NO3, T_S) %>% 
    select(STN_NISKIN, Depth, DIST_OUTFLOW, NH4, NO2, NO3, 
           Oxygen_CTD, PO4, DIC_TA, Temperature, pH, SIL_NO3, T_S) %>%
    distinct() %>% 
    pivot_longer(cols = -c(STN_NISKIN), names_to = "VARIABLES", values_to = "VALUE") %>% 
    mutate(VALUE = na_if(abs(VALUE), Inf)) %>% 
    group_by(STN_NISKIN, VARIABLES) %>% 
    summarise(MEAN_VAL = mean(VALUE)) %>% 
    pivot_wider(names_from = VARIABLES, values_from = MEAN_VAL) %>% 
    arrange(STN_NISKIN) %>% 
    column_to_rownames(var = "STN_NISKIN") %>% as.matrix
  # Normalize
  standarize_metadata <- decostand(metadata_for_rda_0, MARGIN = 2, method = "standardize")
  # Combine normalized numerical data with categorical
  metadata_normalized <- categorical %>% 
    right_join(as.data.frame(standarize_metadata) %>% 
                 rownames_to_column(var = "STN_NISKIN")) %>%
    arrange(STN_NISKIN) %>% 
    column_to_rownames(var = "STN_NISKIN")
  #
  ## Prepare ASV data
  asv_wide_avg_reps_matrix <- mod_asv_long %>% 
    ungroup() %>%  
    select(FeatureID, STN_NISKIN, MEAN_REPS_seq) %>% 
    pivot_wider(names_from = STN_NISKIN, values_from = MEAN_REPS_seq, values_fn = mean, values_fill = 0) %>% 
    column_to_rownames(var = "FeatureID") %>% 
    as.matrix()
  # Transform
  asv_wide_avg_rda_0 <- data.frame(t(asv_wide_avg_reps_matrix)) %>% 
    rownames_to_column(var = "STN_NISKIN") %>% 
    arrange(STN_NISKIN) %>% 
    column_to_rownames(var = "STN_NISKIN") %>% 
    as.matrix
  #
  # CLR normalization
  clr_asv_matrix <- data.frame(compositions::clr(asv_wide_avg_rda_0))
  #
  # setdiff(rownames(metadata_normalized), (rownames(clr_asv_matrix)))
  # RDA step
  rda_output <- rda(clr_asv_matrix ~ ., data = metadata_normalized, na.action = na.exclude)
  # rda_output
  #
  # R-squared value
  adj_rda_out <- RsquareAdj(rda_output)$adj.r.squared #adj_rda_out
  # Run ANOVA
  anova_out <- anova.cca(rda_output, step = 100, by = "term")
  #
  # Get stats
  ## Percentage constrained
  total_amount <- (rda_output$tot.chi)
  constrained <- (rda_output$CCA$tot.chi)
  perc_constratained <- constrained/total_amount
  # Formula
  formula <- as.character(rda_output$call)[2]
  # Create data frame output for RDAs run
  rda_results_df <- data.frame(scores(anova_out)) %>% 
    rownames_to_column(var = "VARIABLES") %>% 
    select(VARIABLES, P_VALUE = `Pr..F.`) %>% 
    add_column(TEST = DATASET,
               PERC_CONSTRAINED = perc_constratained,
               ADJ_R_SQUARED = adj_rda_out,
               FORMULA = formula)
  #
  return(rda_results_df)
}

##
# All
##
all_samples <- table_supp %>% ungroup()
all_samples_list <- as.character(unique(table_supp$STN_NISKIN))
asvs_all_samples <- asv_long_avg_wtax_clean %>%
  filter(STN_NISKIN %in% all_samples_list) %>% ungroup()

all_rda <- run_rda(all_samples, asvs_all_samples,
                   "All samples", features = "include")
# all_rda
print("all_rda")
##
# Above/Below DCM
##
surf_table <- table_supp %>% 
  filter(grepl("Above",DCM_THRESHOLD)) %>% 
  ungroup()
surf_subset <- as.character(unique(surf_table$STN_NISKIN))
surf_asv <- asv_long_avg_wtax_clean %>%
  filter(STN_NISKIN %in% surf_subset) %>% ungroup()
above_dcm_rda <- run_rda(surf_table, surf_asv, "Above DCM", features = "exclude")
print("above_dcm_rda")

below_table <- table_supp %>% 
  filter(!(grepl("Above",DCM_THRESHOLD))) %>% 
  ungroup()
below_subset <- as.character(unique(below_table$STN_NISKIN))
below_asv <- asv_long_avg_wtax_clean %>%
  filter(STN_NISKIN %in% below_subset) %>% ungroup()
below_dcm_rda <- run_rda(below_table, below_asv, "Below DCM", features = "exclude")
# below_dcm_rda
print("below_dcm_rda")
##
# Above/Below O2 minimum
below_02_table <- table_supp %>% 
  filter(!(grepl("Above",DCM_THRESHOLD))) %>% 
  ungroup()
below_02_subset <- as.character(unique(below_02_table$STN_NISKIN))
below_02_asv <- asv_long_avg_wtax_clean %>% 
  filter(STN_NISKIN %in% below_02_subset) %>% ungroup()
below_02_rda <- run_rda(below_02_table, below_02_asv, "Below O2", features = "exclude")
# below_02_rda
print("below_02_rda")

above_02_table <- table_supp %>% 
  filter((grepl("Above",DCM_THRESHOLD))) %>% 
  ungroup()
above_02_subset <- as.character(unique(above_02_table$STN_NISKIN))
above_02_asv <- asv_long_avg_wtax_clean %>%
  filter(STN_NISKIN %in% above_02_subset) %>% ungroup()
above_02_rda <- run_rda(above_02_table, above_02_asv, "Above O2", features = "exclude")
# above_02_rda
print("above_02_rda")
##
# Above/Below MLD
below_DCM_table <- table_supp %>%
  filter(!(grepl("Above", DCM_THRESHOLD))) %>%
  ungroup()
below_DCM_subset <- as.character(unique(below_DCM_table$STN_NISKIN))
# below_MLD_asv <- asv_long_avg_wtax_clean %>% 
#   filter(STN_NISKIN %in% below_MLD_subset) %>% ungroup()
# below_MLD_rda <- run_rda(below_MLD_table, below_MLD_asv, "Below MLD", features = "exclude")
# # below_MLD_rda
# print("below_MLD_rda")
# 
above_DCM_table <- table_supp %>%
  filter((grepl("Above", DCM_THRESHOLD))) %>%
  ungroup()
above_DCM_subset <- as.character(unique(above_DCM_table$STN_NISKIN))
# above_MLD_asv <- asv_long_avg_wtax_clean %>% 
#   filter(STN_NISKIN %in% above_MLD_subset) %>% ungroup()
# above_MLD_rda <- run_rda(above_MLD_table, above_MLD_asv, "Above MLD", features = "exclude")
# print("above_MLD_rda")

##
# Dinoflagellates
dino_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "TSAR-Dinoflagellata") %>% 
  filter(Class != "Syndiniales") %>%
  ungroup()

dino_rda <- run_rda(table_supp, dino_asv, "All dinoflagellates (no Syn)", features = "include")
# dino_rda

aboveDCM_dino_asv <- dino_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
dino_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_dino_asv, "Dinoflagellates (no Syn) above DCM", features = "exclude")

belowDCM_dino_asv <- dino_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
dino_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_dino_asv, "Dinoflagellates (no Syn) below DCM", features = "exclude")

print("Dinos done.")

## 
# Syndiniales
syn_asv <- asv_long_avg_wtax_clean %>%
  filter(Class == "Syndiniales") %>% 
  ungroup()

syndiniales_rda <- run_rda(table_supp, syn_asv, "Syndiniales", features = "include")

aboveDCM_syn_asv <- syn_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
syn_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_syn_asv, "Syndiniales above DCM", features = "exclude")

belowDCM_syn_asv <- syn_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
syn_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_syn_asv, "Syndiniales below DCM", features = "exclude")

print("Syndiniales done.")

##
# Radiolaria
rad_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "TSAR-Radiolaria") %>% 
  ungroup()

rad_rda <- run_rda(table_supp, rad_asv, "Radiolaria", features = "include")

aboveDCM_rad_asv <- rad_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
rad_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_rad_asv, "Radiolaria above DCM", features = "exclude")

belowDCM_rad_asv <- rad_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
rad_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_rad_asv, "Radiolaria below DCM", features = "exclude")

print("Radiolaria done.")

##
# Ciliate
cil_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "TSAR-Ciliophora") %>% 
ungroup()

cil_rda <- run_rda(table_supp, cil_asv, "Ciliate", features = "include")

aboveDCM_cil_asv <- cil_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
cil_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_cil_asv, "Ciliate above DCM", features = "exclude")

belowDCM_cil_asv <- cil_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
cil_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_cil_asv, "Ciliate below DCM", features = "exclude")
print("Ciliates done.")

# Gyrista:
gyr_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "TSAR-Gyrista") %>%
  ungroup()

gyr_rda <- run_rda(table_supp, gyr_asv, "Gyrista", features = "include")

aboveDCM_gyr_asv <- gyr_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
gyr_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_gyr_asv, "Gyrista above DCM", features = "exclude")

belowDCM_gyr_asv <- gyr_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
gyr_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_gyr_asv, "Gyrista below DCM", features = "exclude")
# Class == "Mediophyceae"
# Class == "Bacillariophyceae"
print("Gyrista done.")
#
# Bigyra:
bigy_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "TSAR-Bigyra") %>% 
  ungroup()
big_rda <- run_rda(table_supp, bigy_asv, "Bigyra", features = "include")

aboveDCM_big_asv <- bigy_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
big_rda_aboveDCM <- run_rda(above_DCM_table, aboveDCM_big_asv, "Bigyra above DCM", features = "exclude")

belowDCM_big_asv <- bigy_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
big_rda_belowDCM <- run_rda(below_DCM_table, belowDCM_big_asv, "Bigyra below DCM", features = "exclude")
print("Bigyra done.")
# 
# Haptista
hap_asv <- asv_long_avg_wtax_clean %>%
  filter(Supergroup_simplified == "Haptista") %>%
  ungroup()
hap_subset <- as.character(unique(hap_asv$STN_NISKIN))
hap_table <- table_supp %>% 
  filter(STN_NISKIN %in% hap_subset)
hap_rda <- run_rda(hap_table, hap_asv, "Haptista", features = "include")

aboveDCM_hap_asv <- hap_asv %>% 
  filter(STN_NISKIN %in% above_DCM_subset) %>% ungroup()
hap_subset <- as.character(unique(aboveDCM_hap_asv$STN_NISKIN))
hap_table <- table_supp %>% 
  filter(STN_NISKIN %in% hap_subset)
hap_rda_aboveDCM <- run_rda(hap_table, aboveDCM_hap_asv, "Haptista above DCM", features = "exclude")

belowDCM_hap_asv <- hap_asv %>% 
  filter(STN_NISKIN %in% below_DCM_subset) %>% ungroup()
hap_subset <- as.character(unique(belowDCM_hap_asv$STN_NISKIN))
hap_table <- table_supp %>% 
  filter(STN_NISKIN %in% hap_subset)
hap_rda_belowDCM <- run_rda(hap_table, belowDCM_hap_asv, "Haptista below DCM", features = "exclude")

print("Haptista done.")

# Compile all

df_rda_results <- rbind(all_rda, above_dcm_rda, below_dcm_rda, 
                        above_02_rda, below_02_rda, 
                        # above_MLD_rda, below_MLD_rda, 
                        dino_rda, dino_rda_aboveDCM, dino_rda_belowDCM,
                        syndiniales_rda, syn_rda_aboveDCM, syn_rda_belowDCM, 
                        rad_rda, rad_rda_aboveDCM, rad_rda_belowDCM, 
                        cil_rda, cil_rda_aboveDCM, cil_rda_belowDCM, 
                        gyr_rda, gyr_rda_aboveDCM, gyr_rda_belowDCM,
                        big_rda, big_rda_aboveDCM, big_rda_belowDCM,
                        hap_rda, hap_rda_aboveDCM, hap_rda_belowDCM) %>% 
  pivot_wider(names_from = VARIABLES, values_from = P_VALUE, values_fill = NA)

save(df_rda_results, file = "/scratch/group/hu-lab/GoM-amplicon-analysis/GoM-r-outputs/df_rda_output_06162026.RData")
