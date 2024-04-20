

# Process PLFA-SIP data
# ---------------------

# Script initiation date: 2024-03-09

# In this script, PLFA-SIP data are compiled and processed in order to obtain
# absolute or relative PLFA concentrations and their 13C-isotopic compositions.
# Further on, the PLFAs are assigned to microbial groups and several indicators
# are calculated.
# The ultimate goal is to assess how sulfamethoxazole contamination affects
# the extent to which different microbial groups incorporate native soil
# organic matter (no 13C enrichment) as well as labile glucose (13C-enriched).
# Overall, 13C incorporation is high but still within the linear range
# - Watzinger et al. (2015).




# Define required packages ----
stopifnot(require("tidyverse"),
          require("openxlsx"),
          require("assertthat"),
          require("broom"),
          require("ggtext"),
          require("ggdist"),
          require("geomtextpath"))



# Import data ----

## Sample treatments ----
treatments <-
  openxlsx::read.xlsx(paste0("data/raw_data/",
                             "List of sample numbers and treatments of SMX ",
                             "experiment.xlsx"),
                      sheet = 1) %>%
  rename(sample = "Sample.label") %>%
  rename(soil = Soil) %>%
  rename(glucose_g_c_per_kg = `Glucose.dose.(g.C.kg-1.soil)`) %>%
  rename(smx_mg_per_kg = `SMX.dose.(mg.SMX.kg-1.soil)`) %>%
  select(sample, soil, glucose_g_c_per_kg, smx_mg_per_kg) %>%
  mutate(
    sample_collection = ifelse(nchar(sample) == 3,
                               "initial",
                               "after 29 days incubation"),
    glucose_g_c_per_kg = ifelse(nchar(sample) == 3,
                                NA_real_,
                                glucose_g_c_per_kg)) %>%
  # Add a column with the SMX dose in log scale
  mutate(
    log_smx = ifelse(
      !is.na(smx_mg_per_kg) & smx_mg_per_kg == 0,
      0,
      log10(10000 * smx_mg_per_kg))) %>%
  # Add a column with only the combination of soil type x glucose/SMX dose
  mutate(treatment_per_soil = str_extract(sample, "^[^-]+")) %>%
  # Add a column with only the glucose/SMX dose
  mutate(
    treatment = ifelse(
      nchar(treatment_per_soil) >= 2,
      substr(treatment_per_soil, 2, nchar(treatment_per_soil)),
      NA_character_))


## PLFA samples ----

plfa_batch1 <-
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      startRow = 2,
                      sheet = "raw data")

plfa_batch2 <-
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch2_",
                             "relative.xlsx"),
                      startRow = 2,
                      sheet = "raw data")

plfa_batch3 <-
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_quantification_batch3_",
                             "absolute.xlsx"),
                      startRow = 2,
                      sheet = "raw data")

# Rename duplicated column names
colnames(plfa_batch1) <- make.unique(colnames(plfa_batch1))
colnames(plfa_batch2) <- make.unique(colnames(plfa_batch2))
colnames(plfa_batch3) <- make.unique(colnames(plfa_batch3))


## PLFA blanks ----

blanks <- bind_rows(
  # Batch 1
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      startRow = 2,
                      cols = 1:9,
                      sheet = "Blanks") %>%
    mutate(batch = 1),
  # Batch 2
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch2_",
                             "relative.xlsx"),
                      startRow = 2,
                      cols = 1:9,
                      sheet = "Blanks") %>%
    mutate(batch = 2),
  # Batch 3
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_quantification_batch3_",
                             "absolute.xlsx"),
                      startRow = 2,
                      cols = 1:9,
                      sheet = "Blanks") %>%
    mutate(batch = 3)) %>%
  rename(sample = Sample) %>%
  rename(plfa = PLFA) %>%
  rename(peak_nr = "Peak.nr") %>%
  rename(start_s = "Start.(S)") %>%
  rename(rt_s = "RT(s)") %>%
  rename(width = Width) %>%
  rename(ampl_44 = "Ampl.44") %>%
  rename(area = "Area.All") %>%
  # Delta notation relative to Vienna Pee Dee Belemnite (VPDB)
  rename(d13c_permille = "d13C/d12C.permille") %>%
  # Use average
  filter(sample == "Blank Average") %>%
  filter(!is.na(area)) %>%
  # Harmonise PLFAs
  mutate(plfa = str_replace_all(plfa, "\\s+", ""),
         plfa = case_when(
           plfa == "i17_1w8+10Me16:0" ~ "i17:1w8+10Me16:0",
           plfa == "16:2w6,9c" ~ "18:2w6c",
           .default = plfa)) %>%
  # Note that internal standards 13:0 (to all batches) and 19:0 (to batch 3)
  # have been added to the blanks too, so we should actually ignore those.
  filter(!plfa %in% c("13:0", "19:0")) %>%
  # Create a unique key per PLFA x batch combination to link with the data
  mutate(unique_plfa_batch = paste0(plfa, "_", batch))

# Note that most of the peak areas of blanks are the same in batch 1 and
# batch 2. Maybe, these data were not fully updated?




## C atoms before and after methylation ----

c_methyl_corr <- bind_rows(
  # Batch 1
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      startRow = 2,
                      cols = 1:6,
                      sheet = "MeOH corr") %>%
    mutate(batch = 1),
  # Batch 2
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch2_",
                             "relative.xlsx"),
                      startRow = 2,
                      cols = 1:6,
                      sheet = "MeOH corr") %>%
    mutate(batch = 2),
  # Batch 3
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_quantification_batch3_",
                             "absolute.xlsx"),
                      startRow = 2,
                      cols = 1:6,
                      sheet = "MeOH corr") %>%
    mutate(batch = 3)) %>%
  select(-Sample) %>%
  rename(plfa = PLFA) %>%
  # Number of methyl groups (with one C atom per methyl group)
  # added per fatty acid
  rename(c_atoms_methyl = "C.fraction") %>%
  # Measured δ13C of methanol for methyl group
  rename(d13c_methanol = "d13.MeOH") %>%
  # Number of C atoms before methylation (i.e. from fatty acid)
  rename(c_atoms_fa = "f.FA") %>%
  # Number of C atoms after methylation (i.e. from fatty acid methyl ester)
  rename(c_atoms_fame = "f.FAME") %>%
  # Harmonise PLFAs
  mutate(
    # Remove spaces
    plfa = str_replace_all(plfa, "\\s+", ""),
    # Change "W" to "a" (typo for sample C2-2 in batch 2)
    plfa = str_replace(plfa, "W", "a"),
    plfa = case_when(
      plfa == "16:1w7c+6c" ~ "16:1w7c+16:1w6c",
      plfa == "16:2w6,9c" ~ "18:2w6c",
      plfa == "i17_1w8+10Me16:0" ~ "i17:1w8+10Me16:0",
      plfa == "18:1w7c/9t+8c" ~ "18:1w7c+18:1w9t+18:1w8c",
      plfa == "cy19:0+cy19:0" ~ "cy19:0",
      .default = plfa)) %>%
  group_by(plfa, c_atoms_fa) %>%
  reframe(c_atoms_methyl = unique(c_atoms_methyl),
          d13c_methanol = unique(d13c_methanol),
          c_atoms_fa = unique(c_atoms_fa),
          c_atoms_fame = unique(c_atoms_fame)) %>%
  ungroup() %>%
  arrange(c_atoms_fa)

# Should have one unique value per plfa
assertthat::assert_that(
  length(unique(c_methyl_corr$plfa)) == nrow(c_methyl_corr))

# The sum of columns c_atoms_fa and c_atoms_methyl should equal c_atoms_fame
assertthat::assert_that(
  all(c_methyl_corr %>%
        rowwise() %>%
        mutate(equal =
                 ((c_atoms_fa + c_atoms_methyl) == c_atoms_fame)) %>%
        pull(equal)))




## Sample masses and extracted volume ratios ----

sample_mass_vol <- bind_rows(
  # Batch 1
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      startRow = 1,
                      cols = 1:5,
                      sheet = "plfa nmol") %>%
    mutate(batch = 1),
  # Batch 2
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch2_",
                             "relative.xlsx"),
                      startRow = 1,
                      cols = 1:5,
                      sheet = "plfa nmol") %>%
    mutate(batch = 2),
  # Batch 3
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_quantification_batch3_",
                             "absolute.xlsx"),
                      startRow = 1,
                      cols = 1:5,
                      sheet = "plfa nmol") %>%
    mutate(batch = 3)) %>%
  select(-PLFA) %>%
  rename(sample = Sample) %>%
  # Dry weight of initial soil sample in gram
  rename(soil_dry_weight_g = `Dry.weight(g)`) %>%
  # Volume of total liquid/soil extract (chloroform) in mL
  rename(extract_vol_total_ml = `CHCl3.Input.(mL)`) %>%
  # Volume of liquid/soil extract (chloroform) you continue to use in mL
  # (subsampled volumetrically by taking 3.6 mL from the chloroform phase)
  rename(extract_vol_used_ml = `CHCl3.output.(mL)`) %>%
  distinct(sample,
           soil_dry_weight_g,
           extract_vol_total_ml,
           extract_vol_used_ml) %>%
  arrange(sample)

# Assert that there is one record per sample
assertthat::assert_that(
  length(unique(sample_mass_vol$sample)) == nrow(sample_mass_vol))





## Molecular information ----
# (of the fatty acid methyl esters in g mol-1,
# i.e. molecular weight of fatty acid +
# 1 extra C atom (1 x 12.01) + 2 extra H atoms (2 x 1.01))

molecular_mass <- bind_rows(
  # Batch 1
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      startRow = 1,
                      cols = 1:9,
                      sheet = "plfa nmol") %>%
    mutate(batch = 1),
  # Batch 2
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch2_",
                             "relative.xlsx"),
                      startRow = 1,
                      cols = 1:9,
                      sheet = "plfa nmol") %>%
    mutate(batch = 2),
  # Batch 3
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_quantification_batch3_",
                             "absolute.xlsx"),
                      startRow = 1,
                      cols = 1:9,
                      sheet = "plfa nmol") %>%
    mutate(batch = 3)) %>%
  rename(plfa = PLFA) %>%
  rename(mol_mass = "Mw.PLFA") %>%
  select(plfa, mol_mass) %>%
  # Harmonise PLFAs
  mutate(
    # Remove spaces
    plfa = str_replace_all(plfa, "\\s+", ""),
    # Change "W" to "a" (typo for sample C2-2 in batch 2)
    plfa = str_replace(plfa, "W", "a"),
    plfa = case_when(
      plfa == "16:1w7c+6c" ~ "16:1w7c+16:1w6c",
      plfa == "16:2w6,9c" ~ "18:2w6c",
      plfa == "i17_1w8+10Me16:0" ~ "i17:1w8+10Me16:0",
      plfa == "18:1w7c/9t+8c" ~ "18:1w7c+18:1w9t+18:1w8c",
      plfa == "cy19:0+cy19:0" ~ "cy19:0",
      .default = plfa)) %>%
  distinct(plfa,
           mol_mass) %>%
  arrange(mol_mass)






# Define molar masses
molar_mass_c <- 12.011
molar_mass_h <- 1.00784
molar_mass_o <- 15.999

# Create list of harmonised PLFAs with molecular information
# based on literature research

plfa_molecular_data <- data.frame(
  # PLFAs mostly validated / corrected using Norris et al. (2023)
  plfa = c("12:0",
           "13:0",
           "i14:0",
           "14:0",
           "i15:0",
           "a15:0",
           "15:0",
           "16:1w7c",
           "16:1w6c", # In same peak like previous
           "16:1w5c",
           "cy17:0",
           "i16:0",
           "16:0",
           # This should probably be:
           # "i17:1w9c" or "i17:1w10c"
           # (Norris et al., 2023)
           "i17:1w8",
           # Strange that "i17:1w8" and "10Me16:0" are
           # lumped together in one peak because
           # different masses
           "10Me16:0", # In same peak like previous
           "17:1w8c",
           "i17:0",
           "a17:0",
           "17:0",
           # Probably linoleic acid (18:2w6c)
           "18:2w6c",
           "18:1w9c",
           "18:1w7c",
           # Not in Norris et al. (2023)!
           "18:1w9t", # In same peak like previous
           "18:1w8c", # In same peak like previous
           # Probably lumped together from two omega positions
           # (w9c, w7c or w6c?)
           "cy19:0",
           "10Me17:0",
           "18:0",
           "10Me18:0",
           # Not in Norris et al. (2023)!
           "12Me18:0", # In same peak like previous
           "19:0"),
  # Retrieved using:
  # https://pubchem.ncbi.nlm.nih.gov/
  molecular_formula_fa = c("C12H24O2", "C13H26O2", "C14H28O2", "C14H28O2",
                           "C15H30O2", "C15H30O2", "C15H30O2", "C16H30O2",
                           "C16H30O2", "C16H30O2", "C17H32O2", "C16H32O2",
                           "C16H32O2", "C17H32O2", "C17H34O2", "C17H32O2",
                           "C17H34O2", "C17H34O2", "C17H34O2", "C18H32O2",
                           "C18H34O2", "C18H34O2", "C18H34O2", "C18H34O2",
                           "C19H36O2", "C18H36O2", "C18H36O2", "C19H38O2",
                           "C19H38O2", "C19H38O2"),
  # Based on Joergensen et al. (2021) and validated using other publications
  group_tier1 = c(NA, NA, "gram-positive", NA,
                  "gram-positive", "gram-positive", NA, "gram-negative",
                  NA, "fungi", "gram-negative", "gram-positive",
                  NA, NA, "gram-positive", "gram-negative",
                  "gram-positive", "gram-positive", NA, "fungi",
                  "fungi", "gram-negative", NA, NA,
                  "gram-negative", "gram-positive", NA, "gram-positive",
                  NA, NA),
  group_tier2 = c(NA, NA, "firmicutes", NA,
                  "firmicutes", "firmicutes", NA, NA,
                  NA, "amf", NA, "firmicutes",
                  NA, NA, "actinobacteria", NA,
                  "firmicutes", "firmicutes", NA, "asco- and basidiomycota",
                  "zygomycota", NA, NA, NA,
                  NA, "actinobacteria", NA, "actinobacteria",
                  NA, NA)) %>%
  mutate(group_tier2 = coalesce(group_tier2,
                                group_tier1)) %>%
  # Extract numbers of C, H, and O atoms and calculate molecular mass
  mutate(
    num_c_fame =
      as.numeric(gsub("^.*C(\\d+).*$", "\\1", molecular_formula_fa)) + 1,
    num_h_fame =
      as.numeric(gsub("^.*H(\\d+).*$", "\\1", molecular_formula_fa)) + 2,
    num_o_fame =
      as.numeric(gsub("^.*O(\\d+).*$", "\\1", molecular_formula_fa)),
    # Molecular mass (g mol-1)
    mol_mass_fame = round(
      num_c_fame * molar_mass_c + num_h_fame * molar_mass_h +
      num_o_fame * molar_mass_o, 2))






# Compare with information in original dataframes from this study

# Create function (~ left_join)

join_df_by_plfa <- function(df,
                            df_to_match,
                            col_name_to_match) {

  df_to_match <- df_to_match %>%
    select(plfa, {{col_name_to_match}})

  col_name_new <- paste0(col_name_to_match, "_compar")

  df[[col_name_new]] <- NA

  for (i in seq_len(nrow(df))) {

    row_selected_i <- df_to_match %>%
      filter(plfa == df$plfa[i])

    if (nrow(row_selected_i) == 0) {

      row_selected_i <- df_to_match %>%
        filter(grepl(df$plfa[i], plfa))
    }

    if (nrow(row_selected_i) > 1) {
      cat(df$plfa[i])
      stop
    }

    if (nrow(row_selected_i) == 1) {
      df[[col_name_new]][i] <-
        row_selected_i[[col_name_to_match]]
    }
  }

  return(df)
}

plfa_molecular_data <- plfa_molecular_data %>%
  join_df_by_plfa(df_to_match = molecular_mass,
                  col_name_to_match = "mol_mass") %>%
  mutate(mol_mass_compar = round(mol_mass_compar, 2)) %>%
  join_df_by_plfa(df_to_match = c_methyl_corr,
                  col_name_to_match = "c_atoms_fame") %>%
  mutate(diff_molar_mass = abs(round(abs(mol_mass_compar - mol_mass_fame), 2)),
         diff_c_atoms = abs(c_atoms_fame_compar - num_c_fame))



# Molar mass:
# inconsistencies for:
# · cy17:0
# · cy19:0
# · 10Me16:0 (lumped together with "i17:1w8" for which the info is correct)

assertthat::assert_that(
  all(!is.na(plfa_molecular_data$diff_molar_mass)) &&
  all(plfa_molecular_data %>%
    filter(!plfa %in% c("cy17:0", "cy19:0", "10Me16:0")) %>%
    pull(diff_molar_mass) <= 0.01))

# C atoms: no inconsistencies

assertthat::assert_that(
  all(plfa_molecular_data$diff_c_atoms == 0))


plfa_molecular_data <- plfa_molecular_data %>%
  select(plfa, num_c_fame, molecular_formula_fa, mol_mass_fame,
         group_tier1, group_tier2)


write.table(plfa_molecular_data,
            file = "./data/additional_data/plfa_molecular_data.csv",
            row.names = FALSE,
            na = "",
            sep = ";",
            dec = ".")





## Isotopic signatures of native soil organic carbon and glucose ----

# Isotope ratio (13C / 12C) in Pee Dee Belemnite as reference
isotope_ratio_pdb <- 0.011180

# Native soil organic carbon (unlabelled)

atom_perc_native <- bind_rows(
  # Seibersdorf
  openxlsx::read.xlsx(paste0("data/raw_data/Batch-1611 AMR raw.xlsx"),
                      startRow = 2,
                      cols = 1:10,
                      sheet = 1) %>%
    rename(sample = "Sample.label") %>%
    rename(d13c_permille = "δ.13C.OC") %>%
    rename(soil = "Soil") %>%
    filter(soil %in% c("Seibersdorf")) %>%
    filter(sample %in% c("S-1", "S-2", "S-3")) %>%
    mutate(d13c_permille = as.numeric(d13c_permille)) %>%
    select(soil, d13c_permille),
  # Grabenegg
  openxlsx::read.xlsx(paste0("data/raw_data/Batch-1611 AMR raw.xlsx"),
                      startRow = 45,
                      cols = 1:8,
                      sheet = 1) %>%
    rename(sample = "Sample.label") %>%
    rename(d13c_permille = "δ.13C") %>%
    rename(soil = "Soil") %>%
    filter(soil %in% c("Grabenegg")) %>%
    filter(sample %in% c("G-1", "G-2", "G-3")) %>%
    mutate(d13c_permille = as.numeric(d13c_permille)) %>%
    select(soil, d13c_permille)) %>%
  # Convert to atom%
  mutate(
    atom_perc = (100 * isotope_ratio_pdb * (d13c_permille / 1000 + 1)) /
      (1 + isotope_ratio_pdb * (d13c_permille / 1000 + 1))) %>%
  # Take average
  group_by(soil) %>%
  reframe(atom_perc_min = min(atom_perc),
          atom_perc_max = max(atom_perc),
          atom_perc = mean(atom_perc)) %>%
  ungroup()


# Glucose (labelled)

atom_perc_gluc <-
  openxlsx::read.xlsx(paste0("data/raw_data/Batch-1611 AMR raw.xlsx"),
                      startRow = 84,
                      cols = 1:3,
                      sheet = 1) %>%
  rename(d13c_permille = "δ.13C") %>%
  filter(sample %in% c("Glucose")) %>%
  mutate(d13c_permille = as.numeric(d13c_permille)) %>%
  select(sample, d13c_permille) %>%
  # Convert to atom%
  mutate(
    atom_perc = (100 * isotope_ratio_pdb * (d13c_permille / 1000 + 1)) /
      (1 + isotope_ratio_pdb * (d13c_permille / 1000 + 1))) %>%
  # Take average
  group_by(sample) %>%
  reframe(atom_perc_min = min(atom_perc),
          atom_perc_max = max(atom_perc),
          atom_perc = mean(atom_perc)) %>%
  ungroup()





# Harmonise raw PLFA data ----

plfa_df <- bind_rows(
  # Batch 1
  plfa_batch1 %>%
    # Remove columns 3 to 9
    select(-c(3:9)) %>%
    # This column should be numeric
    mutate(`d13C/d12C.permille.1` =
             ifelse(`d13C/d12C.permille.1` == " -",
                    NA,
                    `d13C/d12C.permille.1`),
           `d13C/d12C.permille.1` = as.numeric(`d13C/d12C.permille.1`)) %>%
    # Pivot longer
    # (one record per replicate measurement)
    pivot_longer(cols = -c(Sample, PLFA),
                 names_to = c(".value", "replicate"),
                 names_pattern = "(.*)\\.(\\d)") %>%
    mutate(batch = 1),
  # Batch 2
  plfa_batch2 %>%
    # Remove columns 3 to 9
    select(-c(3:9)) %>%
    # This column should be numeric
    mutate(`d13C/d12C.permille.1` =
             ifelse(`d13C/d12C.permille.1` == " -",
                    NA,
                    `d13C/d12C.permille.1`),
           `d13C/d12C.permille.1` = as.numeric(`d13C/d12C.permille.1`)) %>%
    # Pivot longer
    # (one record per replicate measurement)
    pivot_longer(cols = -c(Sample, PLFA),
                 names_to = c(".value", "replicate"),
                 names_pattern = "(.*)\\.(\\d)") %>%
    mutate(batch = 2),
  # Batch 3
  plfa_batch3 %>%
    # Remove columns 3 to 9
    select(-c(3:9)) %>%
    # This column should be numeric
    mutate(`d13C/d12C.permille.1` =
             ifelse(`d13C/d12C.permille.1` == " -",
                    NA,
                    `d13C/d12C.permille.1`),
           `d13C/d12C.permille.1` = as.numeric(`d13C/d12C.permille.1`)) %>%
    # Pivot longer
    # (one record per replicate measurement)
    pivot_longer(cols = -c(Sample, PLFA),
                 names_to = c(".value", "replicate"),
                 names_pattern = "(.*)\\.(\\d)") %>%
    mutate(batch = 3)) %>%
  # Sample code
  rename(sample = Sample) %>%
  # (Lumped) phospholipid fatty acid
  rename(plfa = PLFA) %>%
  rename(peak_nr = "Peak.nr") %>%
  # Start retention time (seconds)
  rename(start_s = "Start.(S)") %>%
  rename(rt_s = "RT(s)") %>%
  # Peak width
  rename(width = Width) %>%
  # Peak amplitude at mass-to-charge ratio (m/z) = 44
  # (i.e. m/z for 12C-CO2 produced in the GC-IRMS)
  rename(ampl_44 = "Ampl.44") %>%
  # Peak area (directly linked with the mass-based PLFA concentration)
  rename(area = "Area.All") %>%
  # Delta notation relative to Vienna Pee Dee Belemnite (VPDB)
  # of the measured fatty acid methyl ester (FAME)
  rename(d13c_permille = "d13C/d12C.permille") %>%
  # 12:0 is less reliable (it often disappears) and is less interesting (?)
  filter(plfa != "12") %>%
  # Remove data for 19:0 internal standard in batch 1 and 2 since not reliable
  # (added in methylised form)
  filter(!(plfa == "19:0" & batch %in% c(1, 2))) %>%
  # Harmonise PLFAs
  mutate(
    # Remove spaces
    plfa = str_replace_all(plfa, "\\s+", ""),
    # Change "17" to "17:0"
    plfa = ifelse(plfa == "17",
                  "17:0",
                  plfa),
    # Change "W" to "a" (typo for sample C2-2 in batch 2)
    plfa = str_replace(plfa, "W", "a"),
    plfa = case_when(
      plfa == "16:1w7c+6c" ~ "16:1w7c+16:1w6c",
      plfa == "i17_1w8+10Me16:0" ~ "i17:1w8+10Me16:0",
      plfa == "16:2w6,9c" ~ "18:2w6c",
      plfa == "18:1w7c/9t+8c" ~ "18:1w7c+18:1w9t+18:1w8c",
      plfa == "cy19:0+cy19:0" ~ "cy19:0",
      .default = plfa)) %>%
  # Add treatment information
  left_join(treatments %>%
              select(sample, soil, glucose_g_c_per_kg, smx_mg_per_kg,
                     log_smx, treatment_per_soil, treatment),
            by = "sample") %>%
  # Samples "C1-1" and "C2-2" are unknown
  filter(!is.na(soil)) %>%
  # Create a unique key per PLFA x batch combination to link with the blanks
  mutate(unique_plfa_batch = paste0(plfa, "_", batch)) %>%
  # Create a unique_key per PLFA x sample combination
  mutate(unique_plfa_sample = paste0(plfa, "_", sample)) %>%
  # Create copies with the original data
  mutate(area_orig = area,
         d13c_permille_orig = d13c_permille) %>%
  # Remove δ13C signatures for internal standard 19:0 to avoid confusion
  # (unneeded since added in non-methylised form)
  mutate(
    d13c_permille = ifelse(
      plfa == "19:0",
      NA_real_,
      d13c_permille)) %>%
  # Remove peak areas for internal standard 13:0 to avoid confusion
  # (unneeded since added after the extraction steps)
  mutate(
    area = ifelse(
      plfa == "13:0",
      NA_real_,
      area))




# Process PLFA data ----

## Validate correct peak selection ----

# Check if the correct peaks were selected. Usually, the start retention
# time between different runs should be < 1-2 seconds apart.

# Create a list with the plausible start retention times per PLFA.
# Median retention times per batch per PLFA seem quite stable (< 1-2 sec),
# except for "cy19:0" and for "i17:1w8+10Me16:0" (batch 3) (~ 5 sec)

rt_plausible <-
  plfa_df %>%
  group_by(plfa, batch) %>%
  reframe(median_rt = median(start_s, na.rm = TRUE),
          min_rt = min(start_s, na.rm = TRUE),
          max_rt = max(start_s, na.rm = TRUE),
          diff_rt = ifelse(!is.na(max_rt) & !is.na(min_rt),
                           max_rt - min_rt,
                           NA_real_)) %>%
  ungroup() %>%
  group_by(plfa) %>%
  reframe(stdev = round(sd(median_rt, na.rm = TRUE), 1),
          min_rt = min(min_rt, na.rm = TRUE),
          max_rt = max(max_rt, na.rm = TRUE),
          median_rt = median(median_rt, na.rm = TRUE)) %>%
  ungroup() %>%
  # < 2 seconds difference is okay. Some peaks are more difficult to demarcate.
  # For those, differences < 10, < 8 or < 5 seconds can be tolerated (based on
  # checking of the data)
  mutate(
    min_rt_plaus = case_when(
      plfa %in% c("i17:1w8+10Me16:0") ~ median_rt - 10,
      plfa %in% c("cy19:0") ~ median_rt - 8,
      plfa %in% c("10Me18:0+12Me18:0") ~ median_rt - 5,
      TRUE ~ median_rt - 2),
    max_rt_plaus = case_when(
      plfa %in% c("i17:1w8+10Me16:0") ~ median_rt + 10,
      plfa %in% c("cy19:0") ~ median_rt + 8,
      plfa %in% c("10Me18:0+12Me18:0") ~ median_rt + 5,
      TRUE ~ median_rt + 2)) %>%
  arrange(plfa)

# Add plausible retention times to plfa_df

plfa_df <- plfa_df %>%
  left_join(rt_plausible %>%
              select(plfa, min_rt_plaus, max_rt_plaus, median_rt),
            by = "plfa") %>%
  mutate(diff_rt_implausible =
           ifelse(round(start_s, 1) < round(min_rt_plaus, 1) |
                    round(start_s, 1) > round(max_rt_plaus, 1),
                  round(start_s - median_rt, 1),
                  NA_real_)) %>%
  relocate(diff_rt_implausible, median_rt, .after = plfa) %>%
  # Remove data of records that probably represent the wrong peaks
  mutate(
    area = ifelse(
      is.na(diff_rt_implausible),
      area,
      NA_real_),
    d13c_permille = ifelse(
      is.na(diff_rt_implausible),
      d13c_permille,
      NA_real_))




## Correction δ13C signatures using 13:0 internal standard ----

# The 13:0 internal standard is added in order to validate/correct the
# δ13C signatures of the PLFAs of the samples. The δ13C of the 13:0 used
# in this analysis is -31.5 permille.

d13c_13_0_theor <- -31.5 # permille

# It is added in methylised form (i.e. as fatty acid
# methyl ester) when the other PLFAs are being methylised. Else, the δ13C of
# this 13:0 would be affected by the methylation step.

# However, the δ13C data for 13:0 are assumedly less reliable because of the
# high amplitudes for this fatty acid:

assertthat::assert_that(
  plfa_df %>%
    filter(plfa == "13:0") %>%
    pull(ampl_44) %>%
    min >
    plfa_df %>%
    filter(plfa != "13:0") %>%
    filter(!is.na(ampl_44)) %>%
    pull(ampl_44) %>%
    quantile(0.8))

# This is assumedly caused by the following fact:
# PLFA contents in the soil samples of this experiment were overall low.
# Therefore, the PLFA samples had to be "evaporated" prior to the GC-IRMS
# analysis (in order to concentrate them).
# However, 13:0 tends to evaporate more quickly, leading to overestimated
# GC-IRMS peak areas.

# If the amplitude is high, the δ13C gets shifted. This occurred to a higher
# extent for 13:0 than for the other PLFAs), because of which correction
# becomes less reliable.

# However, we can compare the range after correction with the original range,
# just to see.

# Therefore, create a summary dataframe with δ13C signatures of the 13:0
# internal standard.

internal_std_d13c <- plfa_df %>%
  filter(plfa == "13:0") %>%
  group_by(sample) %>%
  reframe(d13c_permille_internal_std =
            mean(d13c_permille_orig, na.rm = TRUE)) %>%
  ungroup()

assertthat::assert_that(
  all(!is.na(internal_std_d13c$d13c_permille_internal_std)))

summary(internal_std_d13c$d13c_permille_internal_std)

plfa_df <- plfa_df %>%
  left_join(internal_std_d13c,
            by = "sample") %>%
  # Formula from SOP Watzinger and Hood-Nowotny (2019)
  # Basically:
  # R(true sample) / R(true std) =
  #    R(measured sample) / R(measured std)
  # with R = (13C / 12C)
  mutate(
    d13c_permille_corr =
      (d13c_permille + 1000) * (d13c_13_0_theor + 1000) /
      (d13c_permille_internal_std + 1000) - 1000)

# Compare the dataset summaries

plfa_enriched <- plfa_df %>%
  filter(glucose_g_c_per_kg > 0) %>%
  filter(!is.na(d13c_permille))

plfa_nat_ab <- plfa_df %>%
  filter(glucose_g_c_per_kg == 0) %>%
  filter(!is.na(d13c_permille))

summary(plfa_enriched$d13c_permille)
summary(plfa_enriched$d13c_permille_corr)

summary(plfa_nat_ab$d13c_permille)
summary(plfa_nat_ab$d13c_permille_corr)



# Check the density distributions of uncorrected versus corrected data

density_d13c_enriched <- data.frame(
  x = ggdist::density_bounded(plfa_enriched$d13c_permille)$x,
  y = ggdist::density_bounded(plfa_enriched$d13c_permille)$y)

density_d13c_corr_enriched <- data.frame(
  x = ggdist::density_bounded(plfa_enriched$d13c_permille_corr)$x,
  y = ggdist::density_bounded(plfa_enriched$d13c_permille_corr)$y)


density_d13c_nat_ab <- data.frame(
  x = ggdist::density_bounded(plfa_nat_ab$d13c_permille)$x,
  y = ggdist::density_bounded(plfa_nat_ab$d13c_permille)$y)

density_d13c_corr_nat_ab <- data.frame(
  x = ggdist::density_bounded(plfa_nat_ab$d13c_permille_corr)$x,
  y = ggdist::density_bounded(plfa_nat_ab$d13c_permille_corr)$y)


ggplot() +
  theme_minimal() +
  theme(axis.title = element_markdown(size = 10,
                                      hjust = 0.95,
                                      margin = margin(t = 10)),
        axis.text.y = element_blank()) +
  geomtextpath::geom_textline(data = density_d13c_enriched,
                              aes(x = x, y = y),
                              label = "Uncorrected (enriched)",
                              size = 3,
                              linewidth = 1,
                              vjust = 0,
                              hjust = 0.45,
                              col = "#d44102",
                              linecolor = "#d44102") +
  geomtextpath::geom_textline(data = density_d13c_corr_enriched,
                              aes(x = x, y = y),
                              label = "Corrected (enriched)",
                              size = 3,
                              linewidth = 1,
                              linetype = 2,
                              vjust = 0,
                              hjust = 0.85,
                              col = "#d44102",
                              linecolor = "#d44102") +
  geomtextpath::geom_textline(data = density_d13c_nat_ab,
                              aes(x = x, y = 0.12 * y),
                              label = "Uncorrected (nat. ab.)",
                              size = 3,
                              linewidth = 1,
                              vjust = 0,
                              hjust = 0.25,
                              col = "#77203B",
                              linecolor = "#77203B") +
  geomtextpath::geom_textline(data = density_d13c_corr_nat_ab,
                              aes(x = x, y = 0.12 * y),
                              label = "Corrected (nat. ab.)",
                              size = 3,
                              linewidth = 1,
                              linetype = 2,
                              vjust = 0,
                              hjust = 0.75,
                              col = "#77203B",
                              linecolor = "#77203B") +
  labs(x = "**δ<sup>13</sup>C** (‰)",
       y = NULL)


# Conclusion: data ranges for both are okay. However, better to proceed with
# the measured data since the δ13C signatures of 13:0 may be distorted
# (considering their high amplitudes).

plfa_df <- plfa_df %>%
  rename(d13c_permille_corr_dontuse = d13c_permille_corr)




## Validate range of δ13C signatures based on plausible range ----

# Should be > -40 permille
# (ignore internal standard 19:0 for which signatures do not matter)

plfa_df <- plfa_df %>%
  # Remove data of records for which δ13C is too low
  mutate(
    area = ifelse(
      plfa == "19:0" | d13c_permille_orig > -40,
      area,
      NA_real_),
    d13c_permille = ifelse(
      plfa == "19:0" | d13c_permille_orig > -40,
      d13c_permille,
      NA_real_))




## Validate range of δ13C signatures based on variation replicates ----
# (ideally: standard deviation < 0.5 permille)

# Create summary table with sample x plfa combinations with high stdev
d13c_stdev_high <- plfa_df %>%
  group_by(sample, plfa) %>%
  reframe(stdev_d13c = ifelse(any(!is.na(d13c_permille)),
                              round(sd(d13c_permille, na.rm = TRUE), 1),
                              NA_real_)) %>%
  # The 0.5 permille limit seems very optimistic for this dataset.
  # Take a threshold stdev of 8 permille (to flag approximately 5 % of the
  # data)
  filter(stdev_d13c > 8) %>%
  # Make key to label implausible replicate combinations
  mutate(unique_plfa_sample = paste0(plfa, "_", sample))

plfa_df <- plfa_df %>%
  mutate(
    high_d13c_stdev_flag = ifelse(
      unique_plfa_sample %in% d13c_stdev_high$unique_plfa_sample,
      "high stdev sample",
      NA_character_)) %>%
  relocate(high_d13c_stdev_flag, .before = replicate) %>%
  arrange(unique_plfa_sample)

# Usually, standard deviations between treatment replicates (across batches)
# are also high for the sample x plfa combinations that are flagged for having a
# high standard deviation.
# Standard deviations are naturally higher for PLFAs with an uncertain peak
# demarcation (e.g. "10Me18:0+12Me18:0").

# Based on visual checking, these values are potentially off:
# · G6-3 "16:1w7c+16:1w6c" rep 2
# · S6-1 "18:0" rep 1

# However, we can detect and omit outliers (> 3 x stdev). Based on that rule,
# we can detect if there are any outliers. Else, better to keep all data and
# propagate uncertainty.
# Therefore, create another summary table with stats per plfa x treatment x soil

d13c_plausible <- plfa_df %>%
  group_by(treatment_per_soil, plfa) %>%
  reframe(
    mean_d13c = ifelse(
      any(!is.na(d13c_permille)),
      round(mean(d13c_permille, na.rm = TRUE), 1),
      NA_real_),
    median_d13c = ifelse(
      any(!is.na(d13c_permille)),
      round(median(d13c_permille, na.rm = TRUE), 1),
      NA_real_),
    min_d13c = ifelse(
      any(!is.na(d13c_permille)),
      round(min(d13c_permille, na.rm = TRUE), 1),
      NA_real_),
    max_d13c = ifelse(
      any(!is.na(d13c_permille)),
      round(max(d13c_permille, na.rm = TRUE), 1),
      NA_real_),
    stdev_d13c = ifelse(
      any(!is.na(d13c_permille)),
      round(sd(d13c_permille, na.rm = TRUE), 1),
      NA_real_)) %>%
  ungroup() %>%
  mutate(
    min_d13c_plaus = ifelse(
      !is.na(mean_d13c) & !is.na(stdev_d13c),
      mean_d13c - 3 * stdev_d13c,
      mean_d13c),
    max_d13c_plaus = ifelse(
      !is.na(mean_d13c) & !is.na(stdev_d13c),
      mean_d13c + 3 * stdev_d13c,
      mean_d13c)) %>%
  mutate(
    d13c_implausible = ifelse(
      (!is.na(min_d13c) & !is.na(min_d13c_plaus) &
        min_d13c < min_d13c_plaus) |
      (!is.na(max_d13c) & !is.na(max_d13c_plaus) &
        max_d13c > max_d13c_plaus),
      "outlier",
      NA_character_))

assertthat::assert_that(
  all(is.na(d13c_plausible$d13c_implausible)))

# No outliers based on the "3-sigma rule".
# To conclude: standard deviations are overall high, but there is no need
# to remove any data. Propagate uncertainty.






## Correction δ13C signatures for methylation step ----

plfa_df <- plfa_df %>%
  # Copy the column with d13c signatures into a column called
  # "d13c_permille_fame" (since these - measured - δ13C signatures
  # represent those after the methylation step (fatty acid methyl esters),
  # while we want to know the signatures for the original fatty acids)
  mutate(d13c_permille_fame = d13c_permille) %>%
  left_join(c_methyl_corr,
            by = "plfa") %>%
  # Formula:
  # nC(fa) * δ13C(fa) + nC(methyl) * δ13C(methyl) = nC(fame) + δ13C(fame)
  # (with nC referring to the number of carbon atoms)
  mutate(
    d13c_permille =
      (c_atoms_fame * d13c_permille_fame - c_atoms_methyl * d13c_methanol) /
      c_atoms_fa)



## Convert to atom% ----

# Isotope ratio (13C / 12C) in Pee Dee Belemnite as reference
isotope_ratio_pdb <- 0.011180

# Formulas:
# atom% = Rs / (Rs + 1) * 100 %
# Rs = isotopic ratio of sample (i.e. heavy isotope / light isotope = 13C / 12C)
# δ13C / 1000 + 1 = Rs (sample) / isotope_ratio_pdb

# Atom % of PDB is 1.1056 %

plfa_df <- plfa_df %>%
  mutate(
    atom_perc = (100 * isotope_ratio_pdb * (d13c_permille / 1000 + 1)) /
      (1 + isotope_ratio_pdb * (d13c_permille / 1000 + 1)))





## Harmonise "below LOQs" ----
# Measurements with amplitude < 200 mV are considered unreliable.
# Here, we will consider this threshold as a limit of quantification.

ampl_44_loq <- 200 # mV


# Check if there is any "trend" between amplitudes and peak areas

plfa_df %>%
  filter(!plfa %in% c("13:0", "19:0")) %>%
  filter(!is.na(area) & !is.na(d13c_permille)) %>%
  ggplot(aes(x = ampl_44,
             y = area,
             col = as.factor(plfa))) +
  geom_point() +
  geom_vline(xintercept = 200, col = "red") +
  guides(col = "none")

# Based on visual checking, it seems that peak areas and amplitudes are
# linearly related, but differently per PLFA.
# Create a summary dataframe with the slopes per PLFA (intercept = 0).

plfa_loqs <- plfa_df %>%
  filter(!plfa %in% c("13:0")) %>%
  filter(!is.na(area) & (!is.na(d13c_permille) | (plfa == "19:0"))) %>%
  group_by(plfa) %>%
  do(tidy(lm(area ~ ampl_44 + 0, data = .))) %>%
  filter(term == "ampl_44") %>%
  # Column "std.error" shows the standard error associated with the
  # estimated coefficient.
  # Column "statistic" shows the t-statistic for the hypothesis test of
  # whether the coefficient is significantly different from zero.
  # The P-value is sufficient to validate the outcome
  select(plfa,
         slope = estimate,
         p_value = p.value)

# Assert that the P-value is < 0.05
assertthat::assert_that(
  all(plfa_loqs$p_value < 0.05))

# For each PLFA: add LOQs in peak area instead of in ampl_44

plfa_loqs <- plfa_loqs %>%
  mutate(area_loq = ampl_44_loq * slope)

# Add LOQs to dataset

plfa_df <- plfa_df %>%
  left_join(plfa_loqs %>%
              select(plfa, area_loq),
            by = "plfa") %>%
  mutate(
    below_loq = ifelse(
      ampl_44 < ampl_44_loq,
      TRUE,
      FALSE),
    area = ifelse(
      below_loq == TRUE & !is.na(area),
      # Use 50 % of the LOQ (in peak area units)
      0.5 * area_loq,
      area))

# Fraction of below-loq records
length(which(plfa_df$below_loq == TRUE)) / nrow(plfa_df)



## Correct for peak areas measured in blanks ----

plfa_df <- plfa_df %>%
  left_join(blanks %>%
              select(unique_plfa_batch, area) %>%
              rename(area_blank = area),
            by = "unique_plfa_batch") %>%
  mutate(
    area = ifelse(
      !is.na(area_blank) &
        !is.na(area),
      ifelse(
        (area - area_blank) < 0,
        0,
        area - area_blank),
      area))



## Split "lumped PLFAs" into multiple records ----

# Some similar PLFAs give an overlapping peak (similar retention time),
# based on which it is not possible to distinguish them.
# Because of that, the column "plfa" sometimes contains multiple PLFAs
# separated by "+".
# Overall, based on tests, it is plausible to assume that we can just
# distribute the peak areas evenly.
# However, we will here propagate some uncertainty in this assumption,
# by also using an arbitrary minimum of 0.5 and 1.5 x this evenly
# distributed peak area.

plfa_split <- plfa_df %>%
  # Copy the column "plfa" into column "plfa_lumped".
  # Column "plfa_lumped" will contain the original (lumped) PLFAs
  # Column "plfa" will contain the split PLFAs (after the for loop below)
  mutate(plfa_lumped = plfa) %>%
  # Add minimum and maximum of uncertainty range
  mutate(area_min = area,
         area_max = area)

extra_rows <- NULL

for (i in seq_len(nrow(plfa_split))) {

  plfa_split_i <- unlist(strsplit(plfa_split$plfa[i], split = "[+]"))

  area_i <- plfa_split$area[i]

  extra_row_i <- NULL
  extra_row_i_2 <- NULL
  extra_row_i_3 <- NULL

  if (length(plfa_split_i) == 2) {

    # Second option

    extra_row_i <- plfa_split[i, ]

    extra_row_i$plfa[1] <- plfa_split_i[2]
    extra_row_i$area[1] <- 0.5 * area_i
    extra_row_i$area_min[1] <- 0.5 * 0.5 * area_i
    extra_row_i$area_max[1] <- 1.5 * 0.5 * area_i


    # First option

    plfa_split$plfa[i] <- plfa_split_i[1]
    plfa_split$area[i] <- 0.5 * area_i
    plfa_split$area_min[i] <- 0.5 * 0.5 * area_i
    plfa_split$area_max[i] <- 1.5 * 0.5 * area_i

  }

  if (length(plfa_split_i) == 3) {

    # Second option

    extra_row_i_2 <- plfa_split[i, ]

    extra_row_i_2$plfa[1] <- plfa_split_i[2]
    extra_row_i_2$area[1] <- 1/3 * area_i
    extra_row_i_2$area_min[1] <- 0.5 * 1/3 * area_i
    extra_row_i_2$area_max[1] <- 1.5 * 1/3 * area_i

    # Third option

    extra_row_i_3 <- plfa_split[i, ]

    extra_row_i_3$plfa[1] <- plfa_split_i[3]
    extra_row_i_3$area[1] <- 1/3 * area_i
    extra_row_i_3$area_min[1] <- 0.5 * 1/3 * area_i
    extra_row_i_3$area_max[1] <- 1.5 * 1/3 * area_i

    extra_row_i <- bind_rows(extra_row_i_2,
                             extra_row_i_3)

    # First option

    plfa_split$plfa[i] <- plfa_split_i[1]
    plfa_split$area[i] <- 1/3 * area_i
    plfa_split$area_min[i] <- 0.5 * 1/3 * area_i
    plfa_split$area_max[i] <- 1.5 * 1/3 * area_i

  }

  extra_rows <- bind_rows(extra_rows,
                          extra_row_i)

} # End of for loop over rows

plfa_split <- bind_rows(plfa_split,
                        extra_rows) %>%
  arrange(treatment, soil, batch, plfa) %>%
  # Take LOQs into account
  mutate(
    below_loq = ifelse(
      area < area_loq,
      TRUE,
      FALSE),
    area_min = ifelse(
      below_loq == TRUE & !is.na(area),
      0,
      area_min),
    area_max = ifelse(
      below_loq == TRUE & !is.na(area) & area_max < area_loq,
      area_loq,
      area_max),
    area = ifelse(
      below_loq == TRUE & !is.na(area),
      # Use 50 % of the LOQ (in peak area units)
      0.5 * area_loq,
      area))




## Convert peak areas to mole-based concentration units ----

# The 19:0 internal standard is the quantifying standard. It is added before
# the extraction and analysis in order to account for losses during the
# different extraction steps.
# It has to be added in the fatty acid form (not fatty acid methyl ester form).
# For this dataset, 19:0 was inadvertantly added in the methylised form for
# batch 1 and 2, which makes the 19:0 data useless for these two batches.

# Therefore, only relative comparisons of the quantities are possible in
# batch 1 and 2 of this dataset. As such, we will use:
# · relative units (mole-based; n/n) across the three batches
# · absolute units (nmol g-1 dry soil) in batch 3, to also get a rough idea of
#   the absolute concentrations (based on one replicate)

# The fact that we mostly have to resort to relative concentrations is not a
# problem whatsoever, since we are mainly interested in the incorporation of
# 13C-enriched glucose in different fatty acids across treatments.

# The following quantity of 19:0 was initially added (batch 3):

microgram_19_0_per_sample_added <- 22.4

# Create a dataframe with the peak area of 19:0 per sample (batch 3)

internal_std_area <- plfa_split %>%
  filter(batch == 3) %>%
  filter(plfa == "19:0") %>%
  group_by(sample) %>%
  reframe(area = mean(area, na.rm = TRUE))


# Convert units

plfa_split <- plfa_split %>%
  # Add sample masses and extracted volumes
  left_join(sample_mass_vol,
            by = "sample") %>%
  # Add molecular weights
  left_join(plfa_molecular_data %>%
              select(plfa, mol_mass_fame),
            by = "plfa") %>%
  # Add peak areas of quantitative internal standard 19:0
  left_join(internal_std_area %>%
              rename(area_19_0 = area),
            by = "sample") %>%
  # Only relevant for actual validated sample records of batch 3
  mutate(
    area_19_0 = ifelse(
      !is.na(area) &
        batch == 3 &
        !plfa %in% c("13:0", "19:0"),
      area_19_0,
      NA_real_)) %>%
  # Calculate (mole-based) area per g dry soil
  mutate(
    area_mole_based_per_g =
      1000 * (area * (extract_vol_total_ml / extract_vol_used_ml)) /
      soil_dry_weight_g / mol_mass_fame,
    area_mole_based_per_g_min =
      1000 * (area_min * (extract_vol_total_ml / extract_vol_used_ml)) /
      soil_dry_weight_g / mol_mass_fame,
    area_mole_based_per_g_max =
      1000 * (area_max * (extract_vol_total_ml / extract_vol_used_ml)) /
      soil_dry_weight_g / mol_mass_fame) %>%
  # Calculate absolute concentrations (nmol per g dry soil)
  mutate(
    conc_nmol_per_g =
      # The "(extract_vol_total_ml / extract_vol_used_ml)" ratio is actually
      # included in the added µg to peak area ratio of 19:0,
      # since the peak area of 19:0 was likewise measured in the
      # 3.6-mL subsample
      1000 * # 1000 is a unit correction factor (1000 ng per µg)
      (area * (microgram_19_0_per_sample_added / area_19_0)) /
      soil_dry_weight_g / mol_mass_fame,
    conc_nmol_per_g_min =
      # The "(extract_vol_total_ml / extract_vol_used_ml)" ratio is actually
      # included in the added µg to peak area ratio of 19:0,
      # since the peak area of 19:0 was likewise measured in the
      # 3.6-mL subsample
      1000 * (area_min * (microgram_19_0_per_sample_added / area_19_0)) /
      soil_dry_weight_g / mol_mass_fame,
    conc_nmol_per_g_max =
      # The "(extract_vol_total_ml / extract_vol_used_ml)" ratio is actually
      # included in the added µg to peak area ratio of 19:0,
      # since the peak area of 19:0 was likewise measured in the
      # 3.6-mL subsample
      1000 * (area_max * (microgram_19_0_per_sample_added / area_19_0)) /
      soil_dry_weight_g / mol_mass_fame)


# Calculate the sum of the mole-based peak areas per sample, in order to
# calculate the percentage (n/n) of each fatty acid in each sample

mean_per_sample <- plfa_split %>%
  filter(!plfa %in% c("13:0", "19:0")) %>%
  group_by(sample, treatment_per_soil, plfa) %>%
  reframe(
    area_mole_based_per_g = ifelse(
      any(!is.na(area_mole_based_per_g)),
      mean(area_mole_based_per_g, na.rm = TRUE),
      NA_real_),
    area_mole_based_per_g_min = ifelse(
      any(!is.na(area_mole_based_per_g_min)),
      min(area_mole_based_per_g_min, na.rm = TRUE),
      NA_real_),
    area_mole_based_per_g_max = ifelse(
      any(!is.na(area_mole_based_per_g_max)),
      max(area_mole_based_per_g_max, na.rm = TRUE),
      NA_real_))

# There are some NAs due to wrong peak selection for both replicates.
# The NAs are therefore do not mean that the PLFA was absent, but rather
# that that we do not know how much of the PLFA was present because it was
# not (reliably) measured.

length(which(is.na(mean_per_sample$area_mole_based_per_g)))

# In order to make a "correct" sum of the PLFAs per sample (and therefore
# comparable relative concentrations across batches), we'd better gap-fill
# the concentrations by using the range of concentrations reported in other
# repetitions of the same treatment (in other batches).

mean_per_treatment <- mean_per_sample %>%
  group_by(treatment_per_soil, plfa) %>%
  reframe(
    area_mole_based_per_g = ifelse(
      any(!is.na(area_mole_based_per_g)),
      mean(area_mole_based_per_g, na.rm = TRUE),
      NA_real_),
    area_mole_based_per_g_min = ifelse(
      any(!is.na(area_mole_based_per_g_min)),
      min(area_mole_based_per_g_min, na.rm = TRUE),
      NA_real_),
    area_mole_based_per_g_max = ifelse(
      any(!is.na(area_mole_based_per_g_max)),
      max(area_mole_based_per_g_max, na.rm = TRUE),
      NA_real_)) %>%
  rename(area_treatment = area_mole_based_per_g,
         area_treatment_min = area_mole_based_per_g_min,
         area_treatment_max = area_mole_based_per_g_max) %>%
  mutate(key = paste0(plfa, "_", treatment_per_soil))

assertthat::assert_that(
  length(which(is.na(mean_per_treatment$area_treatment))) == 0)



# Gap-fill mean_per_sample with these averages per treatment across batches

mean_per_sample <- mean_per_sample %>%
  mutate(key = paste0(plfa, "_", treatment_per_soil)) %>%
  left_join(mean_per_treatment %>%
              select(-plfa, -treatment_per_soil),
            by = "key") %>%
  select(-key) %>%
  mutate(
    area_mole_based_per_g = coalesce(area_mole_based_per_g,
                                     area_treatment),
    area_mole_based_per_g_min = coalesce(area_mole_based_per_g_min,
                                         area_treatment_min),
    area_mole_based_per_g_max = coalesce(area_mole_based_per_g_max,
                                         area_treatment_max)) %>%
  select(-area_treatment,
         -area_treatment_min,
         -area_treatment_max)

assertthat::assert_that(
  length(which(is.na(mean_per_sample$area_mole_based_per_g))) == 0)



# Calculate sum per sample

total_per_sample <- mean_per_sample %>%
  group_by(sample) %>%
  reframe(area_per_g_total = sum(area_mole_based_per_g),
          area_per_g_total_min = sum(area_mole_based_per_g_min),
          area_per_g_total_max = sum(area_mole_based_per_g_max))


# Add this to the dataset

plfa_split <- plfa_split %>%
  left_join(total_per_sample,
            by = "sample") %>%
  # Multiply with 100 to have the mol fractions as a percentage
  mutate(conc_relative = 100 * area_mole_based_per_g / area_per_g_total,
         conc_relative_min =
           100 * area_mole_based_per_g_min / area_per_g_total_max,
         conc_relative_max =
           100 * area_mole_based_per_g_max / area_per_g_total_min) %>%
  select(-area_mole_based_per_g,
         -area_mole_based_per_g_min,
         -area_mole_based_per_g_max,
         -area_per_g_total,
         -area_per_g_total_min,
         -area_per_g_total_max)




## Assign microbial groups to different fatty acids ----
# Based on Joergensen et al. (2021) (and validated using other publications)

# Firstly, check whether some PLFAs that are not mentioned in Joergensen
# et al. can be ignored.

plfa_split %>%
  filter(!plfa %in% c("13:0", "19:0")) %>%
  filter(!is.na(conc_relative)) %>%
  select(plfa, conc_relative) %>%
  ggplot(aes(x = plfa,
             y = conc_relative)) +
  geom_boxplot(alpha = .5,
               size = 1,
               outlier.size = 5) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  coord_flip()

# Three PLFAs are on average < 1 mol% of the total PLFAs.
# At this point, better not to exclude any PLFA.

plfa_split <- plfa_split %>%
  left_join(plfa_molecular_data %>%
              select(plfa, group_tier1, group_tier2),
            by = "plfa") %>%
  relocate(group_tier1, group_tier2, .after = plfa)



## Assess glucose-derived fraction ----
# (and further variables)

# Mass balance:
# fraction_plfa_gluc = (%13C_plfa_soil_gluc - %13C_plfa_native) /
#                      (%13C_plfa_gluc - %13C_plfa_native)

# With %13C_plfa_gluc: the %13C of a given PLFA derived entirely from labelled
# glucose (i.e. without native soil organic carbon). Since this is hard to
# measure and since it is justifiable to assume that the isotopic fractionation
# from glucose to plfa_gluc is similar to the isotopic fractionation from
# native soil organic carbon to plfa_native, we will simplify this formula to:

# fraction_plfa_gluc = (%13C_plfa_soil_gluc - %13C_plfa_native) /
#                      (%13C_gluc - %13C_native)


# First, derive the %13C_plfa_native, i.e. the % 13C of each PLFA in
# unamended samples ("controls").

atom_perc_plfa_native <- plfa_split %>%
  filter(glucose_g_c_per_kg == 0 &
           smx_mg_per_kg == 0) %>%
  filter(!plfa %in% c("13:0", "19:0")) %>%
  group_by(sample, soil, plfa) %>%
  reframe(
    atom_perc_min = ifelse(
      any(!is.na(atom_perc)),
      min(atom_perc, na.rm = TRUE),
      NA_real_),
    atom_perc_max = ifelse(
      any(!is.na(atom_perc)),
      max(atom_perc, na.rm = TRUE),
      NA_real_),
    atom_perc = ifelse(
      any(!is.na(atom_perc)),
      mean(atom_perc, na.rm = TRUE),
      NA_real_)) %>%
  ungroup() %>%
  group_by(soil, plfa) %>%
  reframe(atom_perc_control_min = min(atom_perc_min, na.rm = TRUE),
          atom_perc_control_max = max(atom_perc_max, na.rm = TRUE),
          atom_perc_control = mean(atom_perc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(key = paste0(soil, "_", plfa))


plfa_split <- plfa_split %>%
  # Remove internal standards because no longer needed now
  filter(!plfa %in% c("13:0", "19:0")) %>%
  # Add % 13C in PLFAs of unamended treatments
  mutate(key = paste0(soil, "_", plfa)) %>%
  left_join(atom_perc_plfa_native %>%
              select(-soil, -plfa),
            by = "key") %>%
  select(-key) %>%
  # Add % 13C in native soil organic carbon
  left_join(atom_perc_native %>%
              rename(atom_perc_soil_min = atom_perc_min) %>%
              rename(atom_perc_soil_max = atom_perc_max) %>%
              rename(atom_perc_soil = atom_perc),
            by = "soil") %>%
  mutate(
    # Fraction of PLFAs derived from 13C-enriched glucose
    frac_gluc_min = ifelse(
      glucose_g_c_per_kg > 0,
      (atom_perc - atom_perc_control_max) /
        (atom_perc_gluc$atom_perc_max - atom_perc_soil_min),
      0),
    frac_gluc_max = ifelse(
      glucose_g_c_per_kg > 0,
      (atom_perc - atom_perc_control_min) /
        (atom_perc_gluc$atom_perc_min - atom_perc_soil_max),
      0),
    frac_gluc = ifelse(
      glucose_g_c_per_kg > 0,
      (atom_perc - atom_perc_control) /
        (atom_perc_gluc$atom_perc - atom_perc_soil),
      0)) %>%
  mutate(
    # Relative concentration (mol%) of PLFAs derived from 13C-enriched glucose
    conc_rel_gluc_min = ifelse(
      glucose_g_c_per_kg > 0,
      frac_gluc_min * conc_relative_min,
      NA_real_),
    conc_rel_gluc_max = ifelse(
      glucose_g_c_per_kg > 0,
      frac_gluc_max * conc_relative_max,
      NA_real_),
    conc_rel_gluc = ifelse(
      glucose_g_c_per_kg > 0,
      frac_gluc * conc_relative,
      NA_real_),
    # Relative concentration (mol%) of PLFAs derived from native SOC
    conc_rel_native_min = ifelse(
      glucose_g_c_per_kg > 0,
      (1 - frac_gluc_max) * conc_relative_min,
      conc_relative_min),
    conc_rel_native_max = ifelse(
      glucose_g_c_per_kg > 0,
      (1 - frac_gluc_min) * conc_relative_max,
      conc_relative_max),
    conc_rel_native = ifelse(
      glucose_g_c_per_kg > 0,
      (1 - frac_gluc) * conc_relative,
      conc_relative))

# Add an index representing the percentage extra PLFAs coming from native soil
# organic carbon (as compared to unamended controls), to see to which extent
# the incorporation of native soil organic carbon-derived C in PLFAs is
# stimulated after the addition of glucose
# (i.e. % extra PLFAs coming from native soil organic carbon)

# First, create a dataframe with the concentrations of PLFAs in unamended
# soil samples:

conc_plfa_native <- plfa_split %>%
  filter(glucose_g_c_per_kg == 0 &
           smx_mg_per_kg == 0) %>%
  group_by(sample, soil, plfa) %>%
  reframe(
    conc_relative_min = ifelse(
      any(!is.na(conc_relative_min)),
      min(conc_relative_min, na.rm = TRUE),
      NA_real_),
    conc_relative_max = ifelse(
      any(!is.na(conc_relative_max)),
      max(conc_relative_max, na.rm = TRUE),
      NA_real_),
    conc_relative = ifelse(
      any(!is.na(conc_relative)),
      mean(conc_relative, na.rm = TRUE),
      NA_real_)) %>%
  ungroup() %>%
  group_by(soil, plfa) %>%
  reframe(conc_rel_control_min = min(conc_relative_min, na.rm = TRUE),
          conc_rel_control_max = max(conc_relative_max, na.rm = TRUE),
          conc_rel_control = mean(conc_relative, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(key = paste0(soil, "_", plfa))


plfa_split <- plfa_split %>%
  # Add % 13C in PLFAs of unamended treatments
  mutate(key = paste0(soil, "_", plfa)) %>%
  left_join(conc_plfa_native %>%
              select(-soil, -plfa),
            by = "key") %>%
  select(-key) %>%
  mutate(
    conc_native_perc_extra_min = ifelse(
      glucose_g_c_per_kg > 0,
      # Use "conc_rel_control" rather than "conc_rel_control_max"
      # to be consistent with conc_native_perc_extra_max (see below)
      (conc_rel_native_min / conc_rel_control - 1) * 100,
      NA_real_),
    conc_native_perc_extra_max = ifelse(
      glucose_g_c_per_kg > 0,
      # Use "conc_rel_control" rather than "conc_rel_control_min" since
      # the latter is often 0 (due to below-LOQs)
      (conc_rel_native_max / conc_rel_control - 1) * 100,
      NA_real_),
    conc_native_perc_extra = ifelse(
      glucose_g_c_per_kg > 0,
      (conc_rel_native / conc_rel_control - 1) * 100,
      NA_real_))


# Reorder columns

plfa_split <- plfa_split %>%
  select(sample,
         soil,
         glucose_g_c_per_kg,
         smx_mg_per_kg,
         log_smx,
         batch,
         treatment_per_soil,
         treatment,
         plfa,
         plfa_lumped,
         mol_mass_fame,
         contains("c_atoms"),
         group_tier1,
         group_tier2,
         replicate,
         soil_dry_weight_g,
         start_s,
         median_rt,
         min_rt_plaus,
         max_rt_plaus,
         diff_rt_implausible,
         rt_s,
         peak_nr,
         width,
         ampl_44,
         below_loq,
         contains("extract_vol"),
         contains("area"),
         contains("d13c"),
         contains("conc_nmol"),
         contains("conc_relative"),
         contains("atom_perc"),
         contains("frac"),
         contains("conc_rel_gluc"),
         contains("conc_rel_native"),
         contains("conc_native_perc_extra"))

write.table(plfa_split,
            file = "./data/final_data/plfa_per_replicate_harmonised.csv",
            row.names = FALSE,
            na = "",
            sep = ";",
            dec = ".")






## Summarise per sample x PLFA ----
# (across two GC-IRMS analysis replicates)

# Create function

summarise_per_group <- function(df,
                                variables_to_summarise) {

  # The dataframe should be grouped already
  assertthat::assert_that(is_grouped_df(df))

  # Retrieve the grouping variables
  grouping_vars <- group_vars(df)

  # Create the "key" column by pasting the values of grouping variables
  # together
  df <- df %>%
    ungroup() %>%
    rowwise() %>%
    mutate(key = paste(c_across(all_of(grouping_vars)),
                       collapse = "_")) %>%
    group_by(across(all_of(grouping_vars)), key)

  assertthat::assert_that(is.character(variables_to_summarise))


  for (var in variables_to_summarise) {

    assertthat::assert_that(
      var %in% names(df))

    var_min <- paste0(var, "_min")
    var_max <- paste0(var, "_max")
    var_stdev <- paste0(var, "_stdev")

    # Option 1: within-repetition variation is indicated in the dataframe
    # (using columns with names ending with "_min" and "_max")

    if (var_min %in% names(df) &
        var_max %in% names(df)) {

      df_summ_var <- df %>%
        reframe(
          !!var_stdev := {
            var_values <- .data[[var]]
            if (sum(!is.na(var_values)) > 1) {
              sd(var_values, na.rm = TRUE)
            } else {
              NA_real_
            }
          },
          !!var_min := ifelse(
            any(!is.na(.data[[var_min]])),
            min(.data[[var_min]], na.rm = TRUE),
            NA_real_
          ),
          !!var_max := ifelse(
            any(!is.na(.data[[var_max]])),
            max(.data[[var_max]], na.rm = TRUE),
            NA_real_
          ),
          !!var := ifelse(
            any(!is.na(.data[[var]])),
            mean(.data[[var]], na.rm = TRUE),
            NA_real_
          )) %>%
        # We are not proceeding with this standard deviation, which only reflects
        # the between-replicate variation (of the estimates), not the within-
        # replicate variation.
        # The uncertainty ranges provided for each estimate already encapsulate
        # the variation. They already gives us an understanding of the range
        # within which the true parameter value is likely to lie.
        select(-{{var_stdev}}) %>%
        relocate({{var}}, .before = {{var_min}})

    } else {

      df_summ_var <- df %>%
        reframe(
          !!var_stdev := {
            var_values <- .data[[var]]
            if (sum(!is.na(var_values)) > 1) {
              sd(var_values, na.rm = TRUE)
            } else {
              NA_real_
            }
          },
          !!var := ifelse(
            any(!is.na(.data[[var]])),
            mean(.data[[var]], na.rm = TRUE),
            NA_real_
          )) %>%
        mutate(
          !!var_min := ifelse(
            !is.na(.data[[var_stdev]]) & !is.na(.data[[var]]),
            .data[[var]] - .data[[var_stdev]],
            NA_real_),
          !!var_max := ifelse(
            !is.na(.data[[var_stdev]]) & !is.na(.data[[var]]),
            .data[[var]] + .data[[var_stdev]],
            NA_real_)) %>%
        # We can proceed with minimum and maximum values of the uncertainty
        # ranges, to be consistent across parameters
        select(-{{var_stdev}}) %>%
        relocate({{var}}, .before = {{var_min}})

    }


    if (which(var == variables_to_summarise) == 1) {

      # For the first variable
      df_summ <- df_summ_var

    } else {

      # For the next variables
      df_summ <- df_summ %>%
        left_join(df_summ_var %>%
                    select(-all_of(grouping_vars)),
                  by = "key")
    }

  } # End of "for loop" over variables

  df_summ <- df_summ %>%
    ungroup() %>%
    select(-key)

  return(df_summ)
}




# Summarise

parameters <- c("conc_relative",
                "conc_nmol_per_g",
                "atom_perc",
                "frac_gluc",
                "conc_rel_gluc",
                "conc_rel_native",
                "conc_native_perc_extra")

grouping_columns <- c("sample",
                      "soil",
                      "glucose_g_c_per_kg",
                      "smx_mg_per_kg",
                      "log_smx",
                      "batch",
                      "treatment_per_soil",
                      "treatment",
                      "plfa",
                      "group_tier1",
                      "group_tier2")


# Check which grouping columns are numeric

numeric_grouping_columns <- plfa_split %>%
  summarise(across(all_of(grouping_columns), is.numeric)) %>%
  unlist()

numeric_grouping_columns <-
  names(numeric_grouping_columns)[numeric_grouping_columns]



plfa_summ <- plfa_split %>%
  # Convert grouping columns to characters
  mutate(across(all_of(numeric_grouping_columns), as.character)) %>%
  group_by(sample,
           soil,
           glucose_g_c_per_kg,
           smx_mg_per_kg,
           log_smx,
           batch,
           treatment_per_soil,
           treatment,
           plfa,
           group_tier1,
           group_tier2) %>%
  summarise_per_group(variables_to_summarise = parameters) %>%
  mutate(across(all_of(numeric_grouping_columns), as.numeric)) %>%
  arrange(plfa,
          soil,
          treatment)



write.table(plfa_summ,
            file = "./data/final_data/plfa_per_sample_harmonised.csv",
            row.names = FALSE,
            na = "",
            sep = ";",
            dec = ".")






# Graphs ----

## Stacked columns absolute concentration (nmol g-1) (batch 3) ----
# (different groups per column, including "non-marker PLFAs")

## Ratio fungi to bacteria ----

## Ratio stress ----

## Canonical Correspondence Analysis ----

## Per group: stacked columns indicating the evolution of sources ----
# (native SOC versus glucose)
# (use facets for graph)

# conc_native_perc_extra




















