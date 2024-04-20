

# Import and process SIP-PLFA data
# ----------------------------

# Script initiation date: 2024-03-09

# In this script, SIP-PLFA data are compiled and processed in order to assess
# how sulfamethoxazole contamination affects the extent to which different
# microbial groups incorporate native soil organic matter (no 13C enrichment)
# as well as labile glucose (13C-enriched)

# Transformation from raw peak areas to PLFA concentrations (in nmol g-1 soil)
# and their δ13C-isotopic signatures happened in advance in Excel and
# followed these steps:

# - Peaks are validated as follows - recheck peaks if:
#   · the stdev between replicate measurements of δ13C are > 0.5 ‰
#   · the peak height (m/z 44) is < 200 mV (values < 100 mV should be removed)
#   · asserting linearity of IRMS (high 13C incorporation but still within
#     okay range - Watzinger et al. (2015))
#   · the start retention time between different runs are > 1 to 2 seconds
#     (since retention times within different runs should be quite stable)

# - Calculate average peak area (at correct retention time and m/z 44) and
#   δ13C across two replicate measurements per sample

# - Correcting δ13C signatures for the added methyl group in each fatty acid
#   methyl ester (fatty acids are converted into fatty acid methyl esters (FAME)
#   to improve their volatility and chromatographic behaviour). With a known
#   δ13C of the methyl group and number of carbon atom per FAME molecule,
#   the correction happens as follows:
#   δ13C_FA = (#C_FAME * δ13C_FAME - #C_MeOH * δ13C_MeOH) / #C_FA

# - Correcting for chloroform (CHCl3) remaining in the final sample

# - Correcting for standards 13:0 and 19:0

# - Blanks for:
#   · 13:0
#   · i15:0
#   · a15:0
#   · 16:0
#   · i13:1w8 + 10Me16:0
#   · a17:0



plfa_conc_batch1 <-
  openxlsx::read.xlsx(paste0("data/raw_data/plfa/",
                             "PLFA_calculation_sheet_deltaC13_batch1_",
                             "relative.xlsx"),
                      sheet = "plfa nmol")






