
# Main script to run the analysis and make the figures for Lu et al. 2016 (Journal of Theoretical Biology).

# Load packages
# Please note that the "plotBy" package is only available from bitbucket.org/remkoduursma/plotby
source("R/load_packages.R")

# Functions for Figs. 6, 7, and 9
source("R/Bootstrap predict_nls.R")
source("R/Bootstrap functions.R")


# 1. Run the analysis, unless the results file already exist.
# NOTE: this takes ca. 3 days on a fast computer.
outf <- "output/data/derived variables.csv"
if(!file.exists(outf)){
  source("run_analysis.R")
}
dvs <- read.csv(outf)

# 2. Make the figures. 
# Output written to output/figures as PDFs.
# Additional simulation output written to output/data
source("R/Fig. 1.R")
source("R/Fig. 2.R")
source("R/Fig. 3.R")
source("R/Fig. 4.R") # ca. 30min on first run
source("R/Fig. 5.R")
source("R/Fig. 6.R")
source("R/Fig. 7.R")
source("R/Fig. 8.R")
source("R/Fig. 9.R")
