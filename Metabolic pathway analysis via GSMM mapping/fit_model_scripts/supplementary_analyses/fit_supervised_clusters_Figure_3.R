# Imports and options -----------------------------------------------------


rm(list = ls())

# Run locally, just call renv
library(renv)
renv::load(".")


library(lme4)
library(lmerTest)
library(emmeans)
library(tictoc)
library(optparse)
library(pbkrtest)
library(MuMIn) # Used for computing the r statistics

emm_options(lmerTest.limit = 5000)
emm_options(pbkrtest.limit = 5000)

# Set up configuration in English
Sys.setenv(LANG = "en") 


# Notes to the format of this file ----------------------------------------

# This was run with R version 4.2.1 and renv version 0.17.3. Some packages were compiled with R 4.2.3

# This script attempts to fit pathway abundance of the pathways shown in Figure 3c to the clusters that were found by the MaAsLin2 pipeline.
# It fits a regression pathway_abundance ~ 1 + Cluster

# Command line parser -----------------------------------------------------

option_list = list(
  make_option(c("-i", "--input"), type="character", 
              help="input file path", metavar="character"),
  make_option(c("-r", "--record"), type="logical", default=FALSE, 
              help="Whether or not to redirect console output to a log file.[default = FALSE]", metavar="character"),
  make_option(c("-f", "--folder_output"), type="character", default=file.path(getwd()), 
              help="output folder [default= current working directory. Results will then be saved in a folder named after the date of the run.", metavar="character"),
  make_option(c("-l", "--log"), type="character", default="log.txt", 
              help="log file name [default= log.txt]", metavar="character"),
  make_option(c("-b", "--beg"), type="double", 
              help="Column index of the first pathway to process. Should be 3 or more.", metavar="character"),
  make_option(c("-e", "--end"), type="double", 
              help="Column index of the last pathway to process. Should be 85 or less", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Set up output files -----------------------------------------------------

folder_output = file.path(opt$folder_output, format(Sys.time(), "%F//"))
dir.create(folder_output)

if (opt$record){
  zz <- file(file.path(paste0(folder_output,"_beg_", opt$beg, "_end_", opt$end, opt$log)), open = "wt")
  sink(zz ,type = "output")
  sink(zz, type = "message")
  on.exit(sink(, type = "output"))
  on.exit(sink(, type = "message"))
}

outfile_em = file.path(folder_output,paste0("_beg_", opt$beg, "_end_", opt$end, "_cluster_results.csv"))



# Prepare the data --------------------------------------------------------

data <- read.csv(opt$input, check.names = FALSE)
rxn_meta <- data # To avoid having to rewrite the script


# Curate the annotations ----
# Change the categorical variables to factors with specified reference levels
cluster_fact <- as.factor(rxn_meta$Cluster)

rxn_factors <- rxn_meta
rxn_factors$Cluster <- cluster_fact


# Prepare the output of the model fitting ---------------------------------

emmeans_res = data.frame(matrix(ncol = 31, nrow = 0)) # emmeans estimates

colnames(emmeans_res) = c( "1-2 est", "1-3 est", "1-4 est", "1-5 est", "2-3 est", "2-4 est", "2-5 est", "3-4 est", "3-5 est", "4-5 est",
                           "1-2 pval", "1-3 pval", "1-4 pval", "1-5 pval", "2-3 pval", "2-4 pval", "2-5 pval", "3-4 pval", "3-5 pval", "4-5 pval",
                           "1-2 std", "1-3 std", "1-4 std", "1-5 std", "2-3 std", "2-4 std", "2-5 std", "3-4 std", "3-5 std", "4-5 std", "r2")


# Model fitting -----------------------------------------------------------

pb <-txtProgressBar(min = opt$b, max = opt$e, style = 3)#
i = opt$b
for (rxn in colnames(rxn_factors[opt$b:opt$e])){

  setTxtProgressBar(pb, i)
  
  tic()
  rxn_col <- c(rxn, "Cluster")
  R_df <- rxn_factors[rxn_col]
  colnames(R_df) <- c("R", "Cluster")
  R_df$R <- as.numeric(R_df$R)
  print(colnames(rxn_factors)[i])
  # Fit the models
  lmm <- lm(R ~ 1 + Cluster,data = R_df)
  r2m <- r.squaredGLMM(lmm)[1, "R2m"]
  r2c <- r.squaredGLMM(lmm)[1, "R2c"]
  results <- summary(lmm, ddf = "Kenward-Roger")
  # Extract fields of interest
  coeff <- data.frame(results$coefficients)
  estimates <- coeff$Estimate
  stds <- coeff$Std..Error
  pvals <- coeff$Pr...t..
  
  # Then compute the pairwise contrasts for the clusters
  # Pvalues will be adjusted later as we will have to adjust for all pathways
  emm_model <- emmeans(lmm, "Cluster", lmer.df = "kenward-roger")
  em_sum <- summary(pairs(emm_model), adjust = "none") 
  em_est <- em_sum$estimate
  em_pval <- em_sum$p.value
  em_se <- em_sum$SE
  
  toc()
  
  em_res <- t(data.frame(c(em_est[1],em_est[2],em_est[3],em_est[4],em_est[5],em_est[6],em_est[7],em_est[8],em_est[9],em_est[10],
                           em_pval[1],em_pval[2],em_pval[3],em_pval[4],em_pval[5],em_pval[6],em_pval[7],em_pval[8],em_pval[9],em_pval[10],
                           em_se[1],em_se[2],em_se[3], em_se[4], em_se[5], em_se[6], em_se[7], em_se[8], em_se[9], em_se[10], r2c)))
  
  
  rownames(em_res) <- rxn
  colnames(em_res) <- colnames(emmeans_res)
  
  emmeans_res <- rbind(emmeans_res, em_res)
  print(i)
  i <- i + 1
  close(pb)
}

write.csv(emmeans_res, outfile_em)


if (opt$record) {
  sink(, type = "output")
  sink(, type = "message")
}
