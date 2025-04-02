# Imports and options -----------------------------------------------------


rm(list = ls())

# Necessary for running on the cluster, comment and adjust for your own system

library(renv, lib = "./Rlib_421")
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

# This script fit one reaction-specific linear mixed effect model (LMEM) or linear model per reaction.
# Each model takes ~2.30mn to fit and there are 7604 such models to fit
# It has been created to run in parallel on our cluster.
# Each run of this script create two files with different estimates and statistics
# The files are indexed by the column indices of the first and last reactions to be fit
# The 'results' file contains estimates and statistics from the LMEM model directly
# The 'emmeans' fils contains estimates and statistics for the pairwise contrasts estimated using emmeans
# The format of the output files is hardcoded in this script
# P-values are collected but not adjusted. Multiple testing adjustment ocurs in the analysis Python notebooks.

# This script also logs an extensive output to help follow what is happening

# The first reaction to process has index 2 as the column with index 1 is the #SampleID column
# The gender of the 2 pregnant mice is curated so that they show as female mice
# The 113 wild mice with unknown gender are eliminated from the analysis as we replace their gender by NA
# Those 113 mice correspond to the 113 cecum samples. Those are thus eliminated from the analysis, hence why the sample type is not used as a predictor in the model.
# See the lmerTest vignette for the default na.action behavior: https://cran.r-project.org/web/packages/lmerTest/lmerTest.pdf
# We also only process the samples for wild mice, SPF mice, and human samples

# Command line parser -----------------------------------------------------

option_list = list(
  make_option(c("-i", "--input"), type="character", 
              help="input file path", metavar="character"),
  make_option(c("-c", "--columns"), type="character", 
              help="file containing the column names referring to the reactions. [default: ../data/processed_files/rxn_names.csv]", default = "../data/processed_files/rxn_names.csv" ,metavar="character"),
  make_option(c("-r", "--record"), type="logical", default=FALSE, 
              help="Whether or not to redirect console output to a log file.[default = FALSE]", metavar="character"),
  make_option(c("-f", "--folder_output"), type="character", default=file.path(getwd()), 
              help="output folder [default= current working directory. Results will then be saved in a folder named after the date of the run.", metavar="character"),
  make_option(c("-l", "--log"), type="character", default="log.txt", 
              help="log file name [default= log.txt]", metavar="character"),
  make_option(c("-b", "--beg"), type="double", 
              help="Column index of the first reaction to process", metavar="character"),
  make_option(c("-e", "--end"), type="double", 
              help="Column index of the last reaction to process", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);




# Set up output files -----------------------------------------------------

folder_output = file.path(opt$folder_output, format(Sys.time(), "%F/"))
dir.create(folder_output)

if (opt$record){
  zz <- file(file.path(paste0(folder_output,"_beg_", opt$beg, "_end_", opt$end, opt$log), open = "wt"))
  sink(zz ,type = "output")
  sink(zz, type = "message")
  on.exit(sink(, type = "output"))
  on.exit(sink(, type = "message"))
}

outfile_res = file.path(folder_output,paste0("_beg_", opt$beg, "_end_", opt$end, "_results.csv"))
outfile_em = file.path(folder_output,paste0("_beg_", opt$beg, "_end_", opt$end, "_emmeans.csv"))



# Prepare the data --------------------------------------------------------

data <- read.csv(opt$input, check.names = FALSE)
rxns <- read.csv(opt$columns)
# Select the columns of interest
rxn_meta <- cbind(data[1], data[,rxns$X0], data$Model, data$Origin, data$Gender, data$Adult_Pups, data$IontorrentChip)
colnames(rxn_meta)<- c("#SampleID", rxns$X0, "Model", "Origin", "Gender", "Adult_Pups", "IontorrentChip")

# Curate the annotations ----

# Replace gender Female_Preg by Female
rxn_meta[rxn_meta == "Female_Preg"] <- "Female"
# Replace unknown values by NA so that they are dropped automatically when fitting the model
# NB: in create_input_dataset.ipynb we checked that this will not result in losing cretain reactions
rxn_meta[rxn_meta == "Unknown"] <- NA

# Select only the samples of interest ---
rxn_meta =  rxn_meta[rxn_meta$Model %in% c("Human", "SPF", "Wild"),]


# Change the categorical variables to factors with specified reference levels
model_fact <- as.factor(rxn_meta$Model)
model_fact <- relevel(model_fact, ref = "Wild")
origin_fact <- as.factor(rxn_meta$Origin)
gender_fact <- as.factor(rxn_meta$Gender)
age_fact <- as.factor(rxn_meta$Adult_Pups)
# Add a factor corresponding to Mouse/Human. TRUE means it is a human
species_fact <- as.factor(model_fact == "Human")
# Add a factor corresponding to the chip
chip_fact <- as.factor(rxn_meta$IontorrentChip)


rxn_factors <- rxn_meta
rxn_factors$Origin <- origin_fact
rxn_factors$Model <- model_fact
rxn_factors$Gender <- gender_fact
rxn_factors$AdultPups <- age_fact
rxn_factors$Species <- species_fact
rxn_factors$Chip <- chip_fact


# Prepare the output of the model fitting ---------------------------------

results_df = data.frame(matrix(ncol = 26, nrow = 0))
emmeans_res = data.frame(matrix(ncol = 9, nrow = 0)) # emmeans estimates

colnames(emmeans_res) = c( "Wild vs SPF est", "Wild vs Human est",  "SPF vs Human est",
                           "Wild vs SPF pval", "Wild vs Human pval",  "SPF vs Human pval",
                           "Wild vs SPF std", "Wild vs Human std",  "SPF vs Human std")
colnames(results_df) <- c("resid std", "Origin std", "Chip std",
                          "Intercept est", "ModelHuman est","ModelSPF est", "MouseGenderMale est",  "HumanGenderMale est", "MouseAdult_PupsPup est",
                          "Intercept std", "ModelHuman std", "ModelSPF std","MouseGenderMale std",  "HumanGenderMale std", "MouseAdult_PupsPup std",
                          "Intercept pval", "ModelHuman pval", "ModelSPF pval","MouseGenderMale pval",  "HumanGenderMale pval", "MouseAdult_PupsPup pval",
                          "fit msg", "conv msg", "israndom", "r2c", "r2m"
)

# Model fitting -----------------------------------------------------------

pb <-txtProgressBar(min = opt$b, max = opt$e, style = 3)#
i = opt$b
for (rxn in colnames(rxn_factors[opt$b:opt$e])){
  setTxtProgressBar(pb, i)
  
  tic()
  rxn_col <- c(rxn, "Model", "Origin", "Adult_Pups", "Gender", "Species", "Chip")
  R_df <- rxn_factors[rxn_col]
  colnames(R_df) <- c("R", "Model", "Origin", "Adult_Pups", "Gender", "Species", "Chip")
  R_df$R <- as.numeric(R_df$R)
  
  # Fit the models
  
  # To do model comparison we cannot used restricted maximum likelihood (REML) and need to set REML = FALSE
  # According to Luke, 2017 (DOI: 10.3758/s13428-016-0809-y)Sattherwaite's and Kenward-Roger estimation of the degress of freedom (DF) is not too anti-conservative when used with REML
  # Kenward-Roger takes longer but is compatible with emmeans so we will use this one
  # Step 1: fit model with and without random effects with ML then run a likelihood ratio test (LRT) to determine whether the random effects are significant
  # Step 2: choose the best model according to LRT, fit it with REML and use Kenward-Roger approximation of DF
  with_ranef <- lmer(R ~ 1 + Species + Species:Model + Species:Gender + Species:Adult_Pups + (1 | Species:Origin) + (1 | Chip), data = R_df, REML = FALSE) # REML = FALSE for LRT
  ranova_res <- ranova(with_ranef)
  
  LRT_res_origin <- ranova_res$`Pr(>Chisq)`[2]
  LRT_res_chip <- ranova_res$`Pr(>Chisq)`[3] # Null Hypothesis: the best model is the one without random effects. Hence we go for random effects if P<=0.05
  if (LRT_res_origin <= 0.05 & LRT_res_chip <= 0.05) { # The random effects on the origin and the chip are significant
    
    lmm <- lmer(R ~ 1 + Species + Species:Model + Species:Gender + Species:Adult_Pups + (1 | Species:Origin) + (1 | Chip), REML = TRUE, data = R_df)
    ran_origin <- as.data.frame(VarCorr(lmm))$sdcor[1]
    ran_chip <- as.data.frame(VarCorr(lmm))$sdcor[2]
    print("random origin AND random chip")
    results <- summary(lmm, ddf = "Kenward-Roger")
    israndom <- "OC"
    r2m <- r.squaredGLMM(lmm)[1, "R2m"]
    r2c <- r.squaredGLMM(lmm)[1, "R2c"]
  } else {
    if (LRT_res_origin <= 0.05 & LRT_res_chip > 0.05) { # The random effects on the origin are significant
      lmm <- lmer(R ~ 1 + Species + Species:Model + Species:Gender + Species:Adult_Pups + (1 | Species:Origin), REML = TRUE, data = R_df)
      ran_origin <- as.data.frame(VarCorr(lmm))$sdcor[1]
      ran_chip <- 0
      print("random origin")
      results <- summary(lmm, ddf = "Kenward-Roger")
      r2m <- r.squaredGLMM(lmm)[1, "R2m"]
      r2c <- r.squaredGLMM(lmm)[1, "R2c"]
      israndom <- "O"
    } else {
      if (LRT_res_origin > 0.05 & LRT_res_chip <= 0.05) { # The random effects on the chip are significant
        lmm <- lmer(R ~ 1 + Species + Species:Model + Species:Gender + Species:Adult_Pups + (1 | Chip), REML = TRUE, data = R_df)
        ran_chip <- as.data.frame(VarCorr(lmm))$sdcor[1]
        ran_origin <- 0
        print("random chip")
        results <- summary(lmm, ddf = "Kenward-Roger")
        r2m <- r.squaredGLMM(lmm)[1, "R2m"]
        r2c <- r.squaredGLMM(lmm)[1, "R2c"]
        israndom <- "C"
      } else {
        # The random effects are not significant
        
        lmm <- lm(R ~ 1 + Species + Species:Model + Species:Gender + Species:Adult_Pups, data = R_df, na.action = na.omit) # This is always a ML estimate, not a REML
        ran_origin <- 0
        ran_chip <- 0
        print(" no random effects")
        results <- summary(lmm) # ddf specification is not necessary here
        r2m <- r.squaredGLMM(lmm)[1, "R2m"]
        r2c <- r.squaredGLMM(lmm)[1, "R2c"]
        israndom <- "fixed"
      }
    }
  }
  
  
  
  # Extract fields of interest
  coeff <- data.frame(results$coefficients)
  estimates <- coeff$Estimate
  stds <- coeff$Std..Error
  pvals <- coeff$Pr...t..

  
  
  # Then compute the pairwise contrasts for the models
  # Do not adjust the pvalues yet, they will be adjusted in the analysis scripts
  emm_model <- emmeans(lmm, "Model", lmer.df = "kenward-roger")
  em_sum <- summary(pairs(emm_model), adjust = "none") 
  em_est <- em_sum$estimate
  em_pval <- em_sum$p.value
  em_se <- em_sum$SE
  
  toc()
  
  
  lmm_res <-  t(data.frame(c(results$sigma, ran_origin, ran_chip,
                             estimates[1],estimates[2], estimates[3], estimates[4], estimates[5], estimates[6],
                             stds[1],stds[2], stds[3], stds[4], stds[5],stds[6],
                             pvals[1],pvals[2], pvals[3], pvals[4], pvals[5],pvals[6], 
                             paste("E",results$fitMsgs), paste("E", results$optinfo$message), israndom, r2c, r2m)))
  
  em_res <- t(data.frame(c(em_est[1],em_est[2],em_est[3],
                           em_pval[1],em_pval[2],em_pval[3],
                           em_se[1],em_se[2],em_se[3])))
  
  
  rownames(lmm_res) <- rxn
  colnames(lmm_res) <- colnames(results_df)
  
  rownames(em_res) <- rxn
  colnames(em_res) <- colnames(emmeans_res)
  
  results_df <- rbind(results_df, lmm_res)
  emmeans_res <- rbind(emmeans_res, em_res)
  print(i)
  i <- i + 1
  close(pb)
}


write.csv(results_df, outfile_res)
write.csv(emmeans_res, outfile_em)


if (opt$record) {
  sink(, type = "output")
  sink(, type = "message")
}
