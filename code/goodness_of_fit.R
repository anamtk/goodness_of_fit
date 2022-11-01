# Graphical posterior predictive check
# April 1, 2022
# Ana Miller-ter Kuile

# this is a hack of the bayesplot functionality to generate posterior
# predictive check graphs - specifically to assess - is the model family and link
# function I've selected appropriate for the data I have, or do I need to consider
# a different link or distribution (e.g. logit vs. cloglog link for binomial data; 
# poisson vs. negative binomial distribution for count data)

# Load packages ---------------------------------------------------------------

# Load packages, here and tidyverse for coding ease, 
package.list <- c("here", "tidyverse", 
                  "coda", "bayesplot",
                  "reshape2", "BayesPostEst")


## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% 
                                 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

theme_set(theme_bw())


# Load data ---------------------------------------------------------------

# we need the model output, where we tracked yrep
GOF_mcmc <- readRDS(here("model_outputs",
                      "GOF_mcmc.RDS"))

#and we also need our original y data
data <- read.csv(here("data",
                      "data.csv"))

#and we need the model with our parameter estimates
model_mcmc <- readRDS(here('model_outputs',
                           'model_mcmc.RDS'))



# Extract observed data from DF -------------------------------------------

#we need to extract our observed data from our dataframe
y <- data %>%
  #set fate to 1-=
  mutate(Fate_class = case_when(Fate_cat == "success" ~ 1,
                                Fate_cat == "failure" ~ 0,
                                TRUE ~ NA_real_)) %>%
  ungroup() %>%
  dplyr::select(Fate_class, Nest_ID) %>%
  #make this type "observed"
  mutate(type = "Observed") %>%
  #match the Nest_ID to the model output numbering
  mutate(Nest_ID = 1:n())


# Get yrep into DF format for graphing ------------------------------------

#extract the yreps, which for this model, which is an array of 
# iterations, nests, visits to nests, or a 3-D matrix
yreps <- GOF_mcmc$sims.list$yrep

#Subset a smaller set of these yreps, here I've selected 500
yreps <- yreps[1:500,,]

#Using the melt function from reshape2 package, turn the 3-D matrix
#into a dataframe with a column for iteration ID, nest ID, and interval ID
yrep <- melt(yreps) %>%
  #make this dataframe have 
  #1. a numbered "iteration" column for iteration #
  #2. a nest ID column 1:n nests
  #3. Interval # corresponds to the survey interval in the 
  ## repeated measures structure of the model
  #4. and then assigning the 1-0 to a fate class
  rename("Iteration" = "Var1",
         "Nest_ID" = "Var2",
         "Interval_num" = "Var3",
         "Fate_class" = "value") %>%
  #indicate these are the simulated data
  mutate(type = "Simulated")


# Graph observed versus simulated -----------------------------------------

#inspired from this: 
#https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1314
#posterior predictive check graphical observation
ggplot() +
  #graph the simulated data density per iteration
  geom_density(data = yrep, aes(x = Fate_class, group = Iteration, fill = type), 
               alpha = 0.2) +
  #graph the observed data density
  geom_density(data = y, aes(x = Fate_class, fill = type), alpha = 0.5)

#more ones than in our dataset are predicted by this model

# AUC ---------------------------------------------------------------------

#this section of script calculates AUC 

#getting a dataset with names that match model (I don't think
# this is necessary, but did it anyway)
nests2 <- data %>%
  mutate(y = case_when(Fate_cat == "success" ~ 1,
                       Fate_cat == "failure" ~ 0,
                       TRUE ~ NA_real_)) %>%
  mutate(Tmax2 = Tmax^2) %>%
  mutate(Trt_cat = factor(Trt_cat, levels = c("U", "B", "H", "HB"))) %>%
  mutate(Tree_sp = factor(Tree_sp, levels = c("PIPO", "Abies", "POTR5",
                                              "JUOC", "PSME"))) %>%
  mutate(Time_groups = factor(Time_groups, 
                              levels = c("oot", "0-3", "4-9", "10+"))) %>%
  mutate(prevStage = factor(prevStage,
                            levels = c("Ne", "Eg"))) %>%
  rename("StageID" = "prevStage",
         "TreatmentID" = "Trt_cat",
         "NestHt" = "Nest_Ht",
         "SpeciesID" = "Tree_sp",
         "InitDay" = "Init_day",
         "Trees50" = "Trees_50",
         "Trees2550" = "Trees_2550",
         "PercPonderosa" = "pPIPO",
         "ForestCV" = "a1000_areacv2",
         "Contag" = "a1000_contag",
         "OpenNm" = "a1000_np1",
         "LandBu" = "a1000_Bu",
         "LandHa" = "a1000_Ha")

colnames(nests2)
# a matrix of the data with all the variables that 
# go into the model
mod_matrix <- model.matrix(y ~  StageID + TreatmentID + 
                             n_tx + Time_groups + NestHt +
                             cosOrientation + SpeciesID +
                             InitDay + Trees50 + Trees2550 + 
                             PercPonderosa +
                             Trees50 * PercPonderosa +
                             Trees2550 * PercPonderosa + 
                             Trees50*Tmax +
                             Trees2550*Tmax + 
                             Tmax +  PPT + ForestCV + 
                             Contag + OpenNm + LandHa + LandBu +
                             LandHa*LandBu + Tmax*PPT,
                           data = nests2)

#a matrix of the data with all the *data* that will go into
# the model
x_data <- as.matrix(model.frame(y ~  StageID + TreatmentID + 
                                  n_tx + Time_groups + NestHt +
                                  cosOrientation + SpeciesID +
                                  InitDay + Trees50 + Trees2550 + 
                                  PercPonderosa +
                                  Trees50 * PercPonderosa +
                                  Trees2550 * PercPonderosa + 
                                  Trees50*Tmax +
                                  Trees2550*Tmax + 
                                  Tmax +  PPT + ForestCV + 
                                  Contag + OpenNm + LandHa + LandBu +
                                  LandHa*LandBu + Tmax*PPT,
                                data = nests2))

#defining the variables that we want to pull out of that for
# the mcmc matrix that matches the model and data matrices
vars <- c('b0', 'b1StageID',
          'b2TreatmentID', "b3TrtTime",
          "b4SpeciesID", 'b'
)
#getting just the simulated list of variables from the mmcmc:
mcmcout <- model_mcmc$sims.list
#selecting those in ou list of variables for the model
mcmcout <- mcmcout[names(mcmcout) %in% vars] 
#making this into a matrix, making sure that the 
# order of the columns matches that of the data and
# model matrices
mcmc_matrix <- as.data.frame(do.call(cbind, mcmcout)) %>%
  dplyr::select(b0, 
                V3, 
                V5:V7,
                V9:V11,
                V13:V16, V21:V40) %>%
  as.matrix()

# using mcmcRocPrcGen, determining the ROC curve for the 
#model
fit_sum1 <- mcmcRocPrcGen(modelmatrix = mod_matrix,
                          modelframe = x_data,
                          mcmcout = mcmc_matrix,
                          curves = TRUE,
                          fullsims = FALSE)

#getting the AUC stat from that output
fit_sum1$area_under_roc

#plotting that ROC curve
ggplot(data = as.data.frame(fit_sum1$roc_dat), 
       aes(x = x, y = y)) +
  geom_line() + 
  geom_abline(intercept = 0, slope = 1, color = "gray") + 
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") +
  annotate(geom = "text", x = 0.75, y = 0.25, label = "AUC = 0.58")


