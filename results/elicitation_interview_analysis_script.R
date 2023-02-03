### Elicitation interview analysis for joint elicitation and modelling paper ###
### v 4.2.1 ###
### 03/02/2023 ###

#### Load packages ####
library(ggplot2)
library(tidyr)
library(nnet) # For multinomial models
library(tidyverse)
library(caret) # For stepwise multinomial model
library(doParallel) # For parallel processing
library(generalhoslem) # Hosmer-Lemeshow goodness-of-fit test

#### Functions ####

## Convert Likert-scale responses to numeric
AgreementScale <- function(x){
  levels(x) <- list("1"="Strongly disagree", "2"="Disagree", "3"="Neither agree nor disagree", "4"="Agree", "5"="Strongly agree")
  x <- x %>% as.character() %>% as.numeric
  return(x)
}

## Calculate p-values from multinomial logistic regression models
p_values <- function(x){
  z <- summary(x)$coefficients / summary(x)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  return(p)
}

## Automatic multinomial model building function
# x = Null model
# y = Dataframe with only columns containing the independent variables for consideration in the model - if y is not a model matrix, add "X" as a suffix to all column names
Multinom_model_build <- function(x,y){
  mod_list <- list()
  z <- 1
  mod_list[[z]] <- x
  
  repeat{
    ## Add in a significant (p<0.05) term
    # Test p-values of all terms when added to the model and save in a dataframe
    p_list <- NULL
    for(i in colnames(y)){
      mod = update(x, as.formula(paste(".~.+", i)))
      p = c(i, ifelse(ncol(p_values(x)) == ncol(p_values(mod)), NA, p_values(mod)[,c((ncol(p_values(x)) + 1):ncol(p_values(mod)))] %>% min()) )
      p_list <- rbind(p_list, p) 
    }
    p_list <- as.data.frame(p_list)
    # Make p-values numeric
    p_list$V2 <- as.character(p_list$V2) %>% as.numeric()
    # Order p-values in ascending size
    p_list <- arrange(p_list, V2)
    # Remove non-significant p-values
    p_list <- subset(p_list, V2<0.05)
    
    # Stop function if there are no terms to add in (i.e. no terms in p_list)
    if(nrow(p_list)==0){
      break
    }
    
    # Select the term with the lowest p-value and add to model
    new_term <- p_list[1,1]
    x <- update(x, as.formula(paste(".~.+", new_term)) )
    z <- z + 1
    mod_list[[z]] <- x
    
    ## Remove terms that are no longer significant (p>=0.05)
    repeat{
      # Extract model p-values into a dataframe
      Mod_results <- p_values(x) %>% t() %>% as.data.frame()
      Mod_results$variable <- rownames(Mod_results)
      rownames(Mod_results) <- NULL
      # Add a column of the lowest P-value for each factor
      Mod_results$P <- apply(Mod_results[,1:(ncol(Mod_results) - 1)], 1, FUN = min)
      # Remove unnecessary columns & rows
      Mod_results <- Mod_results[2:(nrow(Mod_results)),(ncol(Mod_results) - 1):ncol(Mod_results)]
      # Remove factor level from each variable's name
      Mod_results$variable <- gsub("X.*", "X", Mod_results$variable) 
      # Remove non-significant terms
      Sig_vars <- subset(Mod_results[-1,], P<0.05)
      # Create a vector of significant terms
      Sig_vars <- Sig_vars$variable
      # Create a vector of non-significant terms to delete
      Non_sig_vars <- subset(Mod_results[,], P>=0.05)
      Non_sig_vars <- subset(Non_sig_vars, !Non_sig_vars$variable %in% Sig_vars)
      Non_sig_vars <- Non_sig_vars$variable
      Non_sig_vars <- Non_sig_vars[!Non_sig_vars %in% "(Intercept)"]
      # Delete non-significant variables
      for(i in Non_sig_vars){
        x = update(x, as.formula(paste(".~.-", i)) )
      }
      if(length(Non_sig_vars)==0){
        break
      }
      z <- z + 1
      mod_list[[z]] <- x
    }
    # Stop function if the model terms are identical to a previous iteration in the model build
    k <- 0
    for (i in 1:(length(mod_list)-1)){
      if(identical(sort(mod_list[[i]]$coefnames), sort(mod_list[[z]]$coefnames))==TRUE){
        k <- 1
      }
    }
    if(k==1){
      break
    }
  }
  return(x)
}

## Function to merge dataframes in list
list_merge <- function(df1, df2){
  x <- merge(df1, df2, by = "row.names", all = TRUE)
  rownames(x) <- x$Row.names
  x$Row.names <- NULL
  return(x)
}

## Function to caculate standard error
standard_error <- function(x){
  sd(x) / sqrt(length(x))
}

#### Load data ####
Data <- read.csv("./elicitation_interview_data.csv", fileEncoding="UTF-8-BOM")

# Make character columns factor
Data[sapply(Data, is.character)] <- lapply(Data[sapply(Data, is.character)], 
                                           as.factor)
# Make number of breeding cows columns
Data$Beef.cow.size[is.na(Data$Beef.cow.size)] <- 0
Data$Dairy.size[is.na(Data$Dairy.size)] <- 0
Data$N.breeding.cows <- Data$Beef.cow.size + Data$Dairy.size

# Reorder farmer age
levels(Data$Age) <- list("Under 30" = "Under 30", "30 - 39" = "30 - 39", "40 - 49" = "40 - 49", "50 - 59" = "50 - 59", "60 - 69" = "60 - 69", "Over 70" = "Over 70")
# Continuous age
Data$Age_continuous <- Data$Age
levels(Data$Age_continuous) <- list("25" = "Under 30", "35" = "30 - 39", "45" = "40 - 49", "55" = "50 - 59", "65" = "60 - 69", "75" = "Over 70")
Data$Age_continuous <- Data$Age_continuous %>% as.character() %>% as.numeric

# COM-B factors
Data$Psychological_capability <- (as.numeric(Data$I.know.how.to.control.infectious.disease.in.my.cattle) +
                                    as.numeric(Data$I.know.why.it.is.important.to.control.infectious.disease.in.my.cattle) +
                                    as.numeric(Data$I.understand.most.advice.I.receive.about.infectious.disease.in.cattle)) / 3
Data$Physical_opportunity <- ((3 + (3 - as.numeric(Data$I.do.not.have.the.time.to.control.infectious.disease.in.my.cattle))) +
                                (3 + (3 - as.numeric(Data$Controlling.infectious.disease.costs.too.much.money)))) / 2
Data$Social_opportunity <- (as.numeric(Data$Most.farmers.I.know.are.controlling.infectious.disease.in.their.cattle) +
                              as.numeric(Data$My.vet.helps.me.control.infectious.diseases.in.my.cattle) +
                              as.numeric(Data$Government.policy.helps.me.control.infectious.disease.in.my.cattle) +
                              as.numeric(Data$Other.farmers.help.me.achieve.control.of.infectious.disease.in.my.cattle) +
                              (3 + (3 - as.numeric(Data$I.find.it.difficult.to.raise.the.subject.of.infectious.disease.in.my.cattle.with.other.farmers))) +
                              (3 + (3 - as.numeric(Data$I.find.it.difficult.to.raise.the.subject.of.infectious.disease.in.my.cattle.with.vets))) +
                              (3 + (3 - as.numeric(Data$I.find.it.difficult.to.raise.the.subject.of.infectious.disease.in.my.cattle.with.governmental.organisations)))) / 7
Data$Automatic_motivation <- (as.numeric(Data$I.worry.about.getting.infectious.diseases.in.my.cattle) +
                                as.numeric(Data$I.feel.good.about.myself.when.I.control.infectious.disease.in.my.cattle) +
                                as.numeric(Data$I.want.to.control.infectious.diseases.for.the.sake.of.my.herd) +
                                as.numeric(Data$Controlling.infectious.diseases.in.my.cattle.is.part.of.my.routine)) / 4
Data$Reflective_motivation <- (as.numeric(Data$I.have.infectious.disease.goals.that.I.want.to.achieve.for.my.herd) +
                                 as.numeric(Data$I.want.to.control.infectious.disease.in.my.cattle) +
                                 as.numeric(Data$There.are.many.benefits.to.controlling.infectious.disease.in.cattle) +
                                 (3 + (3 - as.numeric(Data$It.is.not.my.responsibility.to.control.infectious.disease.in.my.cattle))) +
                                 as.numeric(Data$I.have.a.plan.for.controlling.infectious.disease.in.my.cattle)) / 5

# Standardise the factors
Psychosocial_COMB_factors <- c("When.dealing.with.farmers.it.is.better.to.be.careful.before.you.trust.them",
                               "I.feel.respected.by.the.government",
                               "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds",
                               "I.trust.other.farmers.I.meet.for.the.first.time",
                               "When.dealing.with.vets.it.is.better.to.be.careful.before.you.trust.them",
                               "When.dealing.with.strangers.it.is.better.to.be.careful.before.you.trust.them",
                               "In.general..one.can.trust.people",
                               "I.feel.respected.by.my.vet",
                               "I.trust.my.vet.s.advice.about.infectious.disease.control.in.my.herd",
                               "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle",
                               "I.feel.respected.by.the.veterinary.profession",
                               "I.trust.vets",
                               "When.dealing.with.the.Government.it.is.better.to.be.careful.before.you.trust.them",
                               "I.trust.beef.farmers",
                               "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession",
                               "My.vet.would.always.tell.me.the.truth.even.if.it.was.not.what.I.wanted.to.hear",
                               "I.trust.dairy.farmers",
                               "I.trust.governmental.organisations",
                               "I.trust.my.neighbours.to.be.controlling.infectious.diseases.in.their.herds",
                               "PP.vet",
                               "PP.vet.community",
                               "PP.neighbouring.farmers",
                               "PP.farming.community",
                               "PP.gov",
                               "PP.cows",
                               "PP.dairy.farmers",
                               "PP.beef.farmers",
                               "Psychological_capability",
                               "Physical_opportunity",
                               "Social_opportunity",
                               "Automatic_motivation",
                               "Reflective_motivation")
Norm_data <- scale(Data[,Psychosocial_COMB_factors]) %>% as.data.frame()
rm(Psychosocial_COMB_factors)
colnames(Norm_data) <- paste(colnames(Norm_data),"norm",sep="_")
Data <- cbind(Data, Norm_data)
Data$June.size_norm <- scale(Data$June.size)

## Create dataset for models
Model_data <- Data[,c("Vaccinate",
                      "June.size", "Vet.visit.freq.yr", "Age_continuous", "Financial.pressures", "N.breeding.cows",
                      "When.dealing.with.farmers.it.is.better.to.be.careful.before.you.trust.them", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds", "I.trust.other.farmers.I.meet.for.the.first.time", "I.trust.beef.farmers", "I.trust.dairy.farmers", "I.trust.my.neighbours.to.be.controlling.infectious.diseases.in.their.herds",
                      "When.dealing.with.vets.it.is.better.to.be.careful.before.you.trust.them", "I.feel.respected.by.my.vet", "I.trust.my.vet.s.advice.about.infectious.disease.control.in.my.herd", "I.feel.respected.by.the.veterinary.profession", "I.trust.vets", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession", "My.vet.would.always.tell.me.the.truth.even.if.it.was.not.what.I.wanted.to.hear",
                      "I.feel.respected.by.the.government", "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle", "When.dealing.with.the.Government.it.is.better.to.be.careful.before.you.trust.them", "I.trust.governmental.organisations",
                      "When.dealing.with.strangers.it.is.better.to.be.careful.before.you.trust.them", "In.general..one.can.trust.people",
                      "PP.neighbouring.farmers", "PP.farming.community", "PP.dairy.farmers", "PP.beef.farmers",
                      "PP.vet", "PP.vet.community",
                      "PP.gov",
                      "PP.cows",
                      "Psychological_capability", "Physical_opportunity", "Social_opportunity", "Automatic_motivation", "Reflective_motivation")]

## Make responses for each model - full dataset
# Multinomial factor response
Model_data$Vaccinate_multinom <- as.factor(Model_data$Vaccinate)
Model_data$Vaccinate_multinom_exc_never <- Model_data$Vaccinate_multinom
levels(Model_data$Vaccinate_multinom) <- list("S1-2" = c("1", "2"), "S3-5" = c("3", "4", "5"), "S6-never" = c("6", "7", "8", "9"))
levels(Model_data$Vaccinate_multinom_exc_never) <- list("S1-2" = c("1", "2"), "S3-5" = c("3", "4", "5"), "S6-8" = c("6", "7", "8"), "Never" = "9")

#### Create bootstrapping data ####
#### Preparation for bootstrapping

# Create a list of 500 bootstrap samples
Boot_list <- list()
set.seed(1341)
for(i in 1:500){
  Boot_sample <- c(sample(1:nrow(Model_data), size = nrow(Model_data), replace = TRUE))
  Boot_dataset <- Model_data[Boot_sample,]
  Boot_list[[i]] <- Boot_dataset
}
rm(Boot_sample, Boot_dataset, i)

# Make random permutations and bootstrapping reproducible
Random_boot_list <- list()
set.seed(1359)
# Create random permutations of the dataset
# Create vector of number of permuted datasets
Perm_data <- 1:10
Ran_boots <- 1:20
for(j in Perm_data){
  # Create dataset with covariates severed from Y
  Ys <- c("Vaccinate", "Vaccinate_multinom")
  Random_model_data <- cbind(Model_data[sample(1:nrow(Model_data)), Ys, drop = FALSE], Model_data[, !(colnames(Model_data) %in% Ys)])
  # Reorder columns to match Model_data
  Random_model_data <- Random_model_data[, colnames(Model_data)]
  
  # Bootstrap each permuted dataset
  for(i in Ran_boots){
    Random_boot_sample <- c(sample(1:nrow(Random_model_data), size = nrow(Random_model_data), replace = TRUE))
    Random_boot_dataset <- Random_model_data[Random_boot_sample,]
    Random_boot_list[[length(Random_boot_list)+1]] <- Random_boot_dataset
  }  
}  
rm(i, j, Ys, Random_boot_sample, Random_model_data, Random_boot_dataset)

#### Multinomial analysis ####
# Select variables for multinomial model
Multinom_mod_data <- Model_data[,c("Vaccinate_multinom",
                                   "June.size", "Vet.visit.freq.yr", "Age_continuous", "Financial.pressures", "N.breeding.cows",
                                   "When.dealing.with.farmers.it.is.better.to.be.careful.before.you.trust.them", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds", "I.trust.other.farmers.I.meet.for.the.first.time", "I.trust.beef.farmers", "I.trust.dairy.farmers", "I.trust.my.neighbours.to.be.controlling.infectious.diseases.in.their.herds",
                                   "When.dealing.with.vets.it.is.better.to.be.careful.before.you.trust.them", "I.feel.respected.by.my.vet", "I.trust.my.vet.s.advice.about.infectious.disease.control.in.my.herd", "I.feel.respected.by.the.veterinary.profession", "I.trust.vets", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession", "My.vet.would.always.tell.me.the.truth.even.if.it.was.not.what.I.wanted.to.hear",
                                   "I.feel.respected.by.the.government", "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle", "When.dealing.with.the.Government.it.is.better.to.be.careful.before.you.trust.them", "I.trust.governmental.organisations",
                                   "When.dealing.with.strangers.it.is.better.to.be.careful.before.you.trust.them", "In.general..one.can.trust.people",
                                   "PP.neighbouring.farmers", "PP.farming.community", "PP.dairy.farmers", "PP.beef.farmers",
                                   "PP.vet", "PP.vet.community",
                                   "PP.gov",
                                   "PP.cows",
                                   "Psychological_capability", "Physical_opportunity", "Social_opportunity", "Automatic_motivation", "Reflective_motivation")]
# Add an X suffix to columns for use in the model building function
colnames(Multinom_mod_data) <- paste(colnames(Multinom_mod_data), "X", sep = "")

# Run model
Multinom_mod <- Multinom_model_build(multinom(Vaccinate_multinomX ~ 1, data = Multinom_mod_data), Multinom_mod_data[,2:ncol(Multinom_mod_data)])

#### Test model fit

# List of variables in full models
Multinom_mod_variables <- all.vars(as.formula(Multinom_mod))

Multinom_mod_fit_data <- Multinom_mod_data[, Multinom_mod_variables]
Multinom_mod_fit_data <- Multinom_mod_fit_data[complete.cases(Multinom_mod_fit_data), ]

# Predict vaccination class from model
Multinom_mod_test_data <- Multinom_mod_fit_data[, Multinom_mod_variables[-1]]

Multinom_mod_fit_data$Predict_vaccinate_class <- predict(Multinom_mod, Multinom_mod_test_data, type = "class")

# What proportion of predictions are correct?
Correct_prediction <- nrow(subset(Multinom_mod_fit_data, Vaccinate_multinomX == Predict_vaccinate_class)) / nrow(Multinom_mod_fit_data)

# Where are incorrect predictions going?
Vaccinate_class_predictions <- with(Multinom_mod_fit_data, table(Predict_vaccinate_class, Vaccinate_multinomX))
# write.csv(Vaccinate_class_predictions, "./Vaccinate_class_predictions_multinomial_model.csv", na = "", row.names = TRUE)

# Predictions from cross-validated model
CV_metrics <- NULL
set.seed(1031)
for(i in 1:10){
  CV_predicted_data <- NULL
  cv_splits <- createFolds(Multinom_mod_data$Vaccinate_multinomX, k=10)
  for(j in cv_splits){
    Model <- multinom(Vaccinate_multinomX ~ June.sizeX + I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattleX + I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herdsX + PP.vetX, data = Multinom_mod_data[-j,])
    CV_predictions <- Model %>% predict(Multinom_mod_data[j,], type = "class")
    CV_predictions <- cbind(Multinom_mod_data[j,]$Vaccinate_multinomX, CV_predictions)
    CV_predicted_data <- rbind(CV_predicted_data, CV_predictions)
  }
  colnames(CV_predicted_data) <- c("Vaccinate_multinom", "CV_predictions")
  CV_metrics <- rbind(CV_metrics, nrow(subset(as.data.frame(CV_predicted_data), Vaccinate_multinom == CV_predictions)) / nrow(as.data.frame(CV_predicted_data)))
}

# Hosmer-Lemeshow goodness-of-fit
HLgof_test <- logitgof(Multinom_mod_data$Vaccinate_multinomX, fitted(Multinom_mod), g = 10, ord = FALSE)
Observed_data <- HLgof_test$observed
Observed_data$Data <- "Observed"
colnames(Observed_data) <- c("Decile", "S1-2", "S3-5", "S6-never", "Data")
Observed_data$Decile <- 1:10
Expected_data <- HLgof_test$expected
Expected_data$Data <- "Expected"
colnames(Expected_data) <- c("Decile", "S1-2", "S3-5", "S6-never", "Data")
Expected_data$Decile <- 1:10
Decile_data <- rbind(Observed_data, Expected_data)
rm(Observed_data, Expected_data)
Decile_data[,2:4] <- round(Decile_data[,2:4])
Decile_data <- pivot_longer(Decile_data, cols = c("S1-2", "S3-5", "S6-never"))
# Rename factors for plot
Decile_data$name <- as.factor(Decile_data$name)
levels(Decile_data$name) <- list("1 - 2" = "S1-2", "3 - 5" = "S3-5", "6 - Never" = "S6-never")

#png(filename="Hosmer_Lemeshow_deciles.png", res=600, width=6000, height=2500)
ggplot(Decile_data, aes(x = Data, y = value)) +
  geom_bar(aes(fill = name), position = "stack", stat = "identity") +
  facet_grid( ~ Decile) +
  scale_fill_manual(values = c("#F0E442", "#E69F00", "#D55E00")) +
  labs(y = "Number of farmers", fill = "Stage of epidemic\nvaccinated") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
#dev.off()

#### Bootstrapping

# Register parallel backend
registerDoParallel()

# Bootstrap model
Multinom_boot_coefs <- foreach(i=1:length(Boot_list)) %dopar% {
  # Load packages
  library(tidyverse)
  library(nnet)
  # Set up variables for models
  Boot_data <- Boot_list[[i]][,c("Vaccinate_multinom",
                                     "June.size", "Vet.visit.freq.yr", "Age_continuous", "Financial.pressures", "N.breeding.cows",
                                     "When.dealing.with.farmers.it.is.better.to.be.careful.before.you.trust.them", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds", "I.trust.other.farmers.I.meet.for.the.first.time", "I.trust.beef.farmers", "I.trust.dairy.farmers", "I.trust.my.neighbours.to.be.controlling.infectious.diseases.in.their.herds",
                                     "When.dealing.with.vets.it.is.better.to.be.careful.before.you.trust.them", "I.feel.respected.by.my.vet", "I.trust.my.vet.s.advice.about.infectious.disease.control.in.my.herd", "I.feel.respected.by.the.veterinary.profession", "I.trust.vets", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession", "My.vet.would.always.tell.me.the.truth.even.if.it.was.not.what.I.wanted.to.hear",
                                     "I.feel.respected.by.the.government", "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle", "When.dealing.with.the.Government.it.is.better.to.be.careful.before.you.trust.them", "I.trust.governmental.organisations",
                                     "When.dealing.with.strangers.it.is.better.to.be.careful.before.you.trust.them", "In.general..one.can.trust.people",
                                     "PP.neighbouring.farmers", "PP.farming.community", "PP.dairy.farmers", "PP.beef.farmers",
                                     "PP.vet", "PP.vet.community",
                                     "PP.gov",
                                     "PP.cows",
                                     "Psychological_capability", "Physical_opportunity", "Social_opportunity", "Automatic_motivation", "Reflective_motivation")]
  # Add an X suffix to columns for use in the model building function
  colnames(Boot_data) <- paste(colnames(Boot_data), "X", sep = "")
  X <- multinom(Vaccinate_multinomX ~ 1, data = Boot_data)
  Y <- Boot_data[,2:ncol(Boot_data)]
  # Run model
  Boot_model <- Multinom_model_build(X, Y)
  # Extract coefs & ps
  Boot_model <- cbind(summary(Boot_model)$coefficients %>% as.data.frame() %>% t(), p_values(Boot_model) %>% as.data.frame() %>% t())
  colnames(Boot_model) <- c(paste(colnames(Boot_model)[1:(ncol(Boot_model)/2)], "Coef", i, sep = "_"), paste(colnames(Boot_model)[((ncol(Boot_model)/2)+1):ncol(Boot_model)], "P", i, sep = "_"))
  return(Boot_model)
}
Multinom_boot_coefs <- Reduce(list_merge, Multinom_boot_coefs)
# Remove intercept row
Multinom_boot_coefs <- Multinom_boot_coefs[!rownames(Multinom_boot_coefs) %in% "(Intercept)",]
# Remove X from rownames
row.names(Multinom_boot_coefs) <- gsub("X", "", row.names(Multinom_boot_coefs))
# Drop P-values
Multinom_boot_coefs <- Multinom_boot_coefs[grepl("Coef", names(Multinom_boot_coefs), fixed=TRUE)]

Multinom_boot_coefs_level_2 <- Multinom_boot_coefs[,c(seq(1,1000,2))]
#write.csv(Multinom_boot_coefs_level_2, "Multinomial_boot_coefs_level_2.csv", row.names = TRUE)
Multinom_boot_coefs_level_3 <- Multinom_boot_coefs[,c(seq(2,1000,2))]
#write.csv(Multinom_boot_coefs_level_3, "Multinomial_boot_coefs_level_3.csv", row.names = TRUE)
#write.csv(Multinom_boot_coefs, "Multinomial_boot_coefs.csv", row.names = TRUE)

## Calculate covariate stability
Multinom_boot_vars <- Multinom_boot_coefs
# Make coefficients 1/0
Multinom_boot_vars[!is.na(Multinom_boot_vars)] <- 1
Multinom_boot_vars[is.na(Multinom_boot_vars)] <- 0
# Remove duplicate columns
col_odd <- seq_len(ncol(Multinom_boot_vars)) %% 2
Multinom_boot_vars <- Multinom_boot_vars[,col_odd == 1]
rm(col_odd)
Multinom_boot_stability <- Multinom_boot_vars
Multinom_boot_stability$Stability <- rowSums(Multinom_boot_stability) / ncol(Multinom_boot_stability)
Multinom_boot_stability <- arrange(Multinom_boot_stability, desc(Stability))
Multinom_boot_stability$Stability_rank <- 1:nrow(Multinom_boot_stability)
Multinom_boot_stability <- Multinom_boot_stability[,c("Stability", "Stability_rank")]
#write.csv(Multinom_boot_stability, "Multinomial_covariate_stabilities.csv", row.names = TRUE)

#### Bootstrapping on randomly permuted y

# Register parallel backend
registerDoParallel()

# Bootstrap model
Multinom_random_boot_coefs <- foreach(i=1:length(Random_boot_list)) %dopar% {
  # Load packages
  library(tidyverse)
  library(nnet)
  # Set up variables for models
  Random_boot_data <- Random_boot_list[[i]][,c("Vaccinate_multinom",
                                               "June.size", "Vet.visit.freq.yr", "Age_continuous", "Financial.pressures", "N.breeding.cows",
                                               "When.dealing.with.farmers.it.is.better.to.be.careful.before.you.trust.them", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds", "I.trust.other.farmers.I.meet.for.the.first.time", "I.trust.beef.farmers", "I.trust.dairy.farmers", "I.trust.my.neighbours.to.be.controlling.infectious.diseases.in.their.herds",
                                               "When.dealing.with.vets.it.is.better.to.be.careful.before.you.trust.them", "I.feel.respected.by.my.vet", "I.trust.my.vet.s.advice.about.infectious.disease.control.in.my.herd", "I.feel.respected.by.the.veterinary.profession", "I.trust.vets", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession", "My.vet.would.always.tell.me.the.truth.even.if.it.was.not.what.I.wanted.to.hear",
                                               "I.feel.respected.by.the.government", "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle", "When.dealing.with.the.Government.it.is.better.to.be.careful.before.you.trust.them", "I.trust.governmental.organisations",
                                               "When.dealing.with.strangers.it.is.better.to.be.careful.before.you.trust.them", "In.general..one.can.trust.people",
                                               "PP.neighbouring.farmers", "PP.farming.community", "PP.dairy.farmers", "PP.beef.farmers",
                                               "PP.vet", "PP.vet.community",
                                               "PP.gov",
                                               "PP.cows",
                                               "Psychological_capability", "Physical_opportunity", "Social_opportunity", "Automatic_motivation", "Reflective_motivation")]
  # Add an X suffix to columns for use in the model building function
  colnames(Random_boot_data) <- paste(colnames(Random_boot_data), "X", sep = "")
  X <- multinom(Vaccinate_multinomX ~ 1, data = Random_boot_data)
  Y <- Random_boot_data[,2:ncol(Random_boot_data)]
  # Run model
  Random_boot_model <- Multinom_model_build(X, Y)
  # Extract coefs & ps
  Random_boot_model <- cbind(summary(Random_boot_model)$coefficients %>% as.data.frame() %>% t(), p_values(Random_boot_model) %>% as.data.frame() %>% t())
  colnames(Random_boot_model) <- c(paste(colnames(Random_boot_model)[1:(ncol(Random_boot_model)/2)], "Coef", i, sep = "_"), paste(colnames(Random_boot_model)[((ncol(Random_boot_model)/2)+1):ncol(Random_boot_model)], "P", i, sep = "_"))
  return(Random_boot_model)
}
Multinom_random_boot_coefs <- Reduce(list_merge, Multinom_random_boot_coefs)
# Remove intercept row
Multinom_random_boot_coefs <- Multinom_random_boot_coefs[!rownames(Multinom_random_boot_coefs) %in% "(Intercept)",]
# Remove X from rownames
row.names(Multinom_random_boot_coefs) <- gsub("X", "", row.names(Multinom_random_boot_coefs))
# Drop P-values
Multinom_random_boot_coefs <- Multinom_random_boot_coefs[grepl("Coef", names(Multinom_random_boot_coefs), fixed=TRUE)]

#write.csv(Multinom_random_boot_coefs, "Multinomial_random_boot_coefs.csv", row.names = TRUE)

# Calculate covariate stability
Multinom_random_boot_vars <- Multinom_random_boot_coefs
# Make coefficients 1/0
Multinom_random_boot_vars[!is.na(Multinom_random_boot_vars)] <- 1
Multinom_random_boot_vars[is.na(Multinom_random_boot_vars)] <- 0
# Remove duplicate columns
col_odd <- seq_len(ncol(Multinom_random_boot_vars)) %% 2
Multinom_random_boot_vars <- Multinom_random_boot_vars[,col_odd == 1]
rm(col_odd)
Multinom_random_boot_stability <- Multinom_random_boot_vars

for(j in Perm_data){
  Multinom_random_boot_stability[,((max(Perm_data)*max(Ran_boots))+j)] <- rowSums(Multinom_random_boot_stability[,((max(Ran_boots)*(j-1))+1):(max(Ran_boots)*j)]) / max(Ran_boots)
}
rm(j)
Multinom_random_boot_stability <- Multinom_random_boot_stability[,201:ncol(Multinom_random_boot_stability)]
# Create table of stability thresholds
Stability_threshold <- NULL
for(i in 1:ncol(Multinom_random_boot_stability)){
  ecdff <- ecdf(Multinom_random_boot_stability[,i])
  Thresholds <- NULL
  for(j in seq(1, 0, -0.01)){
    quant_perm <- quantile(ecdff, probs = j)
    Thresholds <- c(Thresholds, quant_perm)
  }
  Stability_threshold <- cbind(Stability_threshold, Thresholds)
}
Multinomial_stability_thresholds <- cbind(rowMeans(Stability_threshold),
                                          apply(Stability_threshold, 1, standard_error))
colnames(Multinomial_stability_thresholds) <- c("Stability_threshold", "SE")
#write.csv(Multinomial_stability_thresholds, "Multinomial_stability_thresholds.csv")

#### Bootstrapping with 5 most stable covariates for calculating coefficients

# Register parallel backend
registerDoParallel()

# Bootstrap model
Multinom_boot_coefs_5_stable <- foreach(i=1:length(Boot_list)) %dopar% {
  # Load packages
  library(tidyverse)
  library(nnet)
  # Set up variables for models
  Boot_data <- Boot_list[[i]][,c("Vaccinate_multinom", "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle", "Physical_opportunity", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds", "June.size")]
  # Add an X suffix to columns for use in the model building function
  colnames(Boot_data) <- paste(colnames(Boot_data), "X", sep = "")
  X <- multinom(Vaccinate_multinomX ~ 1, data = Boot_data)
  Y <- Boot_data[,2:ncol(Boot_data)]
  # Run model
  Boot_model <- Multinom_model_build(X, Y)
  # Extract coefs & ps
  Boot_model <- cbind(summary(Boot_model)$coefficients %>% as.data.frame() %>% t(), p_values(Boot_model) %>% as.data.frame() %>% t())
  colnames(Boot_model) <- c(paste(colnames(Boot_model)[1:(ncol(Boot_model)/2)], "Coef", i, sep = "_"), paste(colnames(Boot_model)[((ncol(Boot_model)/2)+1):ncol(Boot_model)], "P", i, sep = "_"))
  return(Boot_model)
}
Multinom_boot_coefs_5_stable <- Reduce(list_merge, Multinom_boot_coefs_5_stable)
# Remove intercept row
Multinom_boot_coefs_5_stable <- Multinom_boot_coefs_5_stable[!rownames(Multinom_boot_coefs_5_stable) %in% "(Intercept)",]
# Remove X from rownames
row.names(Multinom_boot_coefs_5_stable) <- gsub("X", "", row.names(Multinom_boot_coefs_5_stable))
# Drop P-values
Multinom_boot_coefs_5_stable <- Multinom_boot_coefs_5_stable[grepl("Coef", names(Multinom_boot_coefs_5_stable), fixed=TRUE)]

Multinom_boot_coefs_5_stable_level_2 <- Multinom_boot_coefs_5_stable[,c(seq(1,1000,2))]
#write.csv(Multinom_boot_coefs_5_stable_level_2, "Multinomial_boot_coefs_5_stable_level_2.csv", row.names = TRUE)
Multinom_boot_coefs_5_stable_level_3 <- Multinom_boot_coefs_5_stable[,c(seq(2,1000,2))]
#write.csv(Multinom_boot_coefs_5_stable_level_3, "Multinomial_boot_coefs_5_stable_level_3.csv", row.names = TRUE)
#write.csv(Multinom_boot_coefs_5_stable, "Multinomial_boot_coefs_5_stable.csv", row.names = TRUE)

# Calculate mean coefficients
Boot_coefs_5_stable_means <- cbind(rowMeans(Multinom_boot_coefs_5_stable_level_2, na.rm = TRUE) %>%
  exp() %>%
  round(digits = 2),
rowMeans(Multinom_boot_coefs_5_stable_level_3, na.rm = TRUE) %>%
  exp() %>%
  round(digits = 2))
colnames(Boot_coefs_5_stable_means) <- c("Mid_OR", "Late_OR")
#write.csv(Boot_coefs_5_stable_means, "Multinomial_boot_coefs_5_stable_means.csv", row.names = TRUE)

#### Investigate clustering of variables of interest ####

#### Stabilised model

## Clustering on psychosocial responses
cluster_colnames <- c("I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle_norm", "Physical_opportunity_norm", "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession_norm", "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds_norm", "June.size_norm")

## K-means clustering
## Use map_dbl to run many models with varying value of k (centers)
set.seed(1612)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = Data[,cluster_colnames], centers = k)
  model$tot.withinss
})
## Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:10,
  tot_withinss = tot_withinss
)
## Plot the elbow plot
#png(filename="Five_covariate_clustering_elbow_plot.png", res=600, width=3000, height=2000)
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:10) +
  labs(y = "Total within-cluster sum of squares") +
  theme_bw()
#dev.off()  

## 3 clusters
set.seed(1343)
k_means_cluster_5_var <- kmeans(x = Data[,cluster_colnames], centers = 3)
# Add cluster to dataset
Data$stab_mod_clust <- k_means_cluster_5_var$cluster %>% as.factor()
Model_data$stab_mod_clust <- k_means_cluster_5_var$cluster %>% as.factor()

#### Plot relationships between clusters & variables
# Plot scenario when vaccination is used & k means cluster
table(Model_data$Vaccinate_multinom, Model_data$stab_mod_clust)
Plot_data <- table(Model_data$Vaccinate_multinom, Model_data$stab_mod_clust) %>% as.data.frame
colnames(Plot_data) <- c("Scenario", "Cluster", "Freq")
Cluster_total <- aggregate(Freq ~ Cluster, data = Plot_data, FUN = sum)
colnames(Cluster_total) <- c("Cluster", "Total")
Plot_data <- merge(Plot_data, Cluster_total, by = "Cluster")
Plot_data$Prop <- Plot_data$Freq / Plot_data$Total
levels(Plot_data$Scenario) <- list("1 - 2" = "S1-2", "3 - 5" = "S3-5", "6 - Never" = "S6-never") 
#png(filename = "Five_covariate_clusters_by_vaccinate_scenario.png", res=600, width=3000, height=2200)
ggplot(data = Plot_data, aes(x = Scenario, y = Prop, fill = Cluster)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(y = "Proportion of farmers in each cluster", x = "Stage of epidemic vaccinated") +
  theme_bw()
#dev.off()
# Save dataframe
#write.csv(Plot_data, "Cluster_prop_by_vaccination_5_stable_terms.csv")

## Plots of variables used in clustering by cluster
# On same plot
Plot_data <- NULL
for(i in colnames(Data[,cluster_colnames])){
  for(j in 1:nlevels(Data$stab_mod_clust)){
    x <- subset(Data, stab_mod_clust==j)
    n <- nrow(x)
    mean <- mean(x[,i])
    sd <- sd(x[,i])
    margin <- qt(0.975, df = n - 1)* sd / sqrt(n)
    LCI <- mean - margin
    UCI <- mean + margin
    Plot_data <- rbind(Plot_data, c(i, j, mean, LCI, UCI))
    colnames(Plot_data) <- c("Covariate", "Cluster", "Mean", "LCI", "UCI")
  }
}
Plot_data <- as.data.frame(Plot_data)
Plot_data[,c("Covariate", "Cluster")] <- lapply(Plot_data[,c("Covariate", "Cluster")], as.factor)
Plot_data[,c("Mean", "LCI", "UCI")] <- lapply(Plot_data[,c("Mean", "LCI", "UCI")], as.numeric)
levels(Plot_data$Covariate) <- list("Trust in\nGovernmental\njudgements\nfor disease\ncontrol" = "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle_norm",
                                    "Physical\nopportunity" = "Physical_opportunity_norm",
                                    "Trust in\nquality of\nadvice from\nveterinarians" = "Farmers.receive.high.quality.veterinary.advice.from.the.veterinary.profession_norm",
                                    "Trust in other\nfarmers to\ncontrol disease" = "I.trust.other.farmers.nationally.to.be.controlling.infectious.diseases.in.their.herds_norm",
                                    "Herd size at\ntime of\noutbreak" = "June.size_norm")
Plot_data <- Plot_data[complete.cases(Plot_data),]
#png(filename = "Five_covariate_clusters_by_covariates.png", res=600, width=3000, height=2200)
ggplot(data = Plot_data, aes(y = Mean)) +
  geom_point(aes(x = Cluster)) +
  geom_errorbar(aes(ymax = UCI, ymin = LCI, x = Cluster), width = 0.25) +
  facet_grid(. ~ Covariate) +
  labs(y = "Normalised scale") +
  theme_bw()
#dev.off()
# Save dataframe
#write.csv(Plot_data, "Cluster_prop_and_CI_by_5_stable_terms.csv")

## Calculate 95% CIs by vaccination
# All
CI_data <- NULL
for(i in levels(Model_data$stab_mod_clust)){
  for(j in levels(Model_data$Vaccinate_multinom_exc_never)){
    n = nrow(subset(Model_data, stab_mod_clust==i & Vaccinate_multinom_exc_never==j))
    denom = nrow(subset(Data, stab_mod_clust==i))
    prop = round(n / denom, digits = 7)
    LCI = round(prop.test(x = n, n = denom, conf.level=.95, correct=FALSE)$conf.int[1], digits = 2)
    UCI = round(prop.test(x = n, n = denom, conf.level=.95, correct=FALSE)$conf.int[2], digits = 2)
    CI_data <- rbind(CI_data,
                     c(j, i, prop, LCI, UCI))
  }
}
CI_data <- as.data.frame(CI_data)
CI_data[,2:4] <- apply(CI_data[,2:4], 2, as.numeric)
colnames(CI_data) <- c("Vaccinate","Cluster",  "Proportion", "LCI", "UCI")
#write.csv(CI_data, "Cluster_prop_and_CI_by_vaccination_5_stable_terms.csv")

#### Stabilised model with only variables above 90% threshold

## Clustering on psychosocial responses
cluster_colnames <- c("I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle_norm", "Physical_opportunity_norm")

## K-means clustering - don't run
## Use map_dbl to run many models with varying value of k (centers)
set.seed(1612)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = Data[,cluster_colnames], centers = k)
  model$tot.withinss
})
## Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:10,
  tot_withinss = tot_withinss
)
## Plot the elbow plot
#png(filename="Two_covariate_clustering_elbow_plot.png", res=600, width=3000, height=2000)
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:20) +
  labs(y = "Total within-cluster sum of squares") +
  theme_bw()
#dev.off() 

## 4 clusters
set.seed(1343)
k_means_cluster_2_var <- kmeans(x = Data[,cluster_colnames], centers = 4)
# Add cluster to dataset
Data$stab_mod_clust_2 <- k_means_cluster_2_var$cluster %>% as.factor()
Model_data$stab_mod_clust_2 <- k_means_cluster_2_var$cluster %>% as.factor()

#### Plot relationships between clusters & variables
# Plot scenario when vaccination is used & k means cluster
table(Model_data$Vaccinate_multinom, Model_data$stab_mod_clust_2)
Plot_data <- table(Model_data$Vaccinate_multinom, Model_data$stab_mod_clust_2) %>% as.data.frame
colnames(Plot_data) <- c("Scenario", "Cluster", "Freq")
Cluster_total <- aggregate(Freq ~ Cluster, data = Plot_data, FUN = sum)
colnames(Cluster_total) <- c("Cluster", "Total")
Plot_data <- merge(Plot_data, Cluster_total, by = "Cluster")
Plot_data$Prop <- Plot_data$Freq / Plot_data$Total
levels(Plot_data$Scenario) <- list("1 - 2" = "S1-2", "3 - 5" = "S3-5", "6 - Never" = "S6-never") 
#png(filename = "Two_covariate_clusters_by_vaccinate_scenario.png", res=600, width=3000, height=2200)
ggplot(data = Plot_data, aes(x = Scenario, y = Prop, fill = Cluster)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733")) +
  labs(y = "Proportion of farmers in each cluster", x = "Stage of epidemic vaccinated") +
  theme_bw()
#dev.off()
# Save dataframe
#write.csv(Plot_data, "Cluster_prop_by_vaccination_2_stable_terms.csv")

# Plots of variables used in clustering by cluster
# On same plot
Plot_data <- NULL
for(i in colnames(Data[,cluster_colnames])){
  for(j in 1:nlevels(Data$stab_mod_clust_2)){
    x <- subset(Data, stab_mod_clust_2==j)
    n <- nrow(x)
    mean <- mean(x[,i])
    sd <- sd(x[,i])
    margin <- qt(0.975, df = n - 1)* sd / sqrt(n)
    LCI <- mean - margin
    UCI <- mean + margin
    Plot_data <- rbind(Plot_data, c(i, j, mean, LCI, UCI))
    colnames(Plot_data) <- c("Covariate", "Cluster", "Mean", "LCI", "UCI")
  }
}
Plot_data <- as.data.frame(Plot_data)
Plot_data[,c("Covariate", "Cluster")] <- lapply(Plot_data[,c("Covariate", "Cluster")], as.factor)
Plot_data[,c("Mean", "LCI", "UCI")] <- lapply(Plot_data[,c("Mean", "LCI", "UCI")], as.numeric)
levels(Plot_data$Covariate) <- list("Trust in Governmental\njudgements for disease control" = "I.trust.governmental.judgements.about.how.to.control.infectious.diseases.in.cattle_norm",
                                    "Physical opportunity" = "Physical_opportunity_norm")
Plot_data <- Plot_data[complete.cases(Plot_data),]
#png(filename = "Two_covariate_clusters_by_covariates.png", res=600, width=3000, height=2200)
ggplot(data = Plot_data, aes(y = Mean)) +
  geom_point(aes(x = Cluster)) +
  geom_errorbar(aes(ymax = UCI, ymin = LCI, x = Cluster), width = 0.25) +
  facet_grid(. ~ Covariate) +
  labs(y = "Normalised scale") +
  theme_bw()
#dev.off()
# Save dataframe
#write.csv(Plot_data, "Cluster_prop_and_CI_by_2_stable_terms.csv")

## Calculate 95% CIs by vaccination
# All
CI_data <- NULL
for(i in levels(Model_data$stab_mod_clust_2)){
  for(j in levels(Model_data$Vaccinate_multinom_exc_never)){
    n = nrow(subset(Model_data, stab_mod_clust_2==i & Vaccinate_multinom_exc_never==j))
    denom = nrow(subset(Data, stab_mod_clust_2==i))
    prop = round(n / denom, digits = 7)
    LCI = round(prop.test(x = n, n = denom, conf.level=.95, correct=FALSE)$conf.int[1], digits = 2)
    UCI = round(prop.test(x = n, n = denom, conf.level=.95, correct=FALSE)$conf.int[2], digits = 2)
    CI_data <- rbind(CI_data,
                     c(j, i, prop, LCI, UCI))
  }
}
CI_data <- as.data.frame(CI_data)
CI_data[,2:4] <- apply(CI_data[,2:4], 2, as.numeric)
colnames(CI_data) <- c("Vaccinate","Cluster",  "Proportion", "LCI", "UCI")
#write.csv(CI_data, "Cluster_prop_and_CI_by_vaccination_2_stable_terms.csv")