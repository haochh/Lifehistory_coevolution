#----------------------------------------------------------------------------#
#                                                                            #
#                   Zhang et al. Code for main analyses                      #
#                   regarding: pairwise PGLS, pPCA, PPA                      #
#                                                                            #
#----------------------------------------------------------------------------#

## load packages 
library(dplyr)
library(purrr)
library(ape)
library(phytools)
library(phylolm)
library(phylopath)
library(FactoMineR)


#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                Mammals                                #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#


#------------------------------------------------------------------------------------------------#
#     1 Preparation: Generate body-size-corrected data for further phylogenetic analyses         #
#------------------------------------------------------------------------------------------------#

### 1.1 Load data & tree
mammal_data <- read.csv("all_mammal_data.csv")
mammal_tree <- read.tree("last_mammal_tree.tre")


### 1.2 Log-transform key life-history traits
# Keep species names
rownames(mammal_data) <- mammal_data$Species
mammal_order <- mammal_data[, c("Species", "Order")]
mammal_trait_data <- mammal_data[, c("BS", "OS", "DT", "AL", "CS", "CF")]
mammal_trait_data <- log(mammal_trait_data)


### 1.3 Body-size correction (using BS) with PGLS
mammal_residual_data <- data.frame(row.names = rownames(mammal_trait_data))
for (trait in names(mammal_trait_data)) {  
  if (trait != "BS") {  
    model <- phylolm(formula = as.formula(paste(trait, "~ BS")), phy = mammal_tree, data = mammal_trait_data, model = "lambda", boot = 1000)  
    mammal_residual_data[[trait]] <- residuals(model)  
  }
}
### 1.4 Standardize data
mammal_data_corrected <- scale(mammal_residual_data)
mammal_data_corrected <- as.data.frame(mammal_data_corrected)

### 1.5 Get the the body-size-corrected data used for phylogenetic analyses
write.csv(mammal_data_corrected, "mammal_data_corrected.csv", row.names = TRUE)

#------------------------------------------------------------------------------------------------#
#             2  Modelling trait evolution to determine the best model for evolution             #
#------------------------------------------------------------------------------------------------#

### 2.1 Functions to fit and compare trait evolution models
fit_models_for_trait <- function(trait, tree, data) {
  
  form <- reformulate("1", response = trait)
  
  models <- c(
    "BM",
    "EB",
    "lambda",
    "OUrandomRoot",
    "OUfixedRoot"
  )
  
  fits <- lapply(models, function(m) {
    
    phylolm(
      formula = form,
      phy = tree,
      data = data,
      model = m
    )
    
  })
  
  names(fits) <- models
  
  tibble(
    trait = trait,
    model = models,
    AIC = sapply(fits, AIC)
  ) %>%
    mutate(
      deltaAIC = AIC - min(AIC)
    ) %>%
    arrange(AIC)
}

### 2.2 Fit models for 5 body-size-corrected traits
traits_to_fit <- c("OS", "DT", "AL", "CS", "CF")

mammal_model_comp <- map_dfr(
  traits_to_fit,
  function(x) {
    fit_models_for_trait(
      trait = x,
      tree = mammal_tree,
      data = mammal_data_corrected
    )
  }
)

### 2.3 Model comparison based on AIC, ΔAIC were computed. OU model turns out to be the best model for most traits
mammal_model_comp


#--------------------------------------------------------------------------------#
#             3  Pairwise PGLS for each trait pair based on OU model             #
#--------------------------------------------------------------------------------#
rownames(mammal_data_corrected) <- rownames(mammal_residual_data)

### 3.1 Fit PGLS for all trait pairs
mammal_pgls1 <- phylolm(formula = DT ~ OS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls2 <- phylolm(formula = CS ~ OS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls3 <- phylolm(formula = CF ~ OS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls4 <- phylolm(formula = AL ~ OS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)
 
mammal_pgls5 <- phylolm(formula = DT ~ CS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls6 <- phylolm(formula = CF ~ CS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls7 <- phylolm(formula = AL ~ CS, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls8 <- phylolm(formula = AL ~ DT, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls9 <- phylolm(formula = CF ~ DT, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)

mammal_pgls10 <- phylolm(formula = AL ~ CF, phy = mammal_tree, data = mammal_data_corrected, model = "OUrandomRoot", boot = 1000)


### 3.2 Summarize PGLS results
save(mammal_pgls1, mammal_pgls2, mammal_pgls3, mammal_pgls4, mammal_pgls5,
     mammal_pgls6, mammal_pgls7, mammal_pgls8, mammal_pgls9, mammal_pgls10,
  file = "mammal_PGLS_results.RData"
)

#-----------------------------------------------------------#
#             4  pPCA analysis; model = lambda              #
#-----------------------------------------------------------#
### 4.1 Prepare the data
mammal_ppca_matrix <- as.matrix(mammal_data_corrected)

# Match species order with phylogenetic tree
mammal_ppca_matrix <- mammal_ppca_matrix[mammal_tree$tip.label,, drop = FALSE]

### 4.2 Phylogenetic PCA
pPCA_mammal <- phyl.pca(mammal_tree, mammal_ppca_matrix, method = "lambda", mode = "cov")

### 4.3 Save pPCA results
save(pPCA_mammal, file = "mammal_pPCA_results.RData")



#--------------------------------------------------------------------------#
#             5  Phylogenetic path analysis (PPA); model = OU              #
#--------------------------------------------------------------------------#

### 5.1 Define candidate causal models
models <- define_model_set(
  Model_0=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT),
  Model_1=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS),
  Model_2=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS),
  Model_3=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS),
  Model_4=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CF~CS),
  Model_5=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF),
  Model_6=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS),
  Model_7=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS),
  Model_8=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF),
  Model_9=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,CF~CS),
  Model_10=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS),
  Model_11=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF),
  Model_12=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS),
  Model_13=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,AL~CF),
  Model_14=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,CF~CS),
  Model_15=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF,CF~CS),
  Model_16=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS),
  Model_17=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF),
  Model_18=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,CF~CS),
  Model_19=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF),
  Model_20=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,CF~CS),
  Model_21=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF,CF~CS),
  Model_22=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF),
  Model_23=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,CF~CS),
  Model_24=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF,CF~CS),
  Model_25=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS,AL~CS),
  Model_26=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,AL~CF),
  Model_27=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,CF~CS),
  Model_28=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF,CF~CS),
  Model_29=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF,CF~CS),
  Model_30=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS),
  Model_31=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS)
  
)

### 5.2 Run the PPA analysis
PPA_mammals <- phylo_path(models, data = mammal_data_corrected, tree = mammal_tree, model = 'OUrandomRoot')

### 5.3 Save PPA results
save(PPA_mammals,file = "mammal_PPA_results.RData")

#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                  Birds                                #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#

#------------------------------------------------------------------------------------------------#
#     1 Preparation: Generate body-size-corrected data for further phylogenetic analyses         #
#------------------------------------------------------------------------------------------------#
### 1.1 Load data & tree
bird_data <- read.csv("all_bird_data.csv")
bird_tree <- read.tree("last_bird_tree.tre")


### 1.2 Log-transform key life-history traits
# Keep species names
rownames(bird_data) <- bird_data$Species
bird_order <- bird_data[, c("Species", "Order")]
bird_trait_data <- bird_data[, c("BS", "OS", "DT", "AL", "CS", "CF")]
bird_trait_data <- log(bird_trait_data)


### 1.3 Body-size correction (using BS) with PGLS
bird_residual_data <- data.frame(row.names = rownames(bird_trait_data))
for (trait in names(bird_trait_data)) {  
  if (trait != "BS") {  
    model <- phylolm(formula = as.formula(paste(trait, "~ BS")), phy = bird_tree, data = bird_trait_data, model = "lambda", boot = 1000)  
    bird_residual_data[[trait]] <- residuals(model)  
  }
}
### 1.4 Standardize data
bird_data_corrected <- scale(bird_residual_data)
bird_data_corrected <- as.data.frame(bird_data_corrected)
### 1.5 Get the the body-size-corrected data used for phylogenetic analyses
write.csv(bird_data_corrected, "bird_data_corrected.csv", row.names = TRUE)

#------------------------------------------------------------------------------------------------#
#             2  Modelling trait evolution to determine the best model for evolution             #
#------------------------------------------------------------------------------------------------#
### 2.1 Functions to fit and compare trait evolution models
fit_models_for_trait <- function(trait, tree, data) {
  
  form <- reformulate("1", response = trait)
  
  models <- c(
    "BM",
    "EB",
    "lambda",
    "OUrandomRoot",
    "OUfixedRoot"
  )
  
  fits <- lapply(models, function(m) {
    
    phylolm(
      formula = form,
      phy = tree,
      data = data,
      model = m
    )
    
  })
  
  names(fits) <- models
  
  tibble(
    trait = trait,
    model = models,
    AIC = sapply(fits, AIC)
  ) %>%
    mutate(
      deltaAIC = AIC - min(AIC)
    ) %>%
    arrange(AIC)
}

### 2.2 Fit models for 5 body-size-corrected traits
traits_to_fit <- c("OS", "DT", "AL", "CS", "CF")

bird_model_comp <- map_dfr(
  traits_to_fit,
  function(x) {
    fit_models_for_trait(
      trait = x,
      tree = bird_tree,
      data = bird_data_corrected
    )
  }
)
### 2.3 Model comparison based on AIC, ΔAIC were computed. 
bird_model_comp

#--------------------------------------------------------------------------------#
#             3  Pairwise PGLS for each trait pair based on OU model             #
#--------------------------------------------------------------------------------#
rownames(bird_data_corrected) <- rownames(bird_residual_data)
### 3.1 Fit PGLS for all trait pairs
bird_pgls1 <- phylolm(formula = DT ~ OS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls2 <- phylolm(formula = CS ~ OS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls3 <- phylolm(formula = CF ~ OS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls4 <- phylolm(formula = AL ~ OS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls5 <- phylolm(formula = DT ~ CS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls6 <- phylolm(formula = CF ~ CS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls7 <- phylolm(formula = AL ~ CS, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls8 <- phylolm(formula = AL ~ DT, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls9 <- phylolm(formula = CF ~ DT, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)

bird_pgls10 <- phylolm(formula = AL ~ CF, phy = bird_tree, data = bird_data_corrected, model = "OUrandomRoot", boot = 1000)


### 3.2 Summarize PGLS results
save(bird_pgls1, bird_pgls2, bird_pgls3, bird_pgls4, bird_pgls5,
     bird_pgls6, bird_pgls7, bird_pgls8, bird_pgls9, bird_pgls10,
     file = "bird_PGLS_results.RData"
)


#-----------------------------------------------------------#
#             4  pPCA analysis; model = lambda              #
#-----------------------------------------------------------#
### 4.1 Prepare tha data
bird_ppca_matrix <- as.matrix(bird_data_corrected)

bird_ppca_matrix <- bird_ppca_matrix[bird_tree$tip.label,, drop = FALSE]

### 4.2 Phylogenetic PCA
pPCA_bird <- phyl.pca(bird_tree, bird_ppca_matrix, method = "lambda", mode = "cov")

### 4.3 Save pPCA results
save(pPCA_bird, file = "bird_pPCA_results.RData")


#--------------------------------------------------------------------------#
#             5  Phylogenetic path analysis (PPA); model = OU              #
#--------------------------------------------------------------------------#
### 5.1 Define candidate causal models
models <- define_model_set(
  Model_0=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT),
  Model_1=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS),
  Model_2=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS),
  Model_3=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS),
  Model_4=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CF~CS),
  Model_5=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF),
  Model_6=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS),
  Model_7=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS),
  Model_8=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF),
  Model_9=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,CF~CS),
  Model_10=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS),
  Model_11=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF),
  Model_12=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS),
  Model_13=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,AL~CF),
  Model_14=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,CF~CS),
  Model_15=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF,CF~CS),
  Model_16=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS),
  Model_17=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF),
  Model_18=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,CF~CS),
  Model_19=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF),
  Model_20=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,CF~CS),
  Model_21=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF,CF~CS),
  Model_22=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF),
  Model_23=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,CF~CS),
  Model_24=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF,CF~CS),
  Model_25=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS,AL~CS),
  Model_26=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,AL~CF),
  Model_27=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,CF~CS),
  Model_28=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF,CF~CS),
  Model_29=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF,CF~CS),
  Model_30=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS),
  Model_31=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS)
  
)

### 5.2 Run the PPA analysis
PPA_birds <- phylo_path(models, data = bird_data_corrected, tree = bird_tree, model = 'OUrandomRoot')

### 5.3 Save PPA results
save(PPA_birds,file = "bird_PPA_results.RData")



#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                  Reptiles                             #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#

#------------------------------------------------------------------------------------------------#
#     1 Preparation: Generate body-size-corrected data for further phylogenetic analyses         #
#------------------------------------------------------------------------------------------------#
### 1.1 Load data & tree
reptile_data <- read.csv("all_reptile_data.csv")
reptile_tree <- read.tree("last_reptile_tree.tre")


### 1.2 Log-transform key life-history traits
# Keep species names
rownames(reptile_data) <- reptile_data$Species
reptile_family <- reptile_data[, c("Species", "Family")]
reptile_trait_data <- reptile_data[, c("BS", "OS", "DT", "AL", "CS", "CF")]
reptile_trait_data <- log(reptile_trait_data)


### 1.3 Body-size correction (using BS) with PGLS
reptile_residual_data <- data.frame(row.names = rownames(reptile_trait_data))
for (trait in names(reptile_trait_data)) {  
  if (trait != "BS") {  
    model <- phylolm(formula = as.formula(paste(trait, "~ BS")), phy = reptile_tree, data = reptile_trait_data, model = "lambda", boot = 1000)  
    reptile_residual_data[[trait]] <- residuals(model)  
  }
}
### 1.4 Standardize data
reptile_data_corrected <- scale(reptile_residual_data)
reptile_data_corrected <- as.data.frame(reptile_data_corrected)

### 1.5 Get the the body-size-corrected data used for phylogenetic analyses
write.csv(reptile_data_corrected, "reptile_data_corrected.csv", row.names = TRUE)

#------------------------------------------------------------------------------------------------#
#             2  Modelling trait evolution to determine the best model for evolution             #
#------------------------------------------------------------------------------------------------#
### 2.1 Functions to fit and compare trait evolution models
fit_models_for_trait <- function(trait, tree, data) {
  
  form <- reformulate("1", response = trait)
  
  models <- c(
    "BM",
    "EB",
    "lambda",
    "OUrandomRoot",
    "OUfixedRoot"
  )
  
  fits <- lapply(models, function(m) {
    
    phylolm(
      formula = form,
      phy = tree,
      data = data,
      model = m
    )
    
  })
  
  names(fits) <- models
  
  tibble(
    trait = trait,
    model = models,
    AIC = sapply(fits, AIC)
  ) %>%
    mutate(
      deltaAIC = AIC - min(AIC)
    ) %>%
    arrange(AIC)
}

### 2.2 Fit models for 5 body-size-corrected traits
traits_to_fit <- c("OS", "DT", "AL", "CS", "CF")

reptile_model_comp <- map_dfr(
  traits_to_fit,
  function(x) {
    fit_models_for_trait(
      trait = x,
      tree = reptile_tree,
      data = reptile_data_corrected
    )
  }
)

### 2.3 Model comparison based on AIC, ΔAIC were computed.
reptile_model_comp


#--------------------------------------------------------------------------------#
#             3  Pairwise PGLS for each trait pair based on OU model             #
#--------------------------------------------------------------------------------#
rownames(reptile_data_corrected) <- rownames(reptile_residual_data)
### 3.1 Fit PGLS for all trait pairs
reptile_pgls1 <- phylolm(formula = DT ~ OS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls2 <- phylolm(formula = CS ~ OS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls3 <- phylolm(formula = CF ~ OS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls4 <- phylolm(formula = AL ~ OS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls5 <- phylolm(formula = DT ~ CS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls6 <- phylolm(formula = CF ~ CS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls7 <- phylolm(formula = AL ~ CS, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls8 <- phylolm(formula = AL ~ DT, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls9 <- phylolm(formula = CF ~ DT, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)

reptile_pgls10 <- phylolm(formula = AL ~ CF, phy = reptile_tree, data = reptile_data_corrected, model = "OUrandomRoot", boot = 1000)


### 3.2 Summarize PGLS results
save(reptile_pgls1, reptile_pgls2, reptile_pgls3, reptile_pgls4, reptile_pgls5,
     reptile_pgls6, reptile_pgls7, reptile_pgls8, reptile_pgls9, reptile_pgls10,
     file = "reptile_PGLS_results.RData"
)


#-----------------------------------------------------------#
#             4  pPCA analysis; model = lambda              #
#-----------------------------------------------------------#
### 4.1 Prepare tha data
reptile_ppca_matrix <- as.matrix(reptile_data_corrected)

# Match species order with phylogenetic tree
reptile_ppca_matrix <- reptile_ppca_matrix[reptile_tree$tip.label,, drop = FALSE]

### 4.2 Phylogenetic PCA
pPCA_reptile <- phyl.pca(reptile_tree, reptile_ppca_matrix, method = "lambda", mode = "cov")

### 4.3 Save pPCA results
save(pPCA_reptile, file = "reptile_pPCA_results.RData")


#--------------------------------------------------------------------------#
#             5  Phylogenetic path analysis (PPA); model = OU              #
#--------------------------------------------------------------------------#
### 5.1 Define candidate causal models
models <- define_model_set(
  Model_0=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT),
  Model_1=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS),
  Model_2=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS),
  Model_3=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS),
  Model_4=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CF~CS),
  Model_5=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF),
  Model_6=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS),
  Model_7=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS),
  Model_8=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF),
  Model_9=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,CF~CS),
  Model_10=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS),
  Model_11=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF),
  Model_12=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS),
  Model_13=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,AL~CF),
  Model_14=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CS,CF~CS),
  Model_15=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~CF,CF~CS),
  Model_16=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS),
  Model_17=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF),
  Model_18=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,CF~CS),
  Model_19=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF),
  Model_20=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,CF~CS),
  Model_21=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CF,CF~CS),
  Model_22=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF),
  Model_23=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,CF~CS),
  Model_24=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CF,CF~CS),
  Model_25=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,CF~CS,AL~CS),
  Model_26=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,AL~CF),
  Model_27=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CS,CF~CS),
  Model_28=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~OS,AL~CF,CF~CS),
  Model_29=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,CS~OS,AL~CS,AL~CF,CF~CS),
  Model_30=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS),
  Model_31=c(CF~OS,DT~OS,DT~CS,AL~DT,CF~DT,AL~OS,AL~CS,AL~CF,CF~CS)
  
)

### 5.2 Run the PPA analysis
PPA_reptiles <- phylo_path(models, data = reptile_data_corrected, tree = reptile_tree, model = 'OUrandomRoot')

### 5.3 Save PPA results
save(PPA_reptiles,file = "reptile_PPA_results.RData")


