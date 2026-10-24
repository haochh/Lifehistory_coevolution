#----------------------------------------------------------------------------#
#                                                                            #
#                   Zhang et al. Code for main analyses                      #
#                   regarding: multivariate OU modelling                     #
#                                                                            #
#----------------------------------------------------------------------------#

library(tidyverse)
library(ape)
library(mvSLOUCH)
library(phylopath)
library(PCMBaseCpp) 
library(igraph)
library(tidyverse)
library(foreach)
library(doParallel)
library(ggpubr)
library(patchwork)


#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                Mammals                                #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#

#-------------------------------------
#      Step 1: Data preparation 
#-------------------------------------

### 1 Scale the maximum tree height to 1
mammal_tree <- read.tree("last_mammal_tree.tre")
tree_height <-mvSLOUCH::phyltree_paths(mammal_tree)$tree_height
mammal_tree_scaled <- mammal_tree
mammal_tree_scaled$edge.length <- mammal_tree_scaled$edge.length/tree_height
# Check the tree height = 1
mvSLOUCH::phyltree_paths(mammal_tree_scaled)$tree_height

### 2 Prepare the data to have the order of columns constrained
mammal_data <- read.csv("mammal_data_corrected.csv")
mammal_data <- mammal_data[,1:6]
mammal_data <- mammal_data %>% column_to_rownames(var = "X")
colnames(mammal_data) <- c("OS", "DT", "AL", "CS", "CF")
mammal_data <- mammal_data[, c("DT", "OS", "AL", "CS", "CF")]
mammal_data_ordered <- mammal_data[mammal_tree_scaled$tip.label, ]
mammal_data_ordered <- as.matrix(mammal_data_ordered)


#--------------------------------------------------------
#      Step 2: Covert all PPA matrix to A matrix
#--------------------------------------------------------
# Our mvOU modelling was based on the same direct/indirect structures used in PPA. 
# We thus convert the PPA matrix for each candidate model to its corresponding A matrix. See Materials and Methods for details.

### 1 Load PPA results for mammal
load("mammal_PPA_results.RData")
length(PPA_mammals[["model_set"]]) # = 32

### 2 Convert PPA matrix to A matrix
A.Matrix <- list()
for (i in 1:32) {
  # Get the name of the candidate model
  model_name <- names(PPA_mammals[["model_set"]][i])
  # Get the PPA matrix
  model <- PPA_mammals[["model_set"]][i][[1]]
  # Re-order the matrix
  model_reordered <- model[colnames(mammal_data_ordered), colnames(mammal_data_ordered)]
  # Transpose the matrix
  model_reordered_t <- t(model_reordered)
  # Conduct the transformation
  model_reordered_t[model_reordered_t == 1] <- NA
  diag(model_reordered_t) <- "+"
  mode(model_reordered_t) <- "character"
  # Get the final A matrix
  A <- model_reordered_t
  # Store them in the list
  A.Matrix[[model_name]] <- A
}

#---------------------------------------------------------------------------------------------------------
#        Step 3: Conducting mvOU modelling to fit mvOU models with varying A matrix; Syy = Diagonal
#---------------------------------------------------------------------------------------------------------

### 1 Initialize the result dataframe
Mammal_SyyDiagonal_OUOU_fit_niter500 <- data.frame()
niter <- 500  

# Set up parallel backend
detectCores()
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("mammal_tree_scaled", "mammal_data_ordered", "A.Matrix"))

Mammal_SyyDiagonal_OUOU_fit_niter500 <- foreach(
  n = 1:niter, 
  .combine = rbind,
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(mammal_tree_scaled, mammal_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", 
                               Syytype="Diagonal",  
                               diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- data.frame()
  for (i in 1:length(A.Matrix)) {
    # Get the A matrix
    model_A <- A.Matrix[[i]]
    # Get the model name
    nameA <- names(A.Matrix)[i]
    # Fit the model
    OUOUmodel <- ouchModel(mammal_tree_scaled, mammal_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", # or "UpperTri" 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    
    # Get the result
    result <- data.frame(
      Model = nameA,
      AIC = OUOUmodel$FinalFound$ParamSummary$aic,
      AICc = OUOUmodel$FinalFound$ParamSummary$aic.c,
      SIC = OUOUmodel$FinalFound$ParamSummary$sic,
      BIC = OUOUmodel$FinalFound$ParamSummary$bic,
      R2 = OUOUmodel$FinalFound$ParamSummary$RSS$R2,
      LogLik = OUOUmodel$FinalFound$ParamSummary$LogLik,
      DOF = OUOUmodel$FinalFound$ParamSummary$dof,
      Group = "Mammal",
      Syytype="Diagonal",
      Iteration = n,
      Seed1 = seed1,
      Seed2 = seed2
    )
    iteration_result <- rbind(iteration_result, result)
  }
  return(iteration_result)
}


stopCluster(cl)
### 2 Save the result dataframe
save(Mammal_SyyDiagonal_OUOU_fit_niter500, "Mammal_SyyDiagonal_OUOU_fit_niter500.RData")

#---------------------------------------------------------------------------------------------------------------------------------
#        Step 4: Model diagnosis: Final the best solution across various candidate models and staring points based on AICc 
#---------------------------------------------------------------------------------------------------------------------------------
### 1 Load the pre-computed mvOU outputs from the Step 3
load("Mammal_SyyDiagonal_OUOU_fit_niter500.RData")
length(unique(Mammal_SyyDiagonal_OUOU_fit_niter500$Model)) # number of models = 32

### 2 Get the best model: Model == "Model_6" & Iteration == 56; AICc = 10551.90
Mammal_SyyDiagonal_OUOU_fit_niter500 %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)


#---------------------------------------------------------------------------------------------------------
#         Step 5: Generate all possible combinations of directional interactions of the best model
#---------------------------------------------------------------------------------------------------------
### 1 The best model is Model_6
A.Matrix[["Model_6"]]

### 2 get all possible combinations of directional interactions of the best model, n = 2^7 = 128
mat <- matrix(
  c("+", "NA", "0", "NA", "0",
    "0", "+", "0", "0", "0",
    "NA", "NA", "+", "0", "0",
    "0", "NA", "0", "+", "0",
    "NA", "NA", "0", "0", "+"),
  nrow = 5, byrow = TRUE,
  dimnames = list(c("DT", "OS", "AL", "CS", "CF"), 
                  c("DT", "OS", "AL", "CS", "CF"))
)

n <- nrow(mat)
non_diag <- matrix(TRUE, n, n)
diag(non_diag) <- FALSE
arrow_positions <- which(mat == "NA" & non_diag, arr.ind = TRUE)
k <- nrow(arrow_positions)
combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))

all_models <- vector("list", nrow(combos))

for (i in 1:nrow(combos)) {
  new_mat <- mat  
  for (j in 1:k) {
    if (combos[i, j]) {  
      row_idx <- arrow_positions[j, "row"]
      col_idx <- arrow_positions[j, "col"]
      new_mat[row_idx, col_idx] <- "0"         
      new_mat[col_idx, row_idx] <- "NA"        
    }
  }
  
  all_models[[i]] <- new_mat
}

### 3 Name the model list
Model_6_list <- all_models

#--------------------------------------------------------------------------------------------------------------------------
#   Step 6:  Constrain the A matrix of model with different directions to reflect the sign of the PPA-OU best model
#---------------------------------------------------------------------------------------------------------------------------

# When there is a positive parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '+'
# When there is a negative parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '-'

### 1 Get the PPA model 
best_model <- best(PPA_mammals)
PPA_OU_mammal <- best_model[["coef"]]

# Change the order of columns and rows to match with A matrix
PPA_OU_mammal  <- PPA_OU_mammal[, c("DT", "OS", "AL", "CS", "CF")]
PPA_OU_mammal  <- PPA_OU_mammal[c("DT", "OS", "AL", "CS", "CF"),]


### 2 Initialize a new list for modified matrices
constrained_Model_6_list <- vector("list", length(Model_6_list))

# Get reference matrix dimensions
ref_rows <- rownames(PPA_OU_mammal)
ref_cols <- colnames(PPA_OU_mammal)

for (idx in seq_along(Model_6_list)) {
  mat <- Model_6_list[[idx]]
  new_mat <- mat  # Create modifiable copy
  
  # Get current matrix's row/column names
  curr_rows <- rownames(mat)
  curr_cols <- colnames(mat)
  
  # Find "NA" positions with names
  na_positions <- which(mat == "NA", arr.ind = TRUE)
  
  for (pos_idx in 1:nrow(na_positions)) {
    i <- na_positions[pos_idx, "row"]
    j <- na_positions[pos_idx, "col"]
    
    row_name <- curr_rows[i]
    col_name <- curr_cols[j]
    
    # Check same position in reference
    ref_val_same <- PPA_OU_mammal[row_name, col_name]
    
    if (ref_val_same > 0) {
      new_mat[i, j] <- "+"
    } else if (ref_val_same < 0) {
      new_mat[i, j] <- "-"
    } else {
      # Check transposed position
      ref_val_trans <- PPA_OU_mammal[col_name, row_name]
      
      if (ref_val_trans > 0) {
        new_mat[i, j] <- "+"
      } else if (ref_val_trans < 0) {
        new_mat[i, j] <- "-"
      }
      # If both zero, "NA" remains unchanged
    }
  }
  
  constrained_Model_6_list[[idx]] <- new_mat
}

#-----------------------------------------------------------------------------------------------------------
#   Step 7:  Fit all directional mvOU models to get the best model reflect the direction of coevolution 
#-----------------------------------------------------------------------------------------------------------

### 1 Prepare the data
mammal_data_ordered <- as.matrix(mammal_data_ordered)

# Rename the candidate model
constrained_Model_6_list_mammal_OU <- constrained_Model_6_list

### 2 Setup for parallel work
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("mammal_tree_scaled", "mammal_data_ordered", "constrained_Model_6_list_mammal_OU"))


### 3 Run the parallel iterations of model fitting: n = 100*128
niter <- 100

result_constrained_Model_6_mammal_OU_niter100 <- foreach(
  n = 1:niter, 
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(mammal_tree_scaled, mammal_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", # Assuming A type as DecomposablePositive
                               Syytype="Diagonal", diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- list()
  for (i in 1:length(constrained_Model_6_list_mammal_OU)) {
    # Get the A matrix
    model_A <- constrained_Model_6_list_mammal_OU[[i]]
    model_A[model_A == "NA"] <- NA
    
    # Fit the model
    OUOUmodel <- ouchModel(mammal_tree_scaled, mammal_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    iteration_result[[i]] <- OUOUmodel
  }
  return(iteration_result)
}

stopCluster(cl)

### 4 Save the results
save(result_constrained_Model_6_mammal_OU_niter100, file = "result_constrained_Model_6_mammal_OU_niter100.RData")

#----------------------------------------------------------------------------------------------
#   Step 8:  Select the best model which deciphers the directional coevolution among traits  
#----------------------------------------------------------------------------------------------
load("result_constrained_Model_6_mammal_OU_niter100.RData")

### 1 Get the summary of AICc for each run
causal_constrained_mammal_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:128) {
    df <- data.frame(
      AICc = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_6_mammal_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Mammal"
    )
    causal_constrained_mammal_OU_dignosis <- rbind(causal_constrained_mammal_OU_dignosis, df)
  }
}
### 2 The best model is when Iteration = 5, Model = 55
causal_constrained_mammal_OU_dignosis %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

causal_constrained_mammal_OU_best <- result_constrained_Model_6_mammal_OU_niter100[[5]][[55]]

### 3 Get the best A matrix
A.causal_constrained_mammal_OU_best <- causal_constrained_mammal_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]


### 4 Get the eigen-decomposition results of the A matrix
half.t.mammal_OU <- rbind(causal_constrained_mammal_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                          causal_constrained_mammal_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])



#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                  Birds                                #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#

#-------------------------------------
#      Step 1: Data preparation 
#-------------------------------------
### 1 Scale the maximum tree height to 1
bird_tree <- read.tree("last_bird_tree.tre")
tree_height <- mvSLOUCH::phyltree_paths(bird_tree)$tree_height
bird_tree_scaled <- bird_tree
bird_tree_scaled$edge.length <- bird_tree_scaled$edge.length/tree_height

# Check the tree height = 1
mvSLOUCH::phyltree_paths(bird_tree_scaled)$tree_height

### 2 Prepare the data to have the order of columns constrained
bird_data <- read.csv("bird_data_corrected.csv")
bird_data <- bird_data[,1:6]
bird_data <- bird_data %>% column_to_rownames(var = "X")
colnames(bird_data) <- c("OS", "DT", "AL", "CS", "CF")
bird_data <- bird_data[, c("DT", "OS", "AL", "CS", "CF")]
bird_data_ordered <- bird_data[bird_tree_scaled$tip.label, ]
bird_data_ordered <- as.matrix(bird_data_ordered)

#--------------------------------------------------------
#      Step 2: Covert all PPA matrix to A matrix
#--------------------------------------------------------
# Our mvOU modelling was based on the same direct/indirect structures used in PPA. 
# We thus convert the PPA matrix for each candidate model to its corresponding A matrix. See Materials and Methods for details.

### 1 Load PPA results for bird
load("bird_PPA_results.RData")
length(PPA_birds[["model_set"]]) # =32

### 2 Convert PPA matrix to A matrix
A.Matrix <- list()
for (i in 1:32) {
  # Get the name of the candidate model
  model_name <- names(PPA_birds[["model_set"]][i])
  # Get the PPA matrix
  model <- PPA_birds[["model_set"]][i][[1]]
  # Re-order the matrix
  model_reordered <- model[colnames(bird_data_ordered), colnames(bird_data_ordered)]
  # Transpose the matrix
  model_reordered_t <- t(model_reordered)
  # Conduct the transformation
  model_reordered_t[model_reordered_t == 1] <- NA
  diag(model_reordered_t) <- "+"
  mode(model_reordered_t) <- "character"
  # Get the final A matrix
  A <- model_reordered_t
  # Store them in the list
  A.Matrix[[model_name]] <- A
}

#---------------------------------------------------------------------------------------------------------
#        Step 3: Conducting mvOU modelling to fit mvOU models with varying A matrix; Syy = Diagonal
#---------------------------------------------------------------------------------------------------------
### 1 Initialize the result dataframe
Bird_SyyDiagonal_OUOU_fit_niter500 <- data.frame()
# Number of iterations
niter <- 500  

# Set up parallel backend
detectCores()
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("bird_tree_scaled", "bird_data_ordered", "A.Matrix"))

Bird_SyyDiagonal_OUOU_fit_niter500 <- foreach(
  n = 1:niter, 
  .combine = rbind,
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(bird_tree_scaled, bird_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", 
                               Syytype="Diagonal",  # or "UpperTri" 
                               diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- data.frame()
  for (i in 1:length(A.Matrix)) {
    # Get the A matrix
    model_A <- A.Matrix[[i]]
    # Get the model name
    nameA <- names(A.Matrix)[i]
    # Fit the model
    OUOUmodel <- ouchModel(bird_tree_scaled, bird_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", # or "UpperTri" 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    
    # Get the result
    result <- data.frame(
      Model = nameA,
      AIC = OUOUmodel$FinalFound$ParamSummary$aic,
      AICc = OUOUmodel$FinalFound$ParamSummary$aic.c,
      SIC = OUOUmodel$FinalFound$ParamSummary$sic,
      BIC = OUOUmodel$FinalFound$ParamSummary$bic,
      R2 = OUOUmodel$FinalFound$ParamSummary$RSS$R2,
      LogLik = OUOUmodel$FinalFound$ParamSummary$LogLik,
      DOF = OUOUmodel$FinalFound$ParamSummary$dof,
      Group = "Bird",
      Syytype="Diagonal",
      Iteration = n,
      Seed1 = seed1,
      Seed2 = seed2
    )
    iteration_result <- rbind(iteration_result, result)
  }
  return(iteration_result)
}


stopCluster(cl)

### 2 Save the result dataframe
save(Bird_SyyDiagonal_OUOU_fit_niter500, "Bird_SyyDiagonal_OUOU_fit_niter500.RData")

#---------------------------------------------------------------------------------------------------------------------------------
#        Step 4: Model diagnosis: Final the best solution across various candidate models and staring points based on AICc 
#---------------------------------------------------------------------------------------------------------------------------------
### 1 Load the pre-computed mvOU outputs from the Step 3
load("Bird_SyyDiagonal_OUOU_fit_niter500.RData")
length(unique(Bird_SyyDiagonal_OUOU_fit_niter500$Model)) # number of models = 32

### 2 Get the best model: Model == "Model_3" & Iteration == 27; AICc = 2198.94
Bird_SyyDiagonal_OUOU_fit_niter500 %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

#---------------------------------------------------------------------------------------------------------
#         Step 5: Generate all possible combinations of directional interactions of the best model
#---------------------------------------------------------------------------------------------------------
### 1 The best model is Model_3
A.Matrix[["Model_3"]]

### 2 get all possible combinations of directional interactions of the best model, n = 2^8 = 64
mat <- matrix(
  c("+", "NA", "0", "NA", "0",
    "0", "+", "0", "0", "0",
    "NA", "0", "+", "NA", "0",
    "0", "0", "0", "+", "0",
    "NA", "NA", "0", "0", "+"),
  nrow = 5, byrow = TRUE,
  dimnames = list(c("DT", "OS", "AL", "CS", "CF"), 
                  c("DT", "OS", "AL", "CS", "CF"))
)

n <- nrow(mat)
non_diag <- matrix(TRUE, n, n)
diag(non_diag) <- FALSE
arrow_positions <- which(mat == "NA" & non_diag, arr.ind = TRUE)
k <- nrow(arrow_positions)
combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))

all_models <- vector("list", nrow(combos))

for (i in 1:nrow(combos)) {
  new_mat <- mat  
  for (j in 1:k) {
    if (combos[i, j]) {  
      row_idx <- arrow_positions[j, "row"]
      col_idx <- arrow_positions[j, "col"]
      new_mat[row_idx, col_idx] <- "0"         
      new_mat[col_idx, row_idx] <- "NA"        
    }
  }
  
  all_models[[i]] <- new_mat
}

### 3 Name the model list
Model_3_list <- all_models

#--------------------------------------------------------------------------------------------------------------------------
#   Step 6:  Constrain the A matrix of model with different directions to reflect the sign of the PPA-OU best model
#---------------------------------------------------------------------------------------------------------------------------

# When there is a positive parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '+'
# When there is a negative parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '-'

### 1 Get the PPA model 
best_model <- best(PPA_birds)
PPA_OU_bird <- best_model[["coef"]]

# Change the order of columns and rows to match with A matrix
PPA_OU_bird  <- PPA_OU_bird[, c("DT", "OS", "AL", "CS", "CF")]
PPA_OU_bird  <- PPA_OU_bird[c("DT", "OS", "AL", "CS", "CF"),]


### 2 Initialize a new list for modified matrices
constrained_Model_3_list <- vector("list", length(Model_3_list))

# Get reference matrix dimensions
ref_rows <- rownames(PPA_OU_bird)
ref_cols <- colnames(PPA_OU_bird)

for (idx in seq_along(Model_3_list)) {
  mat <- Model_3_list[[idx]]
  new_mat <- mat  # Create modifiable copy
  
  # Get current matrix's row/column names
  curr_rows <- rownames(mat)
  curr_cols <- colnames(mat)
  
  # Find "NA" positions with names
  na_positions <- which(mat == "NA", arr.ind = TRUE)
  
  for (pos_idx in 1:nrow(na_positions)) {
    i <- na_positions[pos_idx, "row"]
    j <- na_positions[pos_idx, "col"]
    
    row_name <- curr_rows[i]
    col_name <- curr_cols[j]
    
    # Check same position in reference
    ref_val_same <- PPA_OU_bird[row_name, col_name]
    
    if (ref_val_same > 0) {
      new_mat[i, j] <- "+"
    } else if (ref_val_same < 0) {
      new_mat[i, j] <- "-"
    } else {
      # Check transposed position
      ref_val_trans <- PPA_OU_bird[col_name, row_name]
      
      if (ref_val_trans > 0) {
        new_mat[i, j] <- "+"
      } else if (ref_val_trans < 0) {
        new_mat[i, j] <- "-"
      }
      # If both zero, "NA" remains unchanged
    }
  }
  
  constrained_Model_3_list[[idx]] <- new_mat
}

#-----------------------------------------------------------------------------------------------------------
#   Step 7:  Fit all directional mvOU models to get the best model reflect the direction of coevolution 
#-----------------------------------------------------------------------------------------------------------
### 1 Prepare the data
bird_data_ordered <- as.matrix(bird_data_ordered)
# Rename the candidate model
constrained_Model_3_list_bird_OU <- constrained_Model_3_list

### 2 Setup for parallel work
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("bird_tree_scaled", "bird_data_ordered", "constrained_Model_3_list_bird_OU"))


### 3 Run the parallel iterations of model fitting: n = 100*64
niter <- 100

result_constrained_Model_3_bird_OU_niter100 <- foreach(
  n = 1:niter, 
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(bird_tree_scaled, bird_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", # Assuming A type as DecomposablePositive
                               Syytype="Diagonal", diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- list()
  for (i in 1:length(constrained_Model_3_list_bird_OU)) {
    # Get the A matrix
    model_A <- constrained_Model_3_list_bird_OU[[i]]
    model_A[model_A == "NA"] <- NA
    
    # Fit the model
    OUOUmodel <- ouchModel(bird_tree_scaled, bird_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    iteration_result[[i]] <- OUOUmodel
  }
  return(iteration_result)
}


### 4 Save the results
save(result_constrained_Model_3_bird_OU_niter100, file = "result_constrained_Model_3_bird_OU_niter100.RData")


#----------------------------------------------------------------------------------------------
#   Step 8:  Select the best model which deciphers the directional coevolution among traits  
#----------------------------------------------------------------------------------------------
load("result_constrained_Model_3_bird_OU_niter100.RData")

### 1 Get the summary of AICc for each run
causal_constrained_bird_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:64) {
    df <- data.frame(
      AICc = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_3_bird_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Bird"
    )
    causal_constrained_bird_OU_dignosis <- rbind(causal_constrained_bird_OU_dignosis, df)
  }
}
### The best model is when Iteration = 64, Model = 58
causal_constrained_bird_OU_dignosis %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

causal_constrained_bird_OU_best <- result_constrained_Model_3_bird_OU_niter100[[64]][[58]]

### 3 Get the best A matrix
A.causal_constrained_bird_OU_best <- causal_constrained_bird_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]


### 4 Get the eigen-decomposition results of the A matrix
half.t.bird_OU <- rbind(causal_constrained_bird_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                        causal_constrained_bird_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])

#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#                                                                       #
#                                  Reptiles                             #
#                                                                       #
#-------#-------#-------#-------#-------#-------#-------#-------#-------#
#-------------------------------------
#      Step 1: Data preparation 
#-------------------------------------
### 1 Scale the maximum tree height to 1
reptile_tree <- read.tree("last_reptile_tree.tre")
tree_height <-mvSLOUCH::phyltree_paths(reptile_tree)$tree_height
reptile_tree_scaled <- reptile_tree
reptile_tree_scaled$edge.length <- reptile_tree_scaled$edge.length/tree_height

# Check the tree height = 1
mvSLOUCH::phyltree_paths(reptile_tree_scaled)$tree_height

### 2 Prepare the data to have the order of columns constrained
reptile_data <- read.csv("reptile_data_corrected.csv")
reptile_data <- reptile_data[,1:6]
reptile_data <- reptile_data %>% column_to_rownames(var = "X")
colnames(reptile_data) <- c("OS", "DT", "AL", "CS", "CF")
reptile_data <- reptile_data[, c("DT", "OS", "AL", "CS", "CF")]
reptile_data_ordered <- reptile_data[reptile_tree_scaled$tip.label, ]
reptile_data_ordered <- as.matrix(reptile_data_ordered)

#--------------------------------------------------------
#      Step 2: Covert all PPA matrix to A matrix
#--------------------------------------------------------
# Our mvOU modelling was based on the same direct/indirect structures used in PPA. 
# We thus convert the PPA matrix for each candidate model to its corresponding A matrix. See Materials and Methods for details.

### 1 Load PPA results for reptiles
load("reptile_PPA_results.RData")
length(PPA_reptiles[["model_set"]]) # = 32

### 2 Convert PPA matrix to A matrix
A.Matrix <- list()
for (i in 1:32) {
  # Get the name of the candidate model
  model_name <- names(PPA_reptiles[["model_set"]][i])
  # Get the PPA matrix
  model <- PPA_reptiles[["model_set"]][i][[1]]
  # Re-order the matrix
  model_reordered <- model[colnames(reptile_data_ordered), colnames(reptile_data_ordered)]
  # Transpose the matrix
  model_reordered_t <- t(model_reordered)
  # Conduct the transformation
  model_reordered_t[model_reordered_t == 1] <- NA
  diag(model_reordered_t) <- "+"
  mode(model_reordered_t) <- "character"
  # Get the final A matrix
  A <- model_reordered_t
  # Store them in the list
  A.Matrix[[model_name]] <- A
}

#---------------------------------------------------------------------------------------------------------
#        Step 3: Conducting mvOU modelling to fit mvOU models with varying A matrix; Syy = Diagonal
#---------------------------------------------------------------------------------------------------------

### 1 Initialize the result dataframe
Reptile_SyyDiagonal_OUOU_fit_niter500 <- data.frame()
niter <- 500  

# Set up parallel backend
detectCores()
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("reptile_tree_scaled", "reptile_data_ordered", "A.Matrix"))

Reptile_SyyDiagonal_OUOU_fit_niter500 <- foreach(
  n = 1:niter, 
  .combine = rbind,
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(reptile_tree_scaled, reptile_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", 
                               Syytype="Diagonal",  # or "UpperTri" 
                               diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- data.frame()
  for (i in 1:length(A.Matrix)) {
    # Get the A matrix
    model_A <- A.Matrix[[i]]
    # Get the model name
    nameA <- names(A.Matrix)[i]
    # Fit the model
    OUOUmodel <- ouchModel(reptile_tree_scaled, reptile_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", # or "UpperTri" 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    
    # Get the result
    result <- data.frame(
      Model = nameA,
      AIC = OUOUmodel$FinalFound$ParamSummary$aic,
      AICc = OUOUmodel$FinalFound$ParamSummary$aic.c,
      SIC = OUOUmodel$FinalFound$ParamSummary$sic,
      BIC = OUOUmodel$FinalFound$ParamSummary$bic,
      R2 = OUOUmodel$FinalFound$ParamSummary$RSS$R2,
      LogLik = OUOUmodel$FinalFound$ParamSummary$LogLik,
      DOF = OUOUmodel$FinalFound$ParamSummary$dof,
      Group = "Reptile",
      Syytype="Diagonal",
      Iteration = n,
      Seed1 = seed1,
      Seed2 = seed2
    )
    iteration_result <- rbind(iteration_result, result)
  }
  return(iteration_result)
}


stopCluster(cl)

### 2 Save the result dataframe
save(Reptile_SyyDiagonal_OUOU_fit_niter500, "Reptile_SyyDiagonal_OUOU_fit_niter500.RData")


#---------------------------------------------------------------------------------------------------------------------------------
#        Step 4: Model diagnosis: Final the best solution across various candidate models and staring points based on AICc 
#---------------------------------------------------------------------------------------------------------------------------------
### 1 Load the pre-computed mvOU outputs from the Step 3
load("Reptile_SyyDiagonal_OUOU_fit_niter500.RData")

length(unique(Reptile_SyyDiagonal_OUOU_fit_niter500$Model)) # number of models = 32

### 2 Get the best model: Model == "Model_10" & Iteration == 86; AICc = 6165.886
Reptile_SyyDiagonal_OUOU_fit_niter500 %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

#---------------------------------------------------------------------------------------------------------
#         Step 5: Generate all possible combinations of directional interactions of the best model
#---------------------------------------------------------------------------------------------------------
### 1 The best model is Model_10
A.Matrix[["Model_10"]]


### 2 get all possible combinations of directional interactions of the best model, n = 2^7 = 128
mat <- matrix(
  c("+", "NA", "0", "NA", "0",
    "0", "+", "0", "0", "0",
    "NA", "NA", "+", "NA", "0",
    "0", "0", "0", "+", "0",
    "NA", "NA", "0", "0", "+"),
  nrow = 5, byrow = TRUE,
  dimnames = list(c("DT", "OS", "AL", "CS", "CF"), 
                  c("DT", "OS", "AL", "CS", "CF"))
)

n <- nrow(mat)
non_diag <- matrix(TRUE, n, n)
diag(non_diag) <- FALSE
arrow_positions <- which(mat == "NA" & non_diag, arr.ind = TRUE)
k <- nrow(arrow_positions)
combos <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))

all_models <- vector("list", nrow(combos))

for (i in 1:nrow(combos)) {
  new_mat <- mat  
  for (j in 1:k) {
    if (combos[i, j]) {  
      row_idx <- arrow_positions[j, "row"]
      col_idx <- arrow_positions[j, "col"]
      new_mat[row_idx, col_idx] <- "0"         
      new_mat[col_idx, row_idx] <- "NA"        
    }
  }
  
  all_models[[i]] <- new_mat
}

### Name the model list
Model_10_list <- all_models

#--------------------------------------------------------------------------------------------------------------------------
#   Step 6:  Constrain the A matrix of model with different directions to reflect the sign of the PPA-OU best model
#---------------------------------------------------------------------------------------------------------------------------
# When there is a positive parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '+'
# When there is a negative parameter between 2 traits was identified in the PPA best model, the off-diagonal parameter of the corresponding A matrix will be constrained to '-'

### 1 Get the PPA model 
best_model <- best(PPA_reptiles)
PPA_OU_reptile <- best_model[["coef"]]

# Change the order of columns and rows to match with A matrix
PPA_OU_reptile  <- PPA_OU_reptile[, c("DT", "OS", "AL", "CS", "CF")]
PPA_OU_reptile  <- PPA_OU_reptile[c("DT", "OS", "AL", "CS", "CF"),]

### 2 Initialize a new list for modified matrices
constrained_Model_10_list <- vector("list", length(Model_10_list))  

# Get reference matrix dimensions
ref_rows <- rownames(PPA_OU_reptile)
ref_cols <- colnames(PPA_OU_reptile)

for (idx in seq_along(Model_10_list)) {
  mat <- Model_10_list[[idx]]
  new_mat <- mat  # Create modifiable copy
  
  # Get current matrix's row/column names
  curr_rows <- rownames(mat)
  curr_cols <- colnames(mat)
  
  # Find "NA" positions with names
  na_positions <- which(mat == "NA", arr.ind = TRUE)
  
  for (pos_idx in 1:nrow(na_positions)) {
    i <- na_positions[pos_idx, "row"]
    j <- na_positions[pos_idx, "col"]
    
    row_name <- curr_rows[i]
    col_name <- curr_cols[j]
    
    # Check same position in reference
    ref_val_same <- PPA_OU_reptile[row_name, col_name]
    
    if (ref_val_same > 0) {
      new_mat[i, j] <- "+"
    } else if (ref_val_same < 0) {
      new_mat[i, j] <- "-"
    } else {
      # Check transposed position
      ref_val_trans <- PPA_OU_reptile[col_name, row_name]
      
      if (ref_val_trans > 0) {
        new_mat[i, j] <- "+"
      } else if (ref_val_trans < 0) {
        new_mat[i, j] <- "-"
      }
      # If both zero, "NA" remains unchanged
    }
  }
  
  constrained_Model_10_list[[idx]] <- new_mat
}

#-----------------------------------------------------------------------------------------------------------
#   Step 7:  Fit all directional mvOU models to get the best model reflect the direction of coevolution 
#-----------------------------------------------------------------------------------------------------------
### 1 Prepare the data
reptile_data_ordered <- as.matrix(reptile_data_ordered)

# Rename the candidate model
constrained_Model_10_list_reptile_OU <- constrained_Model_10_list

### 2 Setup for parallel work
ncores <- 5
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, varlist = c("reptile_tree_scaled", "reptile_data_ordered", "constrained_Model_10_list_reptile_OU"))


### 3 Run the parallel iterations of model fitting: n = 100*128
niter <- 100

result_constrained_Model_10_reptile_OU_niter100 <- foreach(
  n = 1:niter, 
  .verbose = T,
  .packages = c("mvSLOUCH", "PCMBaseCpp", "dplyr")
) %dopar% {
  # library necessary pacakges
  library(dplyr)
  library(mvSLOUCH)
  library(PCMBaseCpp)
  
  # Set seed for the starting points model fitting: seed1
  seed1 <- 2025+n
  set.seed(seed1)
  
  # This starting model might take a long time, because we didn't specify any model structure
  OUOUmodel_start <- ouchModel(reptile_tree_scaled, reptile_data_ordered, regimes=NULL, 
                               Atype="DecomposablePositive", # Assuming A type as DecomposablePositive
                               Syytype="Diagonal", diagA=NULL, 
                               maxiter=c(25,50))
  # Set seed for the final model: seed2
  seed2 <- 20250+n
  set.seed(seed2)
  
  # Set data frame for inner layer of iteration
  iteration_result <- list()
  for (i in 1:length(constrained_Model_10_list_reptile_OU)) {
    # Get the A matrix
    model_A <- constrained_Model_10_list_reptile_OU[[i]]
    model_A[model_A == "NA"] <- NA
    
    # Fit the model
    OUOUmodel <- ouchModel(reptile_tree_scaled, reptile_data_ordered, 
                           regimes=NULL, 
                           Atype="Any",  
                           Syytype="Diagonal", 
                           diagA=NULL, 
                           maxiter=c(100,1000),
                           start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), 
                           parameter_signs=list(signsA=model_A)
    )
    iteration_result[[i]] <- OUOUmodel
  }
  return(iteration_result)
}


stopCluster(cl)

### 4 Save the results
save(result_constrained_Model_10_reptile_OU_niter100, file = "result_constrained_Model_10_reptile_OU_niter100.RData")

#----------------------------------------------------------------------------------------------
#   Step 8:  Select the best model which deciphers the directional coevolution among traits  
#----------------------------------------------------------------------------------------------
load("result_constrained_Model_10_reptile_OU_niter100.RData")

### 1 Get the summary of AICc for each run
causal_constrained_reptile_OU_dignosis <- data.frame()
for (n in 1:100) {
  for (i in 1:128) {
    df <- data.frame(
      AICc = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic.c"]],
      AIC = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["aic"]],
      BIC = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["bic"]],
      LogLik = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["LogLik"]],
      R2 = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["RSS"]][["R2"]],
      DOF = result_constrained_Model_10_reptile_OU_niter100[[n]][[i]][["FinalFound"]][["ParamSummary"]][["dof"]],
      Iteration = n,
      Model = i,
      Group = "Reptile"
    )
    causal_constrained_reptile_OU_dignosis <- rbind(causal_constrained_reptile_OU_dignosis, df)
  }
}

### 2 The best model is when Iteration = 52, Model = 83
causal_constrained_reptile_OU_dignosis %>%
  slice_min(order_by = AICc, n = 1, with_ties = TRUE)

causal_constrained_reptile_OU_best <- result_constrained_Model_10_reptile_OU_niter100[[52]][[83]]

### 3 Get the best A matrix
A.causal_constrained_reptile_OU_best <- causal_constrained_reptile_OU_best[["FinalFound"]][["ParamsInModel"]][["A"]]

### 4 Get the eigen-decomposition results of the A matrix
half.t.reptile_OU <- rbind(causal_constrained_reptile_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["directions"]],
                           causal_constrained_reptile_OU_best[["FinalFound"]][["ParamSummary"]][["phyl.halflife"]][["halflives"]])


