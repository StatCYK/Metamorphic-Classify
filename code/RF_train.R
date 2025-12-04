# Load required libraries
library(tidyverse)
library(ranger)
library(caret)
library(pROC)
library(doParallel)
library(foreach)
# Set seed for reproducibility
set.seed(123)

# Read the data
df_all <- read.csv("./data/Metamorphic_classify.csv")

# Prepare training data
df_train <- df_all %>% 
  filter(dataset %in% c("FS", "atlas", "codnas")) %>% 
  select(-1)  # Remove index if first column is index

# Prepare test data (pdb dataset)
df_test <- df_all %>% 
  filter(dataset == "pdb") %>% 
  select(-1)

# Get feature columns for training
X_df <- df_train %>% 
  select(-dataset) %>% 
  select(-1)  # Remove first column of jobname

# Get feature columns for testing
X_test_df <- df_test %>% 
  select(-dataset) %>% 
  select(-1)  # Remove first column of jobname

# Create binary target 
y <- as.integer(df_train$dataset == "FS")
y_test <- as.integer(df_test$dataset == "FS")  # This will be all 0 since test is pdb

# Calculate class weights
calculate_class_weights <- function(y) {
  classes <- unique(y)
  n_samples <- length(y)
  n_classes <- length(classes)
  
  class_counts <- table(y)
  weights <- n_samples / (n_classes * class_counts)
  
  return(list(weights = weights, classes = classes))
}

# Get class weights
class_weight_result <- calculate_class_weights(y)
class_weights <- as.numeric(class_weight_result$weights)

# Create weight_list
alphas <- c(1,2,4,6,8,10)
weight_list <- list()

for (alpha in alphas) {
  weight_list[[as.character(alpha)]] <- c(class_weights[1]/(1 + alpha), alpha * class_weights[2]/(1 + alpha))
}

# Define parameter grid
param_grid_RF <- expand.grid(
  num.trees = 500,
  max.depth = c(5,10,15),
  min.node.size = c(6,8,10,12),  
  class_weight = names(weight_list),
  replace = TRUE,  
  impurity_reduction = c(0, 0.01), 
  stringsAsFactors = TRUE
)

# Function to find best hyperparameters from CV results
find_best_hyperparams_from_cv <- function(cv_results) {
  best_params_list <- list()
  
  for (fold_res in cv_results) {
    if (!inherits(fold_res, "error")) {
      best_params_list[[length(best_params_list) + 1]] <- fold_res$best_params
    }
  }
  
  # Average the best parameters across folds
  if (length(best_params_list) > 0) {
    avg_num_trees <- round(mean(sapply(best_params_list, function(x) x$num.trees)))
    avg_max_depth <- round(mean(sapply(best_params_list, function(x) x$max.depth)))
    avg_min_node_size <- round(mean(sapply(best_params_list, function(x) x$min.node.size)))
    
    # Find most frequent class_weight
    class_weights <- sapply(best_params_list, function(x) as.character(x$class_weight))
    most_freq_weight <- names(sort(table(class_weights), decreasing = TRUE))[1]
    
    avg_impurity_reduction <- mean(sapply(best_params_list, function(x) x$impurity_reduction))
    
    return(list(
      num.trees = avg_num_trees,
      max.depth = avg_max_depth,
      min.node.size = avg_min_node_size,
      class_weight = most_freq_weight,
      replace = TRUE,
      impurity_reduction = avg_impurity_reduction
    ))
  }
  
  return(NULL)
}

# Revised function with parallel computing and threshold balancing
perform_nested_cv <- function(X, y, X_test, y_test, param_grid, weight_list, outer_folds = 5, inner_folds = 5) {
  library(doParallel)
  library(foreach)
  
  # Create outer stratified folds
  outer_skf <- createFolds(factor(y), k = outer_folds, list = TRUE, returnTrain = FALSE)
  
  # Initialize result lists
  train_class0_errors <- numeric(outer_folds)
  train_class1_errors <- numeric(outer_folds)
  test_class0_errors <- numeric(outer_folds)
  test_class1_errors <- numeric(outer_folds)
  test_aucs <- numeric(outer_folds)
  
  # Track misclassification
  misclass_scores_class0 <- numeric(outer_folds)
  misclass_scores_class1 <- numeric(outer_folds)
  misclass_indices_class0 <- numeric(outer_folds)
  misclass_indices_class1 <- numeric(outer_folds)
  
  # Lists to store results for plotting
  y_test_list <- list()
  y_test_prob_list <- list()
  model_list <- list()
  train_confusion_matrices <- list()
  test_confusion_matrices <- list()
  
  # Store all fold results for later analysis
  cv_results <- list()
  
  # Custom function to calculate balanced accuracy
  balanced_accuracy <- function(y_true, y_pred) {
    cm <- table(y_true, y_pred)
    if (nrow(cm) == 2 && ncol(cm) == 2) {
      sensitivity <- cm[2, 2] / sum(cm[2, ])  # TP / (TP + FN)
      specificity <- cm[1, 1] / sum(cm[1, ])  # TN / (TN + FP)
      return((sensitivity + specificity) / 2)
    }
    return(0)
  }
  
  # Outer CV loop - parallelize outer folds
  outer_results <- foreach(fold = 1:outer_folds, 
                           .combine = list, 
                           .multicombine = TRUE,
                           .packages = c("ranger", "pROC", "caret", "foreach", "doParallel"),
                           .errorhandling = "pass") %dopar% {

                             cat(sprintf("\n=== Outer Fold %d/%d ===\n", fold, outer_folds))
                             
                             # Split data for outer fold
                             test_idx <- outer_skf[[fold]]
                             train_idx <- setdiff(1:length(y), test_idx)
                             
                             X_train_outer <- X[train_idx, ]
                             y_train_outer <- y[train_idx]
                             X_test_outer <- X[test_idx, ]
                             y_test_outer <- y[test_idx]
                             
                             # Inner CV for hyperparameter tuning
                             cat("  Performing inner CV for hyperparameter tuning...\n")
                             
                             # Create inner folds
                             inner_skf <- createFolds(factor(y_train_outer), k = inner_folds, list = TRUE, returnTrain = FALSE)
                             
                             # Parallel inner CV grid search
                             inner_grid_results <- foreach(i = 1:nrow(param_grid), 
                                                           .combine = rbind,
                                                           .packages = c("ranger", "pROC", "caret", "foreach", "doParallel"),
                                                           .errorhandling = "remove") %dopar% {
                                                             
                                                             params <- param_grid[i, ]
                                                             weights <- weight_list[[params$class_weight]]
                                                             
                                                             # Evaluate this parameter set across all inner folds
                                                             fold_scores <- numeric(inner_folds)
                                                             
                                                             for (inner_fold in 1:inner_folds) {
                                                               val_idx <- inner_skf[[inner_fold]]
                                                               inner_train_idx <- setdiff(1:length(y_train_outer), val_idx)
                                                               
                                                               X_train_inner <- X_train_outer[inner_train_idx, ]
                                                               y_train_inner <- y_train_outer[inner_train_idx]
                                                               X_val_inner <- X_train_outer[val_idx, ]
                                                               y_val_inner <- y_train_outer[val_idx]
                                                               
                                                               # Train data for inner fold
                                                               inner_data <- cbind(X_train_inner, target = as.factor(y_train_inner))
                                                               
                                                               # Train model with current parameters
                                                               model_inner <- ranger(
                                                                 formula = target ~ .,
                                                                 data = inner_data,
                                                                 num.trees = params$num.trees,
                                                                 max.depth = params$max.depth,
                                                                 min.node.size = params$min.node.size,
                                                                 probability = TRUE,
                                                                 class.weights = weights,
                                                                 replace = params$replace,
                                                                 importance = "none",
                                                                 respect.unordered.factors = "order",
                                                                 seed = 123 + inner_fold,
                                                                 verbose = FALSE
                                                               )
                                                               
                                                               if (!is.null(model_inner)) {
                                                                 # Predict on validation set
                                                                 val_pred <- predict(model_inner, X_val_inner)$predictions[, 2]
                                                                 val_pred_class <- ifelse(val_pred >= 0.5, 1, 0)
                                                                 # Calculate balanced accuracy
                                                                 fold_scores[inner_fold] <- balanced_accuracy(y_val_inner, val_pred_class)
                                                               } else {
                                                                 fold_scores[inner_fold] <- NA
                                                               }
                                                             }
                                                             
                                                             # Return average score for this parameter set
                                                             data.frame(
                                                               num.trees = params$num.trees,
                                                               max.depth = params$max.depth,
                                                               min.node.size = params$min.node.size,
                                                               class_weight = params$class_weight,
                                                               replace = params$replace,
                                                               impurity_reduction = params$impurity_reduction,
                                                               mean_score = mean(fold_scores, na.rm = TRUE),
                                                               stringsAsFactors = FALSE
                                                             )
                                                           }
                             
                             # Find best parameters (highest mean balanced accuracy)
                             best_params_idx <- which.max(inner_grid_results$mean_score)
                             best_params <- inner_grid_results[best_params_idx, ]
                             
                             cat(sprintf("  Best params: trees=%d, depth=%d, min.node=%d, weight=%s, imp.red=%.2f, score=%.3f\n",
                                         best_params$num.trees, best_params$max.depth, best_params$min.node.size,
                                         best_params$class_weight, best_params$impurity_reduction, best_params$mean_score))
                             
                             # Train final model on outer training set with best parameters
                             weights <- weight_list[[best_params$class_weight]]
                             final_data <- cbind(X_train_outer, target = as.factor(y_train_outer))
                             
                             model_final <- ranger(
                               formula = target ~ .,
                               data = final_data,
                               num.trees = best_params$num.trees,
                               max.depth = best_params$max.depth,
                               min.node.size = best_params$min.node.size,
                               probability = TRUE,
                               class.weights = weights,
                               replace = best_params$replace,
                               importance = "impurity",
                               respect.unordered.factors = "order",
                               seed = 123,
                               verbose = FALSE
                             )
                             
                             # Get probabilities
                             y_train_prob <- predict(model_final, X_train_outer)$predictions[, 2]
                             y_test_prob <- predict(model_final, X_test_outer)$predictions[, 2]
                             
                             # Standard predictions (0.5 threshold)
                             y_train_pred <- ifelse(y_train_prob >= 0.5, 1, 0)
                             y_test_pred <- ifelse(y_test_prob >= 0.5, 1, 0)
                             
                             # Calculate metrics
                             test_auc <- auc(roc(response = y_test_outer, predictor = y_test_prob))
                             
                             # Confusion Matrices
                             train_cm <- table(y_train_outer, y_train_pred)
                             test_cm <- table(y_test_outer, y_test_pred)
                             
                             # Class-wise errors
                             train_class0_error <- ifelse(sum(y_train_outer == 0) > 0,
                                                          sum(y_train_pred[y_train_outer == 0] != 0) / sum(y_train_outer == 0),
                                                          NA)
                             train_class1_error <- ifelse(sum(y_train_outer == 1) > 0,
                                                          sum(y_train_pred[y_train_outer == 1] != 1) / sum(y_train_outer == 1),
                                                          NA)
                             
                             test_class0_error <- ifelse(sum(y_test_outer == 0) > 0,
                                                         sum(y_test_pred[y_test_outer == 0] != 0) / sum(y_test_outer == 0),
                                                         NA)
                             test_class1_error <- ifelse(sum(y_test_outer == 1) > 0,
                                                         sum(y_test_pred[y_test_outer == 1] != 1) / sum(y_test_outer == 1),
                                                         NA)
                             
                             # Identify most misclassified points
                             misclass_error <- abs(y_test_prob - y_test_outer)
                             
                             misclass_idx_0 <- NA
                             misclass_score_0 <- NA
                             misclass_idx_1 <- NA
                             misclass_score_1 <- NA
                             
                             # Class 0
                             if (any(y_test_outer == 0)) {
                               idx0 <- which.max(misclass_error[y_test_outer == 0])
                               true_idx0 <- test_idx[which(y_test_outer == 0)[idx0]]
                               misclass_idx_0 <- true_idx0
                               misclass_score_0 <- misclass_error[y_test_outer == 0][idx0]
                             }
                             
                             # Class 1
                             if (any(y_test_outer == 1)) {
                               idx1 <- which.max(misclass_error[y_test_outer == 1])
                               true_idx1 <- test_idx[which(y_test_outer == 1)[idx1]]
                               misclass_idx_1 <- true_idx1
                               misclass_score_1 <- misclass_error[y_test_outer == 1][idx1]
                             }
                             
                             # Return fold results
                             list(
                               fold = fold,
                               model = model_final,
                               best_params = best_params,
                               y_test = y_test_outer,
                               y_test_prob = y_test_prob,
                               test_auc = test_auc,
                               train_cm = train_cm,
                               test_cm = test_cm,
                               #all_val_labels = all_val_labels,
                               train_class0_error = train_class0_error,
                               train_class1_error = train_class1_error,
                               test_class0_error = test_class0_error,
                               test_class1_error = test_class1_error,
                               misclass_idx_0 = misclass_idx_0,
                               misclass_score_0 = misclass_score_0,
                               misclass_idx_1 = misclass_idx_1,
                               misclass_score_1 = misclass_score_1
                             )
                           }
  

  # Process parallel results
  for (fold_res in outer_results) {
    if (inherits(fold_res, "error")) {
      cat("Error in fold:", fold_res$message, "\n")
      next
    }
    
    fold <- fold_res$fold
    
    # Store results
    cv_results[[fold]] <- fold_res
    model_list[[fold]] <- fold_res$model
    y_test_list[[fold]] <- fold_res$y_test
    y_test_prob_list[[fold]] <- fold_res$y_test_prob
    
    test_aucs[fold] <- fold_res$test_auc
    train_confusion_matrices[[fold]] <- fold_res$train_cm
    test_confusion_matrices[[fold]] <- fold_res$test_cm
    
    train_class0_errors[fold] <- fold_res$train_class0_error
    train_class1_errors[fold] <- fold_res$train_class1_error
    test_class0_errors[fold] <- fold_res$test_class0_error
    test_class1_errors[fold] <- fold_res$test_class1_error
    
    misclass_indices_class0[fold] <- fold_res$misclass_idx_0
    misclass_scores_class0[fold] <- fold_res$misclass_score_0
    misclass_indices_class1[fold] <- fold_res$misclass_idx_1
    misclass_scores_class1[fold] <- fold_res$misclass_score_1
  }
  all_val_probs = unlist(y_test_prob_list)
  all_val_labels <- unlist(y_test_list)
  
  # ============ CRITICAL: Find shared optimal threshold ============
  cat("\n=== Finding threshold to strike a balanced type I&II error ===\n")
  best_threshold <- 0.5
  best_balance <- Inf
  
  for (threshold in seq(0.1, 0.9, 0.02)) {
    y_val_pred_thresh <- as.integer(all_val_probs >= threshold)
    
    if (length(unique(y_val_pred_thresh)) > 1 && length(unique(all_val_labels)) > 1) {
      cm <- table(all_val_labels, y_val_pred_thresh)
      if (nrow(cm) == 2 && ncol(cm) == 2) {
        tn <- cm[1, 1]
        fp <- cm[1, 2]
        fn <- cm[2, 1]
        tp <- cm[2, 2]
        
        type1_error <- ifelse((fp + tn) > 0, fp / (fp + tn), 0)
        type2_error <- ifelse((fn + tp) > 0, fn / (fn + tp), 0)
        
        balance_metric <- abs(type1_error - type2_error)
        #print(balance_metric)
        if (balance_metric < best_balance) {
          
          best_balance <- balance_metric
          best_threshold <- threshold
        }
      }
    }
  }
  
  optimal_shared_threshold <- best_threshold
  cat(sprintf("Shared optimal threshold found: %.3f\n", optimal_shared_threshold))
  
  # Apply the shared threshold to get balanced predictions for each fold
  cat("\n=== Applying shared threshold to all folds ===\n")
  
  balanced_test_confusion_matrices <- list()
  balanced_test_class0_errors <- numeric(outer_folds)
  balanced_test_class1_errors <- numeric(outer_folds)
  
  for (fold in 1:outer_folds) {
    y_test <- y_test_list[[fold]]
    y_test_prob <- y_test_prob_list[[fold]]
    
    # Apply shared optimal threshold
    y_test_pred_balanced <- as.integer(y_test_prob >= optimal_shared_threshold)
    
    # Calculate balanced metrics
    balanced_test_cm <- table(y_test, y_test_pred_balanced)
    balanced_test_confusion_matrices[[fold]] <- balanced_test_cm
    
    balanced_test_class0_errors[fold] <- ifelse(sum(y_test == 0) > 0,
                                                sum(y_test_pred_balanced[y_test == 0] != 0) / sum(y_test == 0),
                                                NA)
    balanced_test_class1_errors[fold] <- ifelse(sum(y_test == 1) > 0,
                                                sum(y_test_pred_balanced[y_test == 1] != 1) / sum(y_test == 1),
                                                NA)
    
  }
  
  # Calculate average confusion matrices
  avg_confusion_matrix <- function(matrices) {
    if (length(matrices) == 0) return(NULL)
    
    # Initialize with first matrix
    avg_cm <- matrices[[1]] * 0
    
    # Sum all matrices
    for (cm in matrices) {
      # Ensure same dimensions
      if (all(dim(cm) == dim(avg_cm))) {
        avg_cm <- avg_cm + cm
      }
    }
    
    # Divide by number of folds
    return(avg_cm / length(matrices))
  }
  
  avg_train_cm <- avg_confusion_matrix(train_confusion_matrices)
  avg_test_cm <- avg_confusion_matrix(test_confusion_matrices)
  avg_balanced_test_cm <- avg_confusion_matrix(balanced_test_confusion_matrices)
  
  # ============ Train final model on ALL training data ============
  cat("\n=== Training final model on ALL training data ===\n")
  # Create  folds on full data
  full_skf <-  createFolds(factor(y), k = inner_folds, list = TRUE, returnTrain = FALSE)
  
  
  # Parallel CV grid search on full data
  full_data_results <- foreach(i = 1:nrow(param_grid), 
                               .combine = rbind,
                               .packages = c("ranger", "pROC", "caret", "foreach", "doParallel"),
                               .errorhandling = "remove") %dopar% {
                                 # Custom function to calculate balanced accuracy
                                 balanced_accuracy <- function(y_true, y_pred) {
                                   cm <- table(y_true, y_pred)
                                   if (nrow(cm) == 2 && ncol(cm) == 2) {
                                     sensitivity <- cm[2, 2] / sum(cm[2, ])  # TP / (TP + FN)
                                     specificity <- cm[1, 1] / sum(cm[1, ])  # TN / (TN + FP)
                                     return((sensitivity + specificity) / 2)
                                   }
                                   return(0)
                                 }
                                 params <- param_grid[i, ]
                                 weights <- weight_list[[params$class_weight]]
                                 
                                 # Evaluate this parameter set across all inner folds
                                 fold_scores <- numeric(inner_folds)
                                 
                                 for (inner_fold in 1:inner_folds) {
                                   val_idx <- full_skf[[inner_fold]]
                                   inner_train_idx <- setdiff(1:length(y), val_idx)
                                   
                                   X_train_inner <- X[inner_train_idx, ]
                                   y_train_inner <- y[inner_train_idx]
                                   X_val_inner <- X[val_idx, ]
                                   y_val_inner <- y[val_idx]
                                   
                                   # Train data for inner fold
                                   inner_data <- cbind(X_train_inner, target = as.factor(y_train_inner))
                                   
                                   # Train model with current parameters
                                   model_inner <- ranger(
                                     formula = target ~ .,
                                     data = inner_data,
                                     num.trees = params$num.trees,
                                     max.depth = params$max.depth,
                                     min.node.size = params$min.node.size,
                                     probability = TRUE,
                                     class.weights = weights,
                                     replace = params$replace,
                                     importance = "none",
                                     respect.unordered.factors = "order",
                                     seed = 123 + inner_fold,
                                     verbose = FALSE
                                   )
                                   
                                   if (!is.null(model_inner)) {
                                     # Predict on validation set
                                     val_pred <- predict(model_inner, X_val_inner)$predictions[, 2]
                                     val_pred_class <- ifelse(val_pred >= 0.5, 1, 0)
                                     # Calculate balanced accuracy
                                     fold_scores[inner_fold] <- balanced_accuracy(y_val_inner, val_pred_class)
                                   } else {
                                     fold_scores[inner_fold] <- NA
                                   }
                                 }
                                 
                                 # Return average score for this parameter set
                                 data.frame(
                                   num.trees = params$num.trees,
                                   max.depth = params$max.depth,
                                   min.node.size = params$min.node.size,
                                   class_weight = params$class_weight,
                                   replace = params$replace,
                                   impurity_reduction = params$impurity_reduction,
                                   mean_score = mean(fold_scores, na.rm = TRUE),
                                   stringsAsFactors = FALSE
                                 )
                               }
  
  
  # Find best parameters (highest mean balanced accuracy)
  best_params_idx <- which.max(full_data_results$mean_score)
  best_params <- full_data_results[best_params_idx, ]
  
  cat(sprintf("  Best params: trees=%d, depth=%d, min.node=%d, weight=%s, imp.red=%.2f, score=%.3f\n",
              best_params$num.trees, best_params$max.depth, best_params$min.node.size,
              best_params$class_weight, best_params$impurity_reduction, best_params$mean_score))
  
  # Train final model on outer training set with best parameters
  weights <- weight_list[[best_params$class_weight]]
  final_data <- cbind(X, target = as.factor(y))
  
  final_model <- ranger(
    formula = target ~ .,
    data = final_data,
    num.trees = best_params$num.trees,
    max.depth = best_params$max.depth,
    min.node.size = best_params$min.node.size,
    probability = TRUE,
    class.weights = weights,
    replace = best_params$replace,
    importance = "impurity",
    respect.unordered.factors = "order",
    seed = 123,
    verbose = FALSE
  )
  # Test on pdb dataset with both thresholds
  cat("\n=== Testing on pdb dataset ===\n")
  
  # Get feature importances
  feature_importances <- final_model$variable.importance
  cat("\nTop 10 most important features:\n")
  sorted_importances <- sort(feature_importances, decreasing = TRUE)
  print(head(sorted_importances, 10))
  
  # Create results list
  results <- list(
    # Standard threshold results
    train_class0_mean = mean(train_class0_errors, na.rm = TRUE),
    train_class0_sd = sd(train_class0_errors, na.rm = TRUE),
    train_class1_mean = mean(train_class1_errors, na.rm = TRUE),
    train_class1_sd = sd(train_class1_errors, na.rm = TRUE),
    test_class0_mean = mean(test_class0_errors, na.rm = TRUE),
    test_class0_sd = sd(test_class0_errors, na.rm = TRUE),
    test_class1_mean = mean(test_class1_errors, na.rm = TRUE),
    test_class1_sd = sd(test_class1_errors, na.rm = TRUE),
    test_auc_mean = mean(test_aucs, na.rm = TRUE),
    test_auc_sd = sd(test_aucs, na.rm = TRUE),
    train_confusion_matrices = train_confusion_matrices,
    test_confusion_matrices = test_confusion_matrices,
    avg_train_confusion_matrix = avg_train_cm,
    avg_test_confusion_matrix = avg_test_cm,
    
    # Shared threshold results
    optimal_shared_threshold = optimal_shared_threshold,
    balanced_test_class0_mean = mean(balanced_test_class0_errors, na.rm = TRUE),
    balanced_test_class0_sd = sd(balanced_test_class0_errors, na.rm = TRUE),
    balanced_test_class1_mean = mean(balanced_test_class1_errors, na.rm = TRUE),
    balanced_test_class1_sd = sd(balanced_test_class1_errors, na.rm = TRUE),
    balanced_test_confusion_matrices = balanced_test_confusion_matrices,
    avg_balanced_test_confusion_matrix = avg_balanced_test_cm,
    
    # For plotting
    y_test_list = y_test_list,
    y_test_prob_list = y_test_prob_list,
    model_list = model_list,
    
    # Misclassification info
    misclassified_info = list(
      most_misclassified_class0_idx = if (any(!is.na(misclass_indices_class0))) 
        misclass_indices_class0[which.max(misclass_scores_class0)] else NULL,
      most_misclassified_class0_score = if (any(!is.na(misclass_scores_class0))) 
        max(misclass_scores_class0, na.rm = TRUE) else NULL,
      most_misclassified_class1_idx = if (any(!is.na(misclass_indices_class1))) 
        misclass_indices_class1[which.max(misclass_scores_class1)] else NULL,
      most_misclassified_class1_score = if (any(!is.na(misclass_scores_class1))) 
        max(misclass_scores_class1, na.rm = TRUE) else NULL
    ),
    
    # Final model results
    final_model = if (exists("final_model")) final_model else NULL,
    feature_importances = if (exists("feature_importances")) feature_importances else NULL,
    test_probabilities = if (exists("test_pred_probs")) test_pred_probs else NULL,
    best_hyperparams = best_params,
    # CV results for analysis
    cv_results = cv_results
  )
  
  return(results)
}

# Convert X_df to matrix for consistency
X_matrix <- X_df


# Setup parallel backend
cl <- makePSOCKcluster(80)
registerDoParallel(cl)

# Run the nested CV with test set evaluation
cat("Starting nested CV procedure (5 outer folds, 5 inner folds)...\n")
results <- perform_nested_cv(
  X = X_matrix,
  y = y,
  X_test = X_test_df,
  y_test = y_test,
  param_grid = param_grid_RF,
  weight_list = weight_list,
  outer_folds = 5,
  inner_folds = 5
)


# Stop parallel cluster
stopCluster(cl)
# Print comprehensive performance comparison
cat("\n", rep("=", 70), "\n", sep = "")
cat("PERFORMANCE COMPARISON\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Cross-Validation Results:\n")
cat("Standard threshold (0.5):\n")
cat(sprintf("  Class 0 error: %.3f ± %.3f\n", 
            results$test_class0_mean, results$test_class0_sd))
cat(sprintf("  Class 1 error: %.3f ± %.3f\n", 
            results$test_class1_mean, results$test_class1_sd))
cat(sprintf("  Error difference: %.3f\n", 
            abs(results$test_class0_mean - results$test_class1_mean)))
cat(sprintf("  AUC: %.3f ± %.3f\n", 
            results$test_auc_mean, results$test_auc_sd))

cat(sprintf("\nBalanced threshold (%.3f):\n", results$optimal_shared_threshold))
cat(sprintf("  Class 0 error: %.3f ± %.3f\n", 
            results$balanced_test_class0_mean, results$balanced_test_class0_sd))
cat(sprintf("  Class 1 error: %.3f ± %.3f\n", 
            results$balanced_test_class1_mean, results$balanced_test_class1_sd))
cat(sprintf("  Error difference: %.3f\n", 
            abs(results$balanced_test_class0_mean - results$balanced_test_class1_mean)))

if (!is.null(results$test_accuracy_standard)) {
  cat("\nTest Results on pdb dataset:\n")
  cat(sprintf("  Test accuracy (0.5 threshold): %.3f\n", results$test_accuracy_standard))
  cat(sprintf("  Test accuracy (balanced threshold): %.3f\n", results$test_accuracy_balanced))
  cat(sprintf("  Mean prediction (0.5 threshold): %.3f\n", mean(results$test_predictions_standard)))
  cat(sprintf("  Mean prediction (balanced threshold): %.3f\n", mean(results$test_predictions_balanced)))
}

# Show confusion matrices
cat("\n", rep("=", 70), "\n", sep = "")
cat("CONFUSION MATRICES (Average across folds)\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Standard threshold (0.5) - Average test CM:\n")
print(results$avg_test_confusion_matrix)

cat(sprintf("\nBalanced threshold (%.3f) - Average test CM:\n", results$optimal_shared_threshold))
print(results$avg_balanced_test_confusion_matrix)

# Save results
saveRDS(results, "nested_cv_mia_results.rds")
cat("\nResults saved to 'nested_cv_mia_results.rds'\n")

create_plots <- function(results) {
  library(ggplot2)
  library(gridExtra)
  library(pROC)
  library(RColorBrewer)
  
  # Set up colors for 5 folds
  fold_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")  # ColorBrewer Set1
  fold_names <- c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5")
  
  # Create layout
  par(mfrow = c(1, 2))
  
  # Plot 1: ROC Curves with individual labels
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = "Receiver Operating Characteristic (ROC) Curves")
  
  auc_values <- numeric()
  auc_labels <- character()
  
  # Plot each fold's ROC curve with different colors and store AUC
  for (i in 1:length(results$y_test_list)) {
    if (!is.null(results$y_test_list[[i]]) && !is.null(results$y_test_prob_list[[i]])) {
      roc_obj <- roc(results$y_test_list[[i]], results$y_test_prob_list[[i]])
      auc_val <- auc(roc_obj)
      auc_values <- c(auc_values, auc_val)
      
      # Create label with fold number and AUC
      auc_label <- sprintf("Fold %d (AUC = %.3f)", i, auc_val)
      auc_labels <- c(auc_labels, auc_label)
      
      # Plot ROC curve with fold-specific color
      lines(1 - roc_obj$specificities, roc_obj$sensitivities, 
            col = fold_colors[i], lwd = 2, type = "l")
    }
  }
  
  # Add diagonal line
  abline(0, 1, col = "gray40", lty = 2, lwd = 1.5)
  
  # Add legend with all folds
  if (length(auc_labels) > 0) {
    legend("bottomright", 
           legend = auc_labels,
           col = fold_colors[1:length(auc_labels)],
           lwd = 2,
           cex = 0.8,
           bg = "white")
  }
  
  # Add panel label
  mtext("a", side = 3, line = -2, at = -0.2, cex = 1.5, font = 2)
  
  # Plot 2: Feature Importances
  if (!is.null(results$feature_importances)) {
    importances <- results$feature_importances
    sorted_idx <- order(importances, decreasing = TRUE)
    
    # Get feature names
    feature_names <- names(importances)
    
    # Create bar plot
    bar_colors <- rep("#186A3B", length(importances))
    
    # Adjust margins for better feature name display
    par(mar = c(7, 4, 4, 2) + 0.1)  # Increase bottom margin
    
    barplot(importances[sorted_idx], 
            names.arg = feature_names[sorted_idx],
            las = 2, cex.names = 0.7,
            ylab = "Mean Decrease in Impurity",
            main = "Feature Importances",
            col = bar_colors,
            ylim = c(0, max(importances) * 1.1))
    
    # Add panel label
    mtext("b", side = 3, line = -2, at = -1.5, cex = 1.5, font = 2)
    
    # Reset margins
    par(mar = c(5, 4, 4, 2) + 0.1)
  } else {
    plot(0, 0, type = "n", xlab = "", ylab = "", 
         main = "Feature Importances\n(No final model trained)",
         xlim = c(-1, 1), ylim = c(-1, 1))
    text(0, 0, "Run nested CV first to train final model", cex = 0.8)
  }
  
  par(mfrow = c(1, 1))
}
cat("\n=== Creating plots ===\n")
create_plots(results)


# Print feature importances summary
if (!is.null(results$feature_importances)) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("FEATURE IMPORTANCES (from final model on all training data)\n")
  cat(rep("=", 60), "\n\n", sep = "")
  
  importances <- results$feature_importances
  sorted_importances <- sort(importances, decreasing = TRUE)
  
  for (i in 1:min(10, length(sorted_importances))) {
    cat(sprintf("%2d. %-30s: %.4f\n", i, names(sorted_importances)[i], sorted_importances[i]))
  }
}


# Test on pdb dataset
cat("\n=== Testing on pdb dataset ===\n")
test_pred_probs <- predict(results$final_model, X_test_df)$predictions[, 2]
df_test$pred_prob = test_pred_probs 

# Sort by score (descending) with dplyr
df_test <- df_test %>% arrange(desc(pred_prob))
head(df_test)
