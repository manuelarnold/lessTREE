#' Regularize SEM Trees
#'
#' Regularize SEM trees fitted with semtree. Requires tidyverse and lessSEM.
#'
#' @param tree A semtree object. Only SEM trees fitted with lavaan work.
#' @param regularized optional: Select the parameters that should be regularized.
#' If set to NULL, all parameters with user-defined names will be regularized.
#' @return A list.
#' @export
regularize_semtree <- function(tree, regularized = NULL) {

  # Extracting the models and their parameters ----

  leafs <- getLeafs(tree)
  n_leafs <- length(leafs)

  base_parameters <- c()
  is_valid_parameter <- c()
  parameters <- c()
  models <- c()

  for(leaf in 1:n_leafs){

    model <- leafs[[leaf]]$model

    model_parameters <- get_parameter_labels(model = model,
                                             leaf = leaf)

    base_parameters <- c(base_parameters,
                         unique(model_parameters$base_labels))

    is_valid_parameter <- c(is_valid_parameter,
                            model_parameters$valid_labels[unique(model_parameters$base_labels)]
    )

    models <- c(models,
                lavaan(model = model_parameters$parameter_table,
                       data = lavInspect(object = model, what = "data"))
    )

    # get final parameter names for this model:
    parameters <- c(parameters,
                    names(getLavaanParameters(models[[leaf]])))

  }

  ## Provide some feedback to the user:
  cat("\n")
  cat("Setting up lessTREE. The following parameters will be regularized: ")
  cat(paste0(unique(base_parameters[is_valid_parameter]), collapse = ", "))
  cat("\n")

  if(any(!is_valid_parameter)){
    cat(crayon::red("✖ "),
        "The following parameters will be estimated group-specific,",
                     "but cannot be regularized because their",
                     "names are inconsistent with the naming convention used in lessSEM: ")
    cat(crayon::red(paste0(unique(base_parameters[!is_valid_parameter]), collapse = ", ")))
    cat("\nAll names must start with letters and special characters (e.g., ",
                     "~, -, or >) are not allowed. Consider changing the names of these parameters.\n")
  }

  # Defining the transformations ----

  delta_parameters <- paste0("delta_", parameters)

  # check if only a subset of the parameters should be regularized
  if(!is.null(regularized)){
    regularize <- base_parameters %in% regularized
  }else{
    regularize <- rep(TRUE, length(base_parameters))
  }

  transformation <- paste0(
    "parameters: ",
    paste0(unique(c(base_parameters[is_valid_parameter & regularize],
                    delta_parameters[is_valid_parameter & regularize],
                    parameters[is_valid_parameter & regularize])),
           collapse = ", "),
    "\n\n")

  for(i in 1:length(delta_parameters)){
    # skip parameters with invalid names:
    if(!is_valid_parameter[i] | !regularize[i])
      next
    transformation <- paste0(
      transformation,
      paste0(parameters[i], " = ", base_parameters[i], " + ", delta_parameters[i],
             ";"),
      "\n")
  }

  # return new object which can then be regularized:
  base_tree <- list(
    lavaanModel = models,
    modifyModel = modifyModel(transformations = transformation),
    regularized = delta_parameters[is_valid_parameter & regularize],
    # only those base parameters for regularized parameters; we don't really care about the rest
    base_parameters = unique(base_parameters[is_valid_parameter & regularize]),
    n_leafs = n_leafs,
    lessSEM_object = NULL, # we will save the fitted lessSEM here
    unregularized_tree = tree
  )

  class(base_tree) <- "lessTREE"

  return(base_tree)

}


#' semtree_ridge
#'
#' regularize tree with ridge penalty
#'
#' @param base_tree object returned by the regularize_tree()-function
#' @param lambdas	numeric vector: values for the tuning parameter lambda
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related
#' procedures (currently gist).
#' @param control	used to control the optimizer. This element is generated with
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet for more details.
#' @return regularized SEM-tree
#' @export
semtree_ridge <- function(base_tree,
                       lambdas,
                       method = "glmnet",
                       control = controlGlmnet()){

  base_tree$lessSEM_object <- suppressMessages(lessSEM::ridge(lavaanModel = base_tree$lavaanModel,
                                                   regularized = base_tree$regularized,
                                                   lambdas = lambdas,
                                                   modifyModel = base_tree$modifyModel,
                                                   control = control))

  return(base_tree)
}

#' semtree_lasso
#'
#' regularize tree with ridge penalty
#'
#' @param base_tree object returned by the regularize_tree()-function
#' @param lambdas	numeric vector: values for the tuning parameter lambda
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related
#' procedures (currently gist).
#' @param control	used to control the optimizer. This element is generated with
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet for more details.
#' @return regularized SEM-tree
#' @export
semtree_lasso <- function(base_tree,
                       lambdas,
                       method = "glmnet",
                       control = controlGlmnet()){
  base_tree$lessSEM_object <- suppressMessages(lessSEM::lasso(lavaanModel = base_tree$lavaanModel,
                                                   regularized = base_tree$regularized,
                                                   lambdas = lambdas,
                                                   method = method,
                                                   control = control,
                                                   modifyModel = base_tree$modifyModel))
  return(base_tree)
}

#' semtree_adaptive_lasso
#'
#' regularize tree with ridge penalty
#'
#' @param base_tree object returned by the regularize_tree()-function
#' @param lambdas	numeric vector: values for the tuning parameter lambda
#' @param weights weights used by the adaptive lasso. If set to NULL, lessTREE will
#' construct weights based on a ridge regularized model.
#' @param method which optimizer should be used? Currently implemented are ista
#' and glmnet. With ista, the control argument can be used to switch to related
#' procedures (currently gist).
#' @param control	used to control the optimizer. This element is generated with
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet for more details.
#' @return regularized SEM-tree
#' @export
semtree_adaptive_lasso <- function(base_tree,
                                lambdas,
                                weights = NULL,
                                method = "glmnet",
                                control = controlGlmnet()){

  if(is.null(weights)){
    # we can find some approximate weights using a ridge penalty
    ridge_start <- base_tree$lavaanModel |>
      ridge(regularized = base_tree$regularized,
            lambdas = .01,
            modifyModel = base_tree$modifyModel)

    start <- coef(ridge_start)@estimates[1,]

    # Define adaptive lasso weights:
    weights <- 1/abs(start)
  }

  base_tree$lessSEM_object <- suppressMessages(lessSEM::adaptiveLasso(lavaanModel = base_tree$lavaanModel,
                                                           regularized = base_tree$regularized,
                                                           lambdas = lambdas,
                                                           weights = weights,
                                                           modifyModel = base_tree$modifyModel))

  return(base_tree)
}

#' semtree_adaptive_lasso
#'
#' regularize tree with mixed penalty. See ?lessSEM::mixedPenalty for more details.
#'
#' @param base_tree object returned by the regularize_tree()-function
#' @param control	used to control the optimizer. This element is generated with
#' the controlIsta and controlGlmnet functions. See ?controlIsta and ?controlGlmnet for more details.
#' @return regularized SEM-tree
#' @export
semtree_adaptive_lasso <- function(base_tree,
                               control = controlIsta()){
  return(
    lessSEM::mixedPenalty(lavaanModel = base_tree$lavaanModel,
                          modifyModel = base_tree$modifyModel,
                          control = control)
  )
}

#' print.lessTREE
#'
#' really lazy print for lessTREE...
#' @param x object of class lessTREE
#' @param ... nothing, to be honest
#' @return nothing
#' @export
print.lessTREE <- function(x, ...){
  cat("\nNothing to see here.\n")
}

#' select_final
#'
#' select and return the final parameters from a lessTREE object
#' @param lessTREE object of class lessTREE
#' @param criterion AIC or BIC
#' @return list with summaries
#' @export
select_final <- function(lessTREE, criterion){

  # Extract final parameter values using the AIC or BIC
  if(criterion == "AIC")
    min_at <- which.min(AIC(lessTREE$lessSEM_object)$AIC)
  if(criterion == "BIC")
    min_at <- which.min(BIC(lessTREE$lessSEM_object)$BIC)

  df_coef <- as.data.frame(lessTREE$lessSEM_object@parameters[min_at,])

  # Make output ----

  parameter_average <- df_coef[c("lambda", "alpha", lessTREE$base_parameters)]

  df_delta <- df_coef |>
    pivot_longer(cols = contains("leaf"),
                 names_to = "par_leaf",
                 values_to = "delta") |>
    select(par_leaf, delta) |>
    mutate(parameter = str_replace(par_leaf, pattern = "delta_",
                                   replacement = "")) |>
    mutate(parameter = str_replace(parameter, pattern = "_leaf_",
                                   replacement = "")) |>
    mutate(parameter = str_replace(parameter, pattern = "[:digit:]",
                                   replacement = "")) |>
    mutate(leaf = parse_number(par_leaf)) |>
    select(leaf, parameter, delta)

  list(n_par = BIC(lessTREE$lessSEM_object)$nonZeroParameters[min_at],
       n_leafs = lessTREE$n_leafs,
       lambda = df_coef["lambda"],
       alpha = df_coef["alpha"],
       parameter_average = parameter_average,
       delta = as.data.frame(df_delta),
       internal = lessTREE)

}
