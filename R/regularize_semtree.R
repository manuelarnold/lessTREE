#' Regularize SEM Trees
#'
#' Regularize SEM trees fitted with semtree. Requires tidyverse and lessSEM.
#'
#' @param tree A semtree object. Only SEM trees fitted with lavaan work.
#' @param AIC blabla
#' @return A list.
#' @export

regularize_semtree <- function(tree, penalty = "adaptiveLasso",
                               lambdas = seq(0, 0.25,length.out = 50),
                               criterion = "AIC") {

  # Extracting the models and their parameters ----

  leafs <- getLeafs(tree)
  n_leafs <- length(leafs)

  base_parameters <- c()
  parameters <- c()
  models <- c()

  for(leaf in 1:n_leafs){

    model <- leafs[[leaf]]$model

    base_parameters <- c(base_parameters, names(getLavaanParameters(model)))

    parameter_table <- model@ParTable

    # we have to give labels to those parameters that currently don't have names:
    parameter_labels <- parameter_table$label
    parameter_table$label[parameter_labels == ""] <- paste0(
      parameter_table$lhs[parameter_labels == ""],
      parameter_table$op[parameter_labels == ""],
      parameter_table$rhs[parameter_labels == ""])

    parameter_table$label <- paste0(parameter_table$label, "_leaf_", leaf)

    models <- c(models,
                lavaan(model = parameter_table,
                       data = lavInspect(object = model, what = "data"))
    )

    parameters <- c(parameters,
                    names(getLavaanParameters(models[[leaf]])))

  }


  # Defining the transformations ----

  delta_parameters <- paste0("delta_", parameters)

  transformation <- paste0(
    "parameters: ",
    paste0(unique(c(base_parameters, delta_parameters, parameters)),
           collapse = ", "),
    "\n\n")

  for(i in 1:length(delta_parameters)){
    transformation <- paste0(
      transformation,
      paste0(parameters[i], " = ", base_parameters[i], " + ", delta_parameters[i],
             ";"),
      "\n")
  }


  # Choose penalty ----

  if (penalty == "adaptiveLasso") {

    ridge_start <- models |>
      ridge(regularized = delta_parameters,
            lambdas = .01,
            modifyModel = modifyModel(transformations = transformation))

    start <- coef(ridge_start)@estimates[1,]

    # Define adaptive lasso weights:
    weights <- 1/abs(start)

  lasso_fit <- models |>
    adaptiveLasso(regularized = delta_parameters,
                  lambdas = lambdas,
                  weights = weights,
                  modifyModel = modifyModel(transformations = transformation))
  }

  if (penalty == "lasso") {

    lasso_fit <- models |>
      lasso(regularized = delta_parameters,
                    lambdas = lambdas,
                    modifyModel = modifyModel(transformations = transformation))
  }

  # Extract final parameter values using the AIC
  min_at <- which.min(AIC(lasso_fit)$AIC)
  df_coef <- as.data.frame(lasso_fit@parameters[min_at,])


  # Make output ----

  parameter_labels <- unique(parameter_labels)
  parameter_labels <- parameter_labels[!parameter_labels == ""]

  parameter_average <- df_coef[c("lambda", "alpha", parameter_labels)]

  df_delta <- df_coef |> pivot_longer(cols = contains("leaf"),
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

  list(n_par = length(parameter_labels),
       n_leafs = n_leafs,
       lambda = df_coef["lambda"],
       alpha = df_coef["alpha"],
       parameter_average = parameter_average,
       delta = as.data.frame(df_delta),
       internal = lasso_fit)

}
