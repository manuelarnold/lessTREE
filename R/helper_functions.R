#' get_parameter_labels
#'
#' @param model lavaan model
#' @param leaf integer: number of the leaf
#' @keywords internal
get_parameter_labels <- function(model, leaf){

  ### lessSEM can only regularize parameters that follow the naming convention
  # of variables in C++. That is, the standard labels used by lavaan (e.g.,
  # f=~y1) will not work. We could just change the labels, but this may be
  # inconvenient for users who want to have control over which parameters should
  # be regularized. The following approach is: check if the models has parameter
  # labels given by the user. Check these labels for consistency with C++
  # naming and then provide feedback that only those parameters with explicit
  # labels can be regularized.

  parameter_table <- model@ParTable
  parameter_labels <- parameter_table$label

  # If the user left some parameters unlabeled, we will take care of that here:

  parameter_labels[parameter_labels == ""] <- paste0(
    parameter_table$lhs[parameter_labels == ""],
    parameter_table$op[parameter_labels == ""],
    parameter_table$rhs[parameter_labels == ""])
  # this will assign the standard lavaan labels to the previously unlabeled
  # parameters

  # Next, we will check which of these labels are valid in C++:
  valid_labels <- sapply(parameter_labels, is_valid_name)

  ## Make labels unique
  parameter_table$label <- paste0(parameter_labels,
                                  paste0("_leaf_", leaf))

  # Finally, let's remove all fixed parameters:
  parameter_labels <- parameter_labels[parameter_table$free != 0]
  valid_labels <- valid_labels[parameter_table$free != 0]

  return(list(
    base_labels = parameter_labels, # labels prior to changing the names
    parameter_table = parameter_table, # includes name changes
    valid_labels = valid_labels # indicates if the name can be used in the
    # regularization
  ))

}


#' is_valid_name
#'
#' checks if a given variable name is valid in C++
#' @param x name of the variable to be checked
#' @return boolean
#' @keywords internal
is_valid_name <- function(x){

  return(
    # A valid name must start with a letter or underscore and continue
    # with letters, underscores, or numbers.
    grepl(pattern = "^[a-zA-Z_][a-zA-Z_0-9]*$",
          x = x)
  )
}
