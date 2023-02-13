#' Get Partitioned Sample
#'
#' The function extracts the subsamples from the associated nodes. This function
#' is a workaround for lavaan models because covariates and ID variables are
#' currently not stored for lavaan models.
#'
#' @param tree A semtree object. Only SEM trees fitted with lavaan work.
#' @param covariates A data frame with covariates and ID variables.
#' @param leafs logical value; indicating whether only data from the leaf nodes
#' is to be extracted.
#' @return A list.
#' @export

get_partition <- function(tree, covariates, leafs = FALSE) {

  d <- as.data.frame(lavInspect(object = tree$model, what = "data"))
  covariates <- as.data.frame(covariates)
  ID <- 1:NROW(d)

  # Initialize list with nodes
  nodes <- list(list(
    node_id = 1,
    leaf = tree$caption == "TERMINAL",
    ID = ID,
    model_data = d,
    covariates = covariates,
    n = tree$N
  ))

  nodes <- traverse_tree(nodes = nodes, current_node_ID = 1, tree = tree)

  if (leafs) {
    nodes <- nodes[sapply(nodes, FUN = function(x) x$leaf)]
  }

  nodes

}
