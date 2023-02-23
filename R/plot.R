#' plot.lessTREE
#'
#' Plot the final lessTREE model.
#' @param x lessTREE model
#' @param y not used
#' @param ... additional arguments passed to rpart
#' @export
plot.lessTREE <- function(x, y = NULL, ...){

  no.plot <- FALSE

  ## For simplicity, we will use the BIC at the moment. This should be provided
  # with an argument later on.

  # Note: The function currently only works for a single tuning parameter (lambda)

  best_lambda <- BIC(x$lessSEM_object)$lambda[(which.min(BIC(x$lessSEM_object)$BIC))]

  #### THE FOLLOWING CODE IS ADAPTED FROM THE CODE PROVIDED IN
  #### THE SEMTREE PACKAGE (SEE semtree:::plot.semtree)

  tree <- x$unregularized_tree

  final_pars <- unlist(x$lessSEM_object@parameters[
    x$lessSEM_object@parameters$lambda == best_lambda,
    x$lessSEM_object@parameterLabels])

  final_models <- lessSEM2Lavaan(regularizedSEM = x$lessSEM_object,
                                 lambda = best_lambda)

  # we want to know the node ids of our models. This will
  # allow us to exchange the parameters below.
  leafs <- getLeafs(tree)
  node_ids <- c()
  for(i in 1:length(leafs)){
    node_ids <- c(node_ids, leafs[[i]]$node_id)
  }

  to.rpart.rec <- function(x, xx, node_ids, final_models, final_pars) {
    if (is.null(xx)) {
      num <- 0
      level <- 1
    }
    else {
      num <- sum(2^(length(xx):1 - 1) * xx)
      level <- 2^length(xx)
    }
    num <- num + level
    data <- c()
    if (x$caption == "TERMINAL") {

      # we have to change the estimates here:
      coefs <- coef(final_models[[which(node_ids == x$node_id)]])

      # get current leaf:
      leaf_label <- unique(stringr::str_extract(string = names(coefs),
                                                pattern = "_leaf_[0-9]*$"))

      # remove _leaf_number:
      names(coefs) <- gsub(pattern = "_leaf_[0-9]+$",
                           replacement = "",
                           x = names(coefs))
      coefs <- coefs[unique(names(coefs))]

      # Now, we check if the corresponding deltas are 0
      deltas <- paste0("delta_",names(coefs), leaf_label)
      is_zero <- rep(FALSE, length(deltas))
      names(is_zero) <- deltas

      for(delta in deltas){
        if(delta %in% names(final_pars))
          is_zero[delta] <- final_pars[delta] == 0
      }

      row <- c("<leaf>", x$N, x$N, NA, x$node_id, 0, 0,
               0, num, paste(
                             paste0(ifelse(is_zero, "**", ""),
                                    names(coefs), " = ", round(coefs,
                                                            3),
                                   ifelse(is_zero, "**", "")), collapse = "\n"), "")
      data <- rbind(data, row)
    }
    else {
      crit <- paste("LR=", round(x$lr, 1), "(df=", x$df,
                    ")", sep = "")
      row <- c(x$caption, x$N, x$N, 0, x$node_id, 0, 0,
               0, num, NA, crit)
      data <- rbind(data, row)
      row <- to.rpart.rec(x$left_child, append(xx, 0), node_ids, final_models, final_pars)
      data <- rbind(data, row)
      row <- to.rpart.rec(x$right_child, append(xx, 1), node_ids, final_models, final_pars)
      data <- rbind(data, row)
    }
    return(data)
  }
  l <- list()
  data <- to.rpart.rec(tree, NULL, node_ids, final_models, final_pars)


  l$frame <- data
  l$frame <- data.frame(l$frame, row.names = l$frame[, 9])
  names(l$frame) <- c("var", "n", "wt", "dev", "yval", "complexity",
                      "ncompete", "nsurrogate", "label", "estimates", "crit")
  l$frame[, 1] <- as.factor(l$frame[, 1])
  l$frame$dev <- as.numeric(as.character(l$frame$dev))
  l$frame$yval <- as.numeric(as.character(l$frame$yval))
  l$frame$wt <- as.numeric(as.character(l$frame$wt))
  l$frame$complexity <- as.numeric(as.character(l$frame$complexity))
  l$frame$n <- as.numeric(as.character(l$frame$n))
  l$frame$ncompete <- as.numeric(l$frame$ncompete)
  l$frame$estimates <- as.character(l$frame$estimates)
  l$frame$nsurrogate <- as.numeric(l$frame$nsurrogate)
  l$method <- "anova"
  formatg <- function(x, digits = getOption("digits"), format = paste0("%.",
                                                                       digits, "g")) {
    if (!is.numeric(x))
      stop("'x' must be a numeric vector")
    temp <- sprintf(format, x)
    if (is.matrix(x))
      matrix(temp, nrow = nrow(x))
    else temp
  }
  l$functions$summary <- function(yval, dev, wt, ylevel, digits) {
    paste("  mean=", formatg(yval, digits), ", MSE=", formatg(dev/wt,
                                                              digits), sep = "")
  }
  l$functions$text <- function(yval, dev, wt, ylevel, digits,
                               n, use.n) {
    paste("#", yval, ", N=", n, sep = "")
  }
  class(l) <- "rpart"
  if (no.plot) {
    return(l)
  }
  else {
    rpart.plot::prp(l, left = F, type = 2, roundint = FALSE,
                    node.fun = semtree:::nodeFunSemtree,
                    varlen = 0,
                    ...)
  }
}
