#' Traverse Tree
#'
#' Internal function used to extract the partitioned data set.

traverse_tree <- function(nodes, current_node_ID, tree) {

  if (nodes[[current_node_ID]]$leaf) {
    return(nodes)
  } else {

    parent_node_semtree <- getNodeById(tree, id = length(nodes))
    split_rule <- paste0("nodes[[current_node_ID]]$covariates$",
                         parent_node_semtree$rule$name,
                         parent_node_semtree$rule$relation,
                         parent_node_semtree$rule$value)

    # left side ----
    left_child_semtree <- getNodeById(tree, id = length(nodes) + 1)

    nodes[[left_child_semtree$node_id]] <- list (
      node_id = left_child_semtree$node_id,
      leaf = left_child_semtree$caption == "TERMINAL",
      ID = nodes[[current_node_ID]]$ID[!eval(parse(text = split_rule))],
      model_data = nodes[[current_node_ID]]$model_data[!eval(parse(text = split_rule)), ],
      covariates = nodes[[current_node_ID]]$covariate[!eval(parse(text = split_rule)), ],
      n = left_child_semtree$N
    )

    # recursively continue splitting on the left side
    nodes <- traverse_tree(nodes = nodes, current_node_ID = length(nodes),
                           tree = tree)

    # right side ----

    right_child_semtree <- getNodeById(tree, id = length(nodes) + 1)

    nodes[[right_child_semtree$node_id]] <- list (
      node_id = right_child_semtree$node_id,
      leaf = right_child_semtree$caption == "TERMINAL",
      ID = nodes[[current_node_ID]]$ID[eval(parse(text = split_rule))],
      model_data = nodes[[current_node_ID]]$model_data[eval(parse(text = split_rule)), ],
      covariates = nodes[[current_node_ID]]$covariate[eval(parse(text = split_rule)), ],
      n = right_child_semtree$N
    )

    # recursively continue splitting on the right side
    nodes <- traverse_tree(nodes = nodes, current_node_ID = length(nodes),
                           tree = tree)
  }

  nodes

}
