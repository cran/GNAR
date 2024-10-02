# Missing data weight adjustment and user-specified weight matrix calculation

adjust_weights_missing_data <- function(vts, weight_matrix, stages_tensor, max_stage) {
  W_adjusted = vapply(seq(1:length(vts)), function(x) {adjust_weight_inner(vts, weight_matrix, x)}, rep(0, length(vts)))
  out = normalise_weights(stages_tensor = stages_tensor, max_stage = max_stage, weight_matrix = W_adjusted)
  return (out)
}

normalise_weights <- function(stages_tensor, max_stage, weight_matrix) {
  out = matrix(rep(0, ncol(stages_tensor[[1]])^2), ncol=ncol(stages_tensor[[1]]), nrow=ncol(stages_tensor[[1]]))
  out <- lapply(seq(1:max_stage), function(x) {normalise_weights_inner(stages_tensor, x, as.matrix(stages_tensor[[x]] * weight_matrix))})
  out_sum = recursive_list_sum(out)
  return (out_sum)
}

normalise_weights_inner <- function(stages_tensor, r_stage, weight_matrix_unormalised) {
  aux0 = rowSums(weight_matrix_unormalised)
  aux0[aux0 == 0] = 1
  aux = (1/aux0) * weight_matrix_unormalised
  return (aux)
}

recursive_list_sum <- function(dmat_list) {
  if (length(dmat_list) == 1) {
    return (dmat_list[[1]])
  } else {
    n_dmats = length(dmat_list)
    mid_point = n_dmats %/% 2
    left_split = dmat_list[1:mid_point]
    right_split = dmat_list[(mid_point + 1):n_dmats]
    return (left_split[[1]] + right_split[[1]])
  }
}

adjust_weight_inner <- function(vts, weight_matrix, ith_node) {
  if (is.na(vts[ith_node])) {
    out = rep(0, length(vts))
  } else {
    out = weight_matrix[, ith_node]
  }
  return (out)
}

missing_data_time_step <- function(vts, time_step) {
  if (sum(is.na(vts[time_step, ])) == 0) {
    out = FALSE
  } else {
    out = TRUE
  }
  return (out)
}

get_missing_data_weight_matrices <- function(vts, W, stages_tensor, max_stage) {
  time_steps = nrow(vts)
  W_list <- lapply(seq(1:time_steps), function(x) {adjust_weights_missing_data(vts[x, ], W, stages_tensor, max_stage)})
  return(W_list)
}

adjust_vts_missing_data <- function(vts, W, stages_tensor, max_stage) {
  if (sum(is.na(vts)) == 0) {
    cat("No missing data, weight adjustment is not necessary.")
    vts_edt = vts
  } else {
    cat("Missing values detected, adjusting weight matrices.")
    vts_edt = vts
    vts_edt[is.na(vts_edt)] = 0.0
    W_list = get_missing_data_weight_matrices(vts, W, stages_tensor, max_stage)
  }
  return(list(vts_edt, W_list))
}


