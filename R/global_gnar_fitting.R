get_neighbourhood_regression <- function(stages_tensor, current_lag, t, r_stages, NTS_data, weights_matrix) {
  max_stage = r_stages[current_lag]
  Z = vapply(seq(1:max_stage), function(x) {as.numeric((weights_matrix * stages_tensor[[x]]) %*% NTS_data[t - current_lag, ])}, 
             rep(0, ncol(NTS_data)))
  return(cbind(Z))
}

get_time_step_dmat <- function(stages_tensor, p_lag, t, r_stages, NTS_data, weights_matrix) {
  for (current_lag in 1:p_lag) {
    xvec = NTS_data[t - current_lag, ]
    if (r_stages[current_lag] > 0) {
      Z = get_neighbourhood_regression(stages_tensor, current_lag, t, r_stages, NTS_data, weights_matrix)
      if (current_lag == 1) {
        dmat = cbind(xvec, Z)
      } else {
        aux = cbind(xvec, Z)
        dmat = cbind(dmat, aux)
      }
    } else {
      if (current_lag == 1) {
        dmat = cbind(xvec)
      } else {
        dmat = cbind(dmat, xvec)
      }
    }
  }
  return (as.matrix(dmat))
}

get_time_step_inner <- function(stages_tensor, current_lag, p_lag, t, r_stages, NTS_data, weights_matrix) {
  xvec = NTS_data[t - current_lag, ]
  if (r_stages[current_lag] > 0) {
    Z = get_neighbourhood_regression(stages_tensor, current_lag, t, r_stages, NTS_data, weights_matrix)
    dmat = cbind(xvec, Z)
    } else {
      dmat = cbind(xvec)
    }
  return(as.matrix(dmat))
}

get_time_step_dmat_edt <- function(stages_tensor, p_lag, t, r_stages, NTS_data, weights_matrix){
  d = ncol(NTS_data)
  dmat <- lapply(seq(1:p_lag), function(x) {get_time_step_inner(stages_tensor, x, p_lag, t, r_stages, NTS_data, weights_matrix)})
  dmat_merged <- recursive_dmat_inner(dmat)
  return(dmat_merged)
}

recursive_dmat <- function(dmat_list) {
  if (length(dmat_list) == 1) {
    return (dmat_list[[1]])
  } else {
    n_dmats = length(dmat_list)
    mid_point = n_dmats %/% 2
    left_split = recursive_dmat(dmat_list[1:mid_point])
    right_split = recursive_dmat(dmat_list[(mid_point + 1):n_dmats])
    return (rbind(left_split, right_split))
  }
}

recursive_dmat_inner <- function(dmat_list) {
  if (length(dmat_list) == 1) {
    return (dmat_list[[1]])
  } else {
    n_dmats = length(dmat_list)
    mid_point = n_dmats %/% 2
    left_split = recursive_dmat_inner(dmat_list[1:mid_point])
    right_split = recursive_dmat_inner(dmat_list[(mid_point + 1):n_dmats])
    return (cbind(left_split, right_split))
  }
}

build_linear_model <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  n_steps = nrow(NTS_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(NTS_data), 1), function(x) {return(cbind(NTS_data[x, ]))})
  yvec = recursive_dmat(response_vec)
  if (ar == 'yes') {
    design_mats = lapply(seq(1:n_steps), function(x) { get_time_step_dmat(stages_tensor, p_lag, x + p_lag, r_stages, NTS_data, weights_matrix) })
  } else {
    design_mats = lapply(seq(1:n_steps), function(x) { get_neighbourhood_regression(stages_tensor, p_lag, x + p_lag, 
                                                                                    r_stages, NTS_data, weights_matrix) })
  }
  dmat = recursive_dmat(design_mats)
  
  return (cbind(yvec, dmat))
}

build_linear_model_edt <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  n_steps = nrow(NTS_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(NTS_data), 1), function(x) {return(cbind(NTS_data[x, ]))})
  yvec = recursive_dmat(response_vec)
  if (ar == 'yes') {
    design_mats = lapply(seq(1:n_steps), function(x) { get_time_step_dmat_edt(stages_tensor, p_lag, x + p_lag, r_stages, NTS_data, weights_matrix) })
  } else {
    design_mats = lapply(seq(1:n_steps), function(x) { get_neighbourhood_regression(stages_tensor, p_lag, x + p_lag, 
                                                                                    r_stages, NTS_data, weights_matrix) })
  }
  dmat = recursive_dmat(design_mats)
  return (cbind(yvec, dmat))
}

fit_gnar <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar = 'yes') {
  model_data = build_linear_model_edt(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar)
  gnar_fit <- lm(model_data[, 1] ~ model_data[, 2:ncol(model_data)] + 0)
  # print(summary(gnar_fit))
  return (gnar_fit)
}

global_gnar_fit <- function(vts, network, max_lag, r_stages, weights_matrix, ar = "yes") {
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(r_stages))
  return(fit_gnar(stages_tensor, max_lag, r_stages, vts, weights_matrix, ar))
}

