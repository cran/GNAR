# Edit and fix for the PNACF
# Proper computation with missing values

nacf_edt <- function(h, s, weight_matrix, stages_tensor, nts_data) {
  numerator = 0.0
  denominator = 0.0
  if (s == 0) {
    ar_stage_matrix = diag(ncol(weight_matrix))
    weight_threshold = 0.0
  } else {
    ar_stage_matrix = (weight_matrix * stages_tensor[[s]]) + diag(ncol(stages_tensor[[s]]))
    weight_threshold = sqrt(max(colSums(as.matrix(weight_matrix * weight_matrix))))
  }
  # No missing data case
  if (length(nts_data[is.na(nts_data)]) == 0) {
    xbar = colMeans(nts_data)
    time_steps = nrow(nts_data)
    for (t in 1:time_steps){
      if (t <= (time_steps - h)) {
        numerator = numerator + t(nts_data[t + h, ] - xbar) %*% (ar_stage_matrix) %*% (nts_data[t, ] - xbar)
      }
      denominator = denominator +  t(nts_data[t, ] - xbar) %*% (nts_data[t, ] - xbar)
    }
    out = as.numeric(numerator/((1 + weight_threshold) *  denominator))
  return (out)
  # Missing data case
  } else {
    aux = nts_data
    aux[is.na(aux)] = 0.0
    xbar = colMeans(nts_data, na.rm = TRUE)
    time_steps = nrow(nts_data)
    for (t in 1:time_steps) {
        if (t <= (time_steps - h)) {
          if (is.na(sum(nts_data[t + h, ] + nts_data[t, ]))) {
            W = missing_data_adjust_weights(nts_data[t + h, ] + nts_data[t, ], weight_matrix, stages_tensor, length(stages_tensor))
            ar_stage_matrix = (W * stages_tensor[[s]]) + diag(ncol(stages_tensor[[s]]))
            w_lambda = sqrt(max(colSums(as.matrix(W * W))))
            numerator = numerator + t(aux[t + h, ] - xbar) %*% (ar_stage_matrix) %*% (aux[t, ] - xbar)
            denominator = denominator + (1 + w_lambda) * t(aux[t, ] - xbar) %*% (aux[t, ] - xbar)
          } else {
            numerator = numerator + t(nts_data[t + h, ] - xbar) %*% (ar_stage_matrix) %*% (nts_data[t, ] - xbar)
            denominator = denominator + (1 + weight_threshold) * t(nts_data[t, ] - xbar) %*% (nts_data[t, ] - xbar)
          } 
        } else {
          if (is.na(sum(nts_data[t, ]))) {
            W = missing_data_adjust_weights(nts_data[t, ], weight_matrix, stages_tensor, length(stages_tensor))
            w_lambda = sqrt(max(colSums(as.matrix(W * W))))
            denominator = denominator + (1 + w_lambda) * t(aux[t, ] - xbar) %*% (aux[t, ] - xbar)
          } else {
            denominator = denominator + (1 + weight_threshold) * t(nts_data[t, ] - xbar) %*% (nts_data[t, ] - xbar)
          }
        }
      }
    out = as.numeric(numerator / denominator)
    return(out)
  }
}

pnacf_edt <- function(h, s, weight_matrix, stages_tensor, nts_data) {
  data_net = igraphtoGNAR(graph_from_adjacency_matrix(stages_tensor[[1]], 'undirected'))
  d = ncol(nts_data)
  T_time_steps = nrow(nts_data)
  # check that there are no missing data
  if (length(nts_data[is.na(nts_data)]) == 0) {
    if (h > 1) {
        partial_model <- residuals(GNARfit(vts = nts_data, net = data_net, alphaOrder = (h - 1), betaOrder = rep((s - 1), (h - 1))))
        out = nacf_edt(h, s, weight_matrix, stages_tensor, partial_model)
    } else {
        if (s == 1) {
            out = nacf_edt(h, s,  weight_matrix, stages_tensor, nts_data)
    } else {
        partial_model <- resid_to_matrix(global_gnar_fit(nts_data, data_net, 1, c(s - 1), weight_matrix, ar = "no"), d)
        out = nacf(1, s, weight_matrix, stages_tensor, partial_model)
        }
    }
  } else {
    out = pnacf_missing_data_inner(h, s, weight_matrix, stages_tensor, nts_data, data_net, d)
  }
  return (out)
}

pnacf_missing_data_inner <- function(h, s,  weight_matrix, stages_tensor, nts_data, data_net, d) {
    W_adjsuted_list = get_missing_data_weight_matrices(nts_data, weight_matrix, stages_tensor, length(stages_tensor))
    aux = nts_data
    aux[is.na(aux)] = 0.0
    xbar = colMeans(nts_data, na.rm = TRUE)
    time_steps = nrow(nts_data)
    if (h > 1) {
        partial_model_residuals <- residuals(GNARfit(vts = nts_data, net = data_net, alphaOrder = (h - 1), betaOrder = rep((s - 1), (h - 1))))
        out = nacf(h, s,  weight_matrix, stages_tensor, partial_model_residuals)
    } else if (s == 1) {
        out = nacf(h, s,  weight_matrix, stages_tensor, nts_data)
    } else {
        partial_model <- resid_to_matrix(fit_gnar_tv_weights(stages_tensor, 1, c(s - 1), aux, W_adjsuted_list, ar = "no"), d)
        out = nacf(h, s,  weight_matrix, stages_tensor, partial_model)
    }
    return (out)
}

resid_to_matrix <- function(linear_model_fit, vts_dimension) {
    n_observations = length(linear_model_fit$residuals)
    active_time_steps = n_observations / vts_dimension
    d = vts_dimension
    lm_residuals = linear_model_fit$residuals
    residual_matrix = t(vapply(seq(1:active_time_steps), function(x) {lm_residuals[(d * (x - 1) +1):(x * d)]}, rep(0, d)))
    colnames(residual_matrix) <- vapply(seq(1:d), function(x) {paste0("u", as.character(x))}, "")
    return (residual_matrix)
}

build_linear_model_missing_weights <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix_list, ar = 'yes') {
  n_steps = nrow(NTS_data) - p_lag
  response_vec = lapply(seq((p_lag + 1), nrow(NTS_data), 1), function(x) {return(cbind(NTS_data[x, ]))})
  yvec = recursive_dmat(response_vec)
  if (ar == 'yes') {
    design_mats = lapply(seq(1:n_steps), function(x) { get_time_step_dmat_edt(stages_tensor, p_lag, x + p_lag, r_stages, NTS_data, weights_matrix_list[[x + p_lag]]) })
  } else {
    design_mats = lapply(seq(1:n_steps), function(x) { get_neighbourhood_regression(stages_tensor, p_lag, x + p_lag, 
                                                                                    r_stages, NTS_data, weights_matrix_list[[x + p_lag]]) })
  }
  dmat = recursive_dmat(design_mats)
  return (cbind(yvec, dmat))
}

fit_gnar_tv_weights <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix_list, ar = 'yes') {
  model_data = build_linear_model_missing_weights(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix_list, ar)
  gnar_fit <- lm(model_data[, 1] ~ model_data[, 2:ncol(model_data)] + 0)
  # print(summary(gnar_fit))
  return (gnar_fit)
}

global_gnar_fit_tv_weights <- function(vts, network, max_lag, r_stages, weights_matrix_list, ar = "yes") {
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), max(r_stages))
  return(fit_gnar_tv_weights(stages_tensor, max_lag, r_stages, vts, weights_matrix_list, ar))
}

missing_data_warning <- function(missing_data = TRUE) {
  if (missing_data) {
    warning('Missing values detected, continuing by adjusting the weights accordingly.')
  }
}

all_missing_stop_step <- function(vts, time_step) {
  d = ncol(vts)
  if (d == sum(is.na(vts[time_step, ]))) {
    stop('No realisations for this time-step')
  }
}

#underdetermined_linear_model_stop <- function() {
#  d = ncol(vts)
#  q = (h - 1) * s
#  n = 
#}
