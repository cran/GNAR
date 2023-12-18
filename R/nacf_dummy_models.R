get_y_z_vector <- function(stages_tensor, h_lag, r_stage, NTS_data, weights_tensor) {
  y = cbind(NTS_data[1 + h_lag, ])
  Z = cbind(NTS_data[h_lag, ], (weights_tensor * stages_tensor[[r_stage]]) %*% NTS_data[h_lag, ])
  for (t  in 2:(nrow(NTS_data) - h_lag)) {
    y = rbind(y, cbind(NTS_data[t + h_lag, ]))
    aux = cbind(NTS_data[t + h_lag - 1, ], (weights_tensor * stages_tensor[[r_stage]]) %*% NTS_data[t + h_lag - 1, ])
    Z = rbind(Z, aux)
  }
  return (cbind(y, Z))
}


get_time_dmat <- function(stages_tensor, p_lag, t, r_stages, NTS_data, weights_tensor, ar) {
  if (ar == 'yes') {
  current_lag = 1
  xvec = cbind(NTS_data[t - current_lag, ])
  Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - current_lag, ]
  if (r_stages[current_lag] > 1) {
    for (r in 2:r_stages[current_lag]) {
      Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - current_lag, ]
      Z = cbind(Z, Z_r)
    }
  }
  dmat = cbind(xvec, Z)
  if (p_lag > 1) {
  for (current_lag in 2:p_lag) {
    xvec = cbind(NTS_data[t - current_lag, ])
    Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - current_lag, ]
    if (r_stages[current_lag] > 1) {
        for (r in 2:r_stages[current_lag]) {
          Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - current_lag, ]
          Z = cbind(Z, Z_r)
      }
    }
      aux = cbind(xvec, Z)
      dmat = cbind(dmat, aux)
  }
  }
  } else {
    current_lag = 1
    Z = (weights_tensor * stages_tensor[[1]]) %*% NTS_data[t - 1, ]
    if (r_stages[current_lag] > 1) {
      for (r in 2:r_stages[current_lag]) {
        Z_r = (weights_tensor * stages_tensor[[r]]) %*% NTS_data[t - 1, ]
        Z = cbind(Z, Z_r)
      }
    }
    dmat = Z
  }
  return (dmat)
}


get_sample_dmat <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar) {
  y = cbind(NTS_data[1 + p_lag, ])
  Z = get_time_dmat(stages_tensor, p_lag, 1 + p_lag, r_stages, NTS_data, weights_matrix, ar)
  for (t  in 2:(nrow(NTS_data) - p_lag)) {
    y = rbind(y, cbind(NTS_data[t + p_lag, ]))
    aux = get_time_dmat(stages_tensor, p_lag, t + p_lag, r_stages, NTS_data, weights_matrix, ar)
    Z = rbind(Z, aux)
  }
  return (as.matrix(cbind(y, Z)))
}
  

ls_gnar <- function(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar='yes') {
  dummy_data = get_sample_dmat(stages_tensor, p_lag, r_stages, NTS_data, weights_matrix, ar)
  y = dummy_data[, 1]
  dmat = dummy_data[, 2:ncol(dummy_data)]
  out = lm(y ~ dmat + 0)
  # print(summary(out))
  return(out)
}


ls_gnar_residuals <- function(gnar_lm, vts_dim) {
  realizations =  length(residuals(gnar_lm))
  residual_values = residuals(gnar_lm)
  out = matrix(rep(0, realizations), nrow = realizations / vts_dim, ncol = vts_dim)
  k = 1
  for (i in 1:nrow(out)) {
    for (j in 1:ncol(out)) {
      out[i, j] = residual_values[k]
      k = k + 1
    }
  }
  return (out)
}


get_basic_model_summary <- function(stages_tensor, h_lag, r_stage, NTS_data, weights_tensor, intercept_option="no") {
  dummy_data = get_y_z_vector(stages_tensor, h_lag, r_stage, NTS_data, weights_tensor)
  if (intercept_option == "yes") {
    out = summary(lm(dummy_data[, 1] ~ dummy_data[, 2] + dummy_data[, 3] ))$r.squared
  } else {
    out = summary(lm(dummy_data[, 1] ~ dummy_data[, 2] + dummy_data[, 3] + 0))$r.squared
  }
  return (out)
}


get_model_summaries <- function(stages_tensor, max_lag, max_stage, NTS_data, weights_tensor, intercept_option="no") {
  mod_diag_vals = matrix(rep(0, max_lag * max_stage), nrow=max_lag, ncol=max_stage)
  for (h in 1:max_lag) {
    for (s in 1:max_stage) {
      mod_diag_vals[h, s] = get_basic_model_summary(stages_tensor, h, s, NTS_data, weights_tensor, intercept_option)
    }
  }
  return (mod_diag_vals)
 }


get_nacf_values <- function(max_lag, max_stage, weight_matrix, stages_tensor, nts_data) {
  nacf_corr_vals = matrix(rep(0, max_lag * max_stage), nrow=max_lag, ncol=max_stage)
  for (s in 1:max_stage) {
    for (h in 1:max_lag) {
      nacf_corr_vals[h, s] = nacf(h, s, weight_matrix, stages_tensor, nts_data)
    }
  }
  return (nacf_corr_vals)
}


get_pnacf_values <- function(max_lag, max_stage, weight_matrix, stages_tensor, nts_data) {
  nacf_corr_vals = matrix(rep(0, max_lag * max_stage), nrow=max_lag, ncol=max_stage)
  for (s in 1:max_stage) {
    for (h in 1:max_lag) {
      nacf_corr_vals[h, s] = pnacf(h, s, weight_matrix, stages_tensor, nts_data)
    }
  }
  return (nacf_corr_vals)
}


get_planet_coordinates <- function(max_lag, max_stage) {
  y = c()
  x = c()
  for (stage in 1:max_stage) {
    dummy_radius = stage
    aux = get_circle_points_polar(dummy_radius, max_lag)
    y = c(y, aux[, 2])
    x = c(x, aux[, 1])
    }
  return (cbind(x, y))
}


get_dummy_circle_points <- function(radius, num_points) {
  k = 0
  if ((num_points %% 2 ) != 0) {
    num_points = num_points + 1
    k = 1
  }
  aux_points = round((num_points + 2) / 2)
  x = seq(-radius, radius, (2 * radius) / (aux_points - 1))
  x = -sort(-x)
  xtop = c()
  xbot = c()
  ytop = c()
  ybot = c()
  for (i in 1:length(x)) {
    aux0 = sqrt(radius*radius - x[i] * x[i])
    aux1 = -1 * aux0
    if (aux0 == aux1) {
      xtop = c(xtop, x[i])
      ytop = c(ytop, aux0)
    } else {
      xtop = c(xtop, x[i])
      xbot = c(x[i], xbot)
      ytop = c(ytop, aux0)
      ybot = c(aux1, ybot)
    }
  }
  xout = c(xtop, xbot)
  y = c(ytop, ybot)
  return (cbind(xout, y)[1:(length(xout) - k), ])
}


get_circle_points_polar <- function(radius, num_points) {
  x = rep(0, num_points)
  y = rep(0, num_points)
  if (num_points == 2) {
    x[1] = radius
    y[1] = 0.0
    x[2] = - radius
    y[2] = 0.0
  } else {
  for (i in 1:num_points) {
    theta_val = (i - 1) * 2 * pi / num_points
    # print(theta_val)
    x[i] = radius * cos(theta_val)
    y[i] = radius * sin(theta_val)
  }
  }
  return (cbind(x, y))
}


get_m_circles <-function(centre_coords, m_num_points, max_stage) {
  x = rep(0, length(centre_coords[, 1]) * (m_num_points + 1))
  y = rep(0, length(centre_coords[, 2]) * (m_num_points + 1))
  max_lag = length(centre_coords[, 1]) / max_stage
  base_radi = (1 / 4) * 2 * pi / max_lag
  base_cirlce = get_circle_points_polar(base_radi, m_num_points)
  k = 1
  r_aux = 1
  for (i in 1:length(centre_coords[, 1])) {
    x[k] = centre_coords[i, 1]
    y[k] = centre_coords[i, 2]
    k = k + 1
    for (j in 1:m_num_points) {
      x[k] = base_cirlce[j, 1] + centre_coords[i, 1]
      y[k] = base_cirlce[j, 2] + centre_coords[i, 2]
      k = k + 1
    }
  }
  return(cbind(x, y))
}

