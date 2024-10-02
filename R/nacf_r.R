nacf <- function(h, s, weight_matrix, stages_tensor, nts_data) {
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


pnacf <- function(h, s,  weight_matrix, stages_tensor, nts_data) {
  data_net = igraphtoGNAR(graph_from_adjacency_matrix(stages_tensor[[1]], 'undirected'))
  if (h > 1) {
    partial_model <- residuals(GNARfit(vts = nts_data, net = data_net, alphaOrder = (h - 1), betaOrder = rep((s - 1), (h - 1))))
    out = nacf(h, s,  weight_matrix, stages_tensor, partial_model)
  } else {
    if (s == 1) {
        out = nacf(h, s,  weight_matrix, stages_tensor, nts_data)
    } else {
      # partial_model <- residuals(GNARfit(vts = nts_data, net = data_net, alphaOrder = 1, betaOrder = rep((s - 1), 1)))
      partial_model <- ls_gnar_residuals((ls_gnar(stages_tensor, 1, c(s - 1), nts_data, weight_matrix, ar='no')), ncol(nts_data))
      out = nacf(h, s,  weight_matrix, stages_tensor, partial_model)
    }
  }
  return (out)
}


r_stage_pair_selector <- function(r, t, s, nts_data, stage_adjacency) {
  y = c()
  z = c()
  aux = stage_adjacency
  for (p in 1:ncol(nts_data)) {
    for (q in p:ncol(nts_data)) {
      if (aux[p, q]  == 1) {
        y = c(y, nts_data[t, p])
        z = c(z, nts_data[s, q])
        y = c(y, nts_data[t, q])
        z = c(z, nts_data[s, p])
      }
    }
  }
  return (matrix(c(y, z), ncol=2, nrow=length(y)))
}


r_stage_lag_h_pair_builder <- function(r, h, stages_tensor, nts_data) {
  merge_by_rows = r_stage_pair_selector(r, 1 + h, 1, nts_data, stages_tensor[[r]])
  merge_by_columns = merge_by_rows
  for (t in 1:nrow(nts_data) - h) {
    aux = r_stage_pair_selector(r, t + h, t, nts_data, stages_tensor[[r]])
    merge_by_columns = cbind(merge_by_columns, aux)
    merge_by_rows = rbind(merge_by_rows, aux)
  }
  return (list(merge_by_rows, merge_by_columns))
}


get_r_stage_adjacency_matrix <- function(St_r_1, St_1, connection_matrix) {
  out = matrix(rep(0, nrow(St_1)^2), nrow=nrow(St_1), ncol=ncol(St_1))
  out_connection = connection_matrix
  for (p in 1:(ncol(St_1)-1)) {
    for (q in (p+1):ncol(St_1)) {
      if (connection_matrix[p, q] == 0) {
        if (get_r_stage_neighbor(St_r_1, St_1, p, q) == 1) {
          out[p, q] = 1
          out[q, p] = 1
          out_connection[p, q] = 1
          out_connection[q, p] = 1
        }
      }
    }
  }
  return (list(out, out_connection))
}


get_r_stage_neighbor <- function(St_r_1, St_1, p, q) {
  out = 0
  for (l in 1:ncol(St_1)) {
    if (St_1[q, l] == 1) {
      if (St_r_1[p, l] == 1 && St_r_1[p, q] == 0) {
        out = 1
        break
      }
    }
  }
  return (out)
}


get_k_stages_adjacency_tensor <- function(St_1, r) {
  out = list(St_1)
  connection_matrix = St_1
  if (r == 1) {
    return(list(connection_matrix))
  } else {
  for (s in 2:r) {
    aux = get_r_stage_adjacency_matrix(out[[s - 1]], St_1, connection_matrix)
    St_s = aux[[1]]
    connection_matrix = aux[[2]]
    out = c(out, list(St_s))
  }
  return (out)
  }
}


plt_mean <- function(x, size) {
  y = x[1:size]
  for(i in 1:size){
    y[i] = mean(x[1:i])
  }
  plot(seq(1, length(y), 1), y, 'l')
}


normalize_weights <- function(stages_tensor, max_stage, weight_matrix) {
  out = matrix(rep(0, ncol(stages_tensor[[1]])^2), ncol=ncol(stages_tensor[[1]]), nrow=ncol(stages_tensor[[1]]))
  for (r in 1:max_stage) {
    aux = as.matrix(stages_tensor[[r]] * weight_matrix)
    for (p in 1:nrow(aux)) {
      aux0 = rowSums(aux)
      aux0[aux0 == 0] = 1
      aux = (1/aux0) * aux
    }
    out = out + aux
  }
  return (out)
}


latent_emp_covariance <- function(vts) {
  xmean = colMeans(vts)
  cov_mat = (vts[2, ] - xmean) %*% t(vts[1, ] - xmean)
  for (t in 2:nrow(vts) - 1) {
    cov_mat = cov_mat + (vts[t + 1, ] - xmean) %*% t(vts[t, ] - xmean)
  }
  out = (1 / nrow(vts) * cov_mat)
  # * (matrix(rep(1, ncol(vts) * ncol(vts)), nrow= ncol(vts), ncol = ncol(vts)) - diag(ncol(vts))) +
  diag(ncol(vts))
  return(out)
}


weights_matrix <- function(network, max_r_stage) {
  W = distances(GNARtoigraph(network))
  if (!is.finite(max(W))) { 
    cat("Warning the graph is not fully connected, adjusting by removing non-connected nodes from r-stage adjacency sets.")
    W[!is.finite(W)] = 0
    }
  adj_matrix = as.matrix(GNARtoigraph(network))
  if (missing(max_r_stage)) {
    dummy_distances = distances(graph_from_adjacency_matrix(adj_matrix))
    if (!is.finite(max(dummy_distances))) {
      dummy_distances[!is.finite(dummy_distances)] = 0
    }
    max_r_stage = max(dummy_distances)
  }
  W_norm = normalize_weights(get_k_stages_adjacency_tensor(adj_matrix, max_r_stage), max_r_stage, W)
  return (W_norm)
}


missing_data_adjust_weights <- function(vts, weight_matrix, stages_tensor, max_stage) {
  for (j in 1:ncol(weight_matrix)) {
    if (is.na(vts[j])) {
      weight_matrix[, j] = 0.0
    }
  }
  out = normalize_weights(stages_tensor = stages_tensor, max_stage = max_stage, weight_matrix = weight_matrix)
  return (out)
}
