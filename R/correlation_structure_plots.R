stage_beta <- function(beta_coeffs, r_stage) {
  out = sum(abs(beta_coeffs[, r_stage]))
  return (out)
}

active_node_scale <- function(beta_coeffs, stages_tensor, weight_matrix) {
  active_r = ncol(beta_coeffs)
  out = matrix(rep(0, nrow(weight_matrix) * ncol(weight_matrix)), nrow = nrow(weight_matrix), ncol = ncol(weight_matrix))
  aux = matrix(rep(0, nrow(weight_matrix) * ncol(weight_matrix)), nrow = nrow(weight_matrix), ncol = ncol(weight_matrix))
  for (r in 1:active_r) {
    beta = stage_beta(beta_coeffs, r)
    # out = out + 1 / r * beta * (weight_matrix * stages_tensor[[r]])
    out = out + beta * (weight_matrix * stages_tensor[[r]])
    aux = aux + stages_tensor[[r]]
  }
  out = normalize_weights(list(aux), 1, out)
  return (out)
}

global_relevance_scale <- function(node_i, node_j, weight_matrix) {
  aux_mat = weight_matrix
  # aux_mat[aux_mat > 0.0] = 1
  return (abs(sum(aux_mat[, node_i]) - sum(aux_mat[, node_j])))
}

global_relevance_scale_matrix <- function(weight_matrix) {
  out = matrix(rep(0, nrow(weight_matrix) * ncol(weight_matrix)), nrow = nrow(weight_matrix), ncol = ncol(weight_matrix))
  for (i in 1:ncol(weight_matrix)) {
    for (j in i:nrow(weight_matrix)){
      if (i == j) {
        out[i, j] = NA
      } else {
        aux = global_relevance_scale(i, j, weight_matrix)
        out[i, j] = aux
        out[j, i] = aux
      }
    }
  }
  max_val = max(out[!is.na(out)])
  out = 1 / max_val * out
  return (out)
}

get_stoplight_grid <- function(vts_dim) {
  x = c()
  y = c()
  for (i in 1:vts_dim) {
    x = c(x, rep(i, vts_dim))
    y = c(y, seq(vts_dim, 1, -1))
  }
  out = cbind(x, y)
  return (out)
}

stoplight_data <- function(top_diag, bot_diagonal) {
  out = matrix(rep(0, nrow(top_diag) * ncol(bot_diagonal)), nrow = ncol(top_diag), ncol = nrow(bot_diagonal))
  for (i in 1:ncol(top_diag)) {
    for (j in 1:nrow(bot_diagonal)) {
      if (i == j) {
        out[i, j] = NA
      } else if (i < j) {
        out[i, j] = top_diag[i, j]
      } else {
        out[i, j] = bot_diagonal[i, j]
      }
    }
  }
  return (out)
}

stoplight_data_frame <- function(top_diag, bot_diag, block_names) {
  data_mat = stoplight_data(top_diag, bot_diag)
  point_grid = get_stoplight_grid(nrow(top_diag))
  out = matrix(rep(0, length(point_grid) * 2), nrow = nrow(point_grid), ncol = 4)
  colnames(out) <- c('x', 'y', 'scale_val', 'type')
  curr_column = 1
  row = 1
  k = 1
  for (i in 1:nrow(point_grid)) {
    out[i, 3] = data_mat[row, curr_column]
    if (row <= curr_column) {
      # local relevance
      out[i, 4] = 0
    } else {
      # global relevance
      out[i, 4] = 1
    }
    if (i %% nrow(top_diag) == 0) {
      curr_column = curr_column + 1
      row = 1
    } else {
      row = row + 1
    }
  }
  out[, 1] = point_grid[, 1]
  out[, 2] = point_grid[, 2]
  out = data.frame(out)
  top_data = out[out$type == 0, ]
  bot_data = out[out$type == 1, ]
  colnames(top_data) <- c('x', 'y', block_names[1], 'type')
  colnames(bot_data) <- c('x', 'y', block_names[2], 'type')
  return (list(top_data, bot_data))
}

get_stoplight_local_relevance <- function(stoplight_grid, local_relevance_matrix, vts_dim) {
  aux = rep(0, nrow(stoplight_grid))
  for (i in 1:nrow(stoplight_grid)) {
    if (stoplight_grid[i, 1] + stoplight_grid[i, 2] == vts_dim + 1) {
      aux[i] = NA
    } else {
      node_i = stoplight_grid[i, 1] 
      node_j = vts_dim - stoplight_grid[i, 2] + 1
      aux[i] = local_relevance_matrix[node_j, node_i]
    }
  }
  out = cbind(stoplight_grid, aux)
  colnames(out) <- c("x", "y", "loc_relevance")
  return (out)
}

get_correlated_nodes <- function(r_star, adjancency_mat) {
  corr_tens = get_k_stages_adjacency_tensor(adjancency_mat, 2 * r_star)
  out = corr_tens[[1]]
  for (i in 2:(2 * r_star)) {
    if (i <= r_star) {
      out = out + 1 / i * corr_tens[[i]]
    } else {
      out = out + 1 / (2 * i) * corr_tens[[i]]
    }
  }
  return (out)
}

active_node_plot <- function(vts, network, max_lag, r_stages) {
  beta_coeff_mat = beta_coeffs(vts, network, max_lag, r_stages)
  r_star = max(r_stages)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), r_star)
  W_norm = weights_matrix(network, r_star)
  data_mat = active_node_scale(beta_coeff_mat, stages_tensor, W_norm)
  dta_plot = stoplight_data_frame(data_mat, data_mat, block_names = c('Influence', 'Influence'))
  data_for_plot = rbind(dta_plot[[1]], dta_plot[[2]])
  ggplot() +
    geom_raster(data = data_for_plot, aes(x = data_for_plot$x, y = data_for_plot$y, fill = data_for_plot$Influence)) +
    scale_fill_viridis(option='viridis') +
    theme_void() +
    labs(title = paste0("Active r-stage Neighbour Influence for a GNAR(", as.character(max_lag), " [", toString(r_stages), "])")) +
    theme(plot.title = element_text(hjust=0.5, size=30),
          legend.position = 'bottom',
          legend.key.size = unit(2, 'cm'),
          legend.title = element_blank()
    )
          # legend.title = element_text(size=30))
}

beta_coeffs <- function(vts, network, max_lag, r_stages) {
  r_star = max(r_stages)
  fit_model <- GNARfit(vts, network, alphaOrder = max_lag, betaOrder = r_stages)
  stages_tensor = get_k_stages_adjacency_tensor(as.matrix(GNARtoigraph(network)), r_star)
  out = c()
  counter = 1
  for (r in 1:length(r_stages)) {
    aux = c()
    for (sr in 1:r_stages[[r]]) {
      aux = c(aux, fit_model$mod$coefficients[[sr + counter]])
    } 
    if (r_stages[[r]] < r_star) {
      aux = c(aux, rep(0, r_star - r_stages[[r]]))
    }
    counter = counter + 1 + r_stages[[r]]
    out = rbind(out, aux)
  }
  return (out)
}

node_relevance_hierarchy <- function(network, r_star, node_names) {
  W = weights_matrix(network, r_star)
  weight_sum = colSums(W)
  relevance_scale_base = max(weight_sum)
  relevance_scale = 1 / relevance_scale_base * weight_sum
  nodes = node_names
  data_aux = data.frame(nodes, relevance_scale)
  colnames(data_aux) <- c('Node', 'Relevance')
  data_aux = data_aux[order(data_aux$Relevance), ]
  rownames(data_aux) <- NULL
  return (data_aux)
}

node_relevance_plot <- function(network, r_star, node_names, node_label_size = 2) {
  if (missing(node_names)) {
    node_names = as.character(seq(1, length(network$edges), 1))
  }
  data_plot <- node_relevance_hierarchy(network, r_star, node_names)
  relevance_plot <- ggplot(data_plot, aes(x = .data$Node, y = .data$Relevance)) + geom_col(aes(x = seq(length(node_names), 1, -1), fill = .data$Relevance)) +
    scale_fill_viridis(option = 'viridis') +
    coord_flip() +
    theme_void() +
    labs(title = paste("Node Relevance Plot for max r-stage equal to", as.character(r_star))) +
    theme(
          plot.title = element_text(hjust=0.5, size=30),
          axis.title.x = element_blank(),
          axis.text.y=element_blank(),  #remove y axis labels
          legend.position = 'bottom',
          axis.ticks.y=element_blank(),
          legend.key.size = unit(2, 'cm'),
          legend.title = element_blank()
    ) +
    annotate("text", x=seq(length(node_names), 1, -1), y=-0.01, label=data_plot$Node, col="black", size = node_label_size)
  return (list(relevance_plot, data_plot))
}

local_relevance_plot <- function(network, r_star) {
  adj_matrix = as.matrix(GNARtoigraph(network))
  correlated_nodes = get_correlated_nodes(r_star, adj_matrix)
  data_frames <- stoplight_data_frame(correlated_nodes, correlated_nodes, block_names = c('Local', 'Local'))
  data_for_plot = rbind(data_frames[[1]], data_frames[[2]])
  ggplot() +
    geom_raster(data = data_for_plot, aes(x = data_for_plot$x, y = data_for_plot$y, fill = data_for_plot$Local)) +
    scale_fill_viridis(option='viridis') +
    theme_void() +
    labs(title = paste("Conditional Correlation Heat-Map for Max r-stage equal to", as.character(r_star))) +
    theme(plot.title = element_text(hjust=0.5, size=30),
          legend.position = 'bottom', 
          legend.key.size = unit(2, 'cm'),
          legend.title = element_blank())
}

cross_correlation_plot <- function(h, vts) {
  data_mat = ccf(h, vts)
  dta_plot <- stoplight_data_frame(data_mat, data_mat, block_names = c('Cross_Corr', 'Cross_Corr'))
  data_for_plot = rbind(dta_plot[[1]], dta_plot[[2]])
  ggplot() +
    geom_raster(data = data_for_plot, aes(x = data_for_plot$x, y = data_for_plot$y, fill = data_for_plot$Cross_Corr)) +
    scale_fill_viridis(option='viridis') +
    theme_void() +
    labs(title = paste0(as.character(h), "-lag Cross-Correlation Heat-Map")) +
    theme(plot.title = element_text(hjust=0.5, size=30),
          legend.position = 'bottom',
          legend.key.size = unit(2, 'cm'),
          legend.title = element_blank())
}

ccf <-function(h, vts) {
  aux = 0.0
  xbar = colMeans(vts)
  right = c()
  for (i in 1:ncol(vts)) {
    right = c(right, sd(vts[, i]))
  }
  right = diag(1 / right)
  for (t in 1:(nrow(vts) - h)) {
    aux = aux + (vts[t + h, ] - xbar) %*% t(vts[t, ] - xbar)
  }
  aux = 1 / nrow(vts) * aux
  return (right %*% aux %*% right)
}
