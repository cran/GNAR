#CorrBit plot code
#Works for nacf values and/or fitted models and/or significance
#The plan is to exploit the aes capabilities of ggplot to not have to color the graph by ourselves, just need to add the necessary x y 
#coordinates for the circles and the scatter plot takes care of the rest
#The correlation is the color aesthetic, blue = -1 and orange = +1
#The R2 or another model diagnostic will correspond to the size of the circle


plot_corbit <-function(corbit_data, max_lag, max_stage, viridis_color_option="viridis", size_option="absolute_val", wagner="no") {
  corbit_data = rbind(c(0, 0, 0, 0), corbit_data)
  if (wagner == "yes") {
    label_points = get_circle_points_polar(2 * max_stage + 1, max_lag)
    num_samples = (length(corbit_data[, 1]) - 1) / (max_lag * max_stage) - 1
    capt_text = paste0("Wagner plot with ", num_samples, " time frames or samples, max lag")
  } else {
    label_points = get_circle_points_polar(max_stage + 1, max_lag)
    capt_text = "Corbit plot with max lag"
  }
  #corbit_data = rbind(c(0, -1/3, -1, 0.1), corbit_data)
  aux_size = 4
  if (size_option == "absolute_val") {
    aux_size = 3
    colnames(corbit_data) <- c("x", "y", colnames(corbit_data)[3], paste0("abs(", colnames(corbit_data)[3], ")"))
  }
  p <- ggplot(corbit_data, aes(x = .data$x, y= .data$y))
  p + geom_point(aes(color = .data[[colnames(corbit_data)[3]]], size = abs(.data[[colnames(corbit_data)[3]]]))) +
        scale_color_viridis(option=viridis_color_option) +
        # geom_circle(aes(x0=0, y0=0, r = 1), linetype='dashed', color='blue',
        #            fill='yellow', lwd=1, inherit.aes=FALSE) +
      coord_fixed() + 
      theme_void() +
      theme(axis.text.x=element_blank(), #remove x axis labels
              axis.ticks.x=element_blank(), #remove x axis ticks
              axis.text.y=element_blank(),  #remove y axis labels
              axis.ticks.y=element_blank(), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(), 
              plot.caption = element_text(hjust = 0.5, size = 14, face = "bold")
        ) + 
      annotate("text", x=label_points[, 1], y=label_points[, 2], label=rownames(data.frame(label_points)), col="black", size = 5) +
      labs(caption = paste(colnames(corbit_data)[3], capt_text, max_lag, "and max path length", max_stage, sep = " " ), 
             color = colnames(corbit_data)[3],
             size = colnames(corbit_data)[4])
}


get_corbit_data <-function(max_lag, max_stage,  W_norm, stages_tensor, nts_data, size_option = "absolute_val", partial = "no") {
  if (partial == "yes") {
    nacf_corbit = c(get_pnacf_values(max_lag, max_stage, W_norm, stages_tensor, nts_data))
    if (size_option != "absolute_val") {
      r2_corbit = c(get_model_summaries(stages_tensor, max_lag, max_stage, nts_data, W_norm, intercept_option="no"))
      model_diagnostics = matrix(c(nacf_corbit, r2_corbit), nrow=length(nacf_corbit), ncol = 2)
    } else {
      model_diagnostics = matrix(c(nacf_corbit, rep(0, length(nacf_corbit))), nrow=length(nacf_corbit), ncol = 2)
    }
    corbit_data = data.frame(cbind(get_planet_coordinates(max_lag, max_stage), model_diagnostics))
    colnames(corbit_data) <- c("x", "y", "pnacf", "R2")
  } else {
    nacf_corbit = c(get_nacf_values(max_lag, max_stage, W_norm, stages_tensor, nts_data))
    if (size_option != "absolute_val") {
      r2_corbit = c(get_model_summaries(stages_tensor, max_lag, max_stage, nts_data, W_norm, intercept_option="no"))
      model_diagnostics = matrix(c(nacf_corbit, r2_corbit), nrow=length(nacf_corbit), ncol = 2)
    } else {
      model_diagnostics = matrix(c(nacf_corbit, rep(0, length(nacf_corbit))), nrow=length(nacf_corbit), ncol = 2)
    }
    corbit_data = data.frame(cbind(get_planet_coordinates(max_lag, max_stage), model_diagnostics))
    colnames(corbit_data) <- c("x", "y", "nacf", "R2")
  }
  return (corbit_data)
}
  

get_list_mean <- function(data_list, target_index, target_column) {
  aux = 0.0
  for (i in 1:length(data_list)){
    aux = data_list[[i]][target_index, target_column] + aux
  }
  return (aux / length(data_list))
}


get_wagner_data <- function(max_lag, max_stage, W_list, stages_list, nts_samples, same_net="yes", partial = "no", centre_option="mean") {
  nacf_samples = list()
  if (same_net == "yes") {
  for (i in 1:length(nts_samples)) {
    nacf_samples <- c(nacf_samples, list(get_corbit_data(max_lag, max_stage, W_list[[1]], stages_list[[1]], 
                                                         nts_samples[[i]], partial = partial)))
  } } else {
    for (i in 1:length(nts_samples)) {
      nacf_samples <- c(nacf_samples, list(get_corbit_data(max_lag, max_stage, W_list[[i]], stages_list[[i]], 
                                                           nts_samples[[i]], partial = partial)))
    }
  }
  centre_coordinates = nacf_samples[[1]][, 1:2]
  point_coordinates = get_m_circles(centre_coordinates, length(nts_samples), max_stage)
  nacf_values = rep(0, length(point_coordinates[, 1]))
  k = 1
  aux_counter = 1
  for (stage in 1:max_stage) {
    for (lag in 1:max_lag) {
      nacf_values[k] = get_list_mean(nacf_samples, aux_counter, 3)
      k = k + 1
        for(j in 1:length(nacf_samples)) {
          nacf_values[k] = nacf_samples[[j]][aux_counter, 3]
          k = k + 1
        }
        aux_counter = aux_counter + 1
        }
  }
  model_diagnostics = matrix(c(nacf_values, abs(nacf_values)), nrow=length(nacf_values), ncol = 2)
  wagner_data = data.frame(cbind(point_coordinates, model_diagnostics))
  if (partial == "yes") {
    colnames(wagner_data) <- c("x", "y", "pnacf", "R2")
  } else {
    colnames(wagner_data) <- c("x", "y", "nacf", "R2")
  }
  wagner_data$x = 2 * wagner_data$x
  wagner_data$y = 2 * wagner_data$y
  return (wagner_data)
}


corbit_plot <- function(vts, net, max_lag, max_stage, weight_matrix, viridis_color_option="viridis", size_option="absolute_val", partial="no", wagner="no") {
  dummy_net = GNARtoigraph(net)
  adj_mat = as.matrix(dummy_net)
  net_result = graph_from_adjacency_matrix(adj_mat, 'undirected')
  stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
  if (missing(weight_matrix)) {
    W_norm = weights_matrix(net)
  } else {
    if (!is.finite(max(distances(weight_matrix)))) { 
      warning("Warning the graph is not fully connected, adjusting by removing non-connected nodes from r-stage adjacency sets.\n")
      weight_matrix[!is.finite(weight_matrix)] = 0
    }
    W_norm = normalize_weights(stages_tensor, max_stage, weight_matrix)
  }
  corbit_data = get_corbit_data(max_lag, max_stage, W_norm, stages_tensor, vts, size_option, partial = partial)
  plot_corbit(corbit_data, max_lag, max_stage,  viridis_color_option, size_option, wagner=wagner)
}


get_weight_matrices <- function(network_list, max_stage, weight_matrices, same_net="no") {
  # Same net case is for time frames and/or subsamples
  if (same_net == "yes") {
    out = list()
    dummy_net = GNARtoigraph(network_list[[1]])
    adj_mat = as.matrix(dummy_net)
    net_result = graph_from_adjacency_matrix(adj_mat, 'undirected')
    stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
    if (missing(weight_matrices)) {
      weight_matrix = distances(net_result)
      W_norm = normalize_weights(stages_tensor, max_stage, weight_matrix)
    } else {
      W_norm = normalize_weights(stages_tensor, max_stage, weight_matrices[[1]])
    }
    out = list(stages_tensor, W_norm)
  } else {
    out = list()
    wlist = list()
    stages_list = list()
    for (i in 1:length(network_list)) {
      dummy_net = GNARtoigraph(network_list[[i]])
      adj_mat = as.matrix(dummy_net)
      net_result = graph_from_adjacency_matrix(adj_mat, 'undirected')
      stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
      if (missing(weight_matrices)) {
        weight_matrix = distances(net_result)
        W_norm = normalize_weights(stages_tensor, max_stage, weight_matrix)
      } else {
        W_norm = normalize_weights(stages_tensor, max_stage, weight_matrices[[i]])
      }
      stages_list = c(stages_list, list(stages_tensor))
      wlist = c(wlist, list(W_norm))
    }
    out = list(stages_list, wlist)
  }
  return (out)
}


wagner_plot <- function(vts_frames, network_list, max_lag, max_stage, weight_matrices, frame_names, same_net="no", viridis_color_option="viridis", size_option="absolute_val", partial="no", wagner="yes") {
  if (missing(frame_names)) {
    cat("No name for each covariate and/or time slice, using generic names")
    frame_names = as.character(seq(1:length(vts_frames)))
  }
  if (length(frame_names) != length(vts_frames)) {
    cat("Missing name for slice and/or covariate, using generic names")
    frame_names = as.character(seq(1:length(vts_frames)))
  }
  wagner_stages_weights = get_weight_matrices(network_list = network_list, max_stage = max_stage, weight_matrices = weight_matrices, same_net=same_net)
  stages_tensor_list = wagner_stages_weights[[1]]
  weight_matrices_list = wagner_stages_weights[[2]]
  if (same_net == "yes") {
    stages_tensor_list = list()
    weight_matrices_list = list()
    for (i in 1:length(vts_frames)) {
      stages_tensor_list = c(stages_tensor_list, list(wagner_stages_weights[[1]]))
      weight_matrices_list = c(weight_matrices_list, list(wagner_stages_weights[[2]]))
    }
  }
  wagner_plot_data = get_wagner_data(max_lag, max_stage, weight_matrices_list, stages_tensor_list, vts_frames, same_net=same_net, partial = partial)
  ggarrange(plot_corbit(wagner_plot_data, max_lag, max_stage,  viridis_color_option, size_option, wagner=wagner), 
            plot_wagner_legend(frame_names, partial=partial), widths = c(2, 0.5),
            ncol = 2, nrow = 1)
}


get_covariate_weight_matrix <- function(weight_matrix, covariate_criteria=NULL, covariate_column, vts) {
  covariate_matrix = c()
  for (i in 1:nrow(vts)) {
    if (vts[i, covariate_column] == covariate_criteria) {
      aux = rep(1, nrow(vts))
    } else {
      aux = rep(0, nrow(vts))
    }
    covariate_matrix = cbind(covariate_matrix, aux)
  }
  colnames(covariate_matrix) = vts[, 1]
  return (weight_matrix * covariate_matrix)
}


plot_wagner_legend <- function(frame_names, partial="no") {
  if (partial=="yes") {
    cap_text = "mean pnacf"
  } else {
    cap_text = "mean nacf"
  }
  frame_caption = get_frame_names_caption(frame_names)
  circle_points = data.frame(rbind(get_circle_points_polar(1, length(frame_names)), c(0,0)))
  colnames(circle_points) <- c("x", "y")
  label_points = get_circle_points_polar(1 - 0.15, length(frame_names))
  p <- ggplot(circle_points, aes(x = .data$x, y = .data$y))
  p + geom_point() +
    coord_fixed() + 
    theme_void() +
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          plot.caption = element_text(hjust = 0.5, size = 14, face = "italic")) +
    annotate("text",  x=label_points[, 1], y=label_points[, 2], label=as.character(seq(1:length(frame_names))), col="black") +
    # annotate("text", x=0.0, y=0.0, label=".", color="black", size=12) + 
    labs(caption = paste0("\n\n\n", "Each circle is associated to a", "\n", "covariate and/or time slice","\n" ,"the circle at the center is the", "\n", cap_text, 
                          " value", "\n\n\n", frame_caption))
}


get_frame_names_caption <- function(frame_names) {
  out = ""
  for (i in 1:length(frame_names)) {
    out = paste0(out, as.character(i), ": ", frame_names[i], "\n")
  }
  return(out)
}
