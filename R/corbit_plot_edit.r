corbit_plot <- function(vts, net, max_lag, max_stage, weight_matrix, viridis_color_option="viridis", size_option="absolute_val", partial="no", r_corbit="no", line_segment="no", rectangular_plot="no") {
  if (rectangular_plot == "yes") {
    missing_data_warning(missing_data = (sum(is.na(vts)) != 0))
    corbit_plot_rectangular_two(vts=vts, net=net, max_lag=max_lag, max_stage=max_stage, weight_matrix=weight_matrix, viridis_color_option=viridis_color_option, size_option="absolute_val", partial=partial, r_corbit="no")
  } else if (rectangular_plot == "square") {
    missing_data_warning(missing_data = (sum(is.na(vts)) != 0))
    corbit_plot_rectangular(vts=vts, net=net, max_lag=max_lag, max_stage=max_stage, weight_matrix=weight_matrix, viridis_color_option=viridis_color_option, size_option="absolute_val", partial=partial, r_corbit="no")
  } else {
    dummy_net = GNARtoigraph(net)
    adj_mat = as.matrix(dummy_net)
    stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
    if (missing(weight_matrix)) {
      W_norm = weights_matrix(net)
    } else {
      W_norm = weight_matrix
    }
    missing_data_warning(missing_data = (sum(is.na(vts)) != 0))
    corbit_data = get_corbit_data(max_lag, max_stage, W_norm, stages_tensor, vts, size_option, partial = partial)
    corbit_values = matrix(corbit_data[, 3], nrow = max_lag, ncol = max_stage)
    colnames(corbit_values) <- vapply(c(1:max_stage), function(x) {paste0("r-stage: ", as.character(x))}, "")
    if (line_segment == "no") {
      plot_corbit(corbit_data = corbit_data, max_lag = max_lag, max_stage = max_stage, viridis_color_option = viridis_color_option, size_option = size_option, 
      r_corbit=r_corbit)
      invisible(corbit_values)
    } else {
      plot_corbit_line_segment(corbit_data = corbit_data, max_lag = max_lag, max_stage = max_stage, viridis_color_option = viridis_color_option, size_option = size_option, 
      r_corbit=r_corbit)
      invisible(corbit_values)
    }
  }
}


r_corbit_plot <- function(vts_frames, network_list, max_lag, max_stage, weight_matrices, frame_names, same_net="no", viridis_color_option="viridis", size_option="absolute_val", partial="no", r_corbit="yes") {
  if (missing(frame_names)) {
    cat("No name for each covariate and/or time slice, using generic names")
    frame_names = as.character(seq(1:length(vts_frames)))
  }
  if (length(frame_names) != length(vts_frames)) {
    cat("Missing name for slice and/or covariate, using generic names")
    frame_names = as.character(seq(1:length(vts_frames)))
  }
  r_corbit_stages_weights = get_weight_matrices(network_list = network_list, max_stage = max_stage, weight_matrices = weight_matrices, same_net=same_net)
  stages_tensor_list = r_corbit_stages_weights[[1]]
  weight_matrices_list = r_corbit_stages_weights[[2]]
  if (same_net == "yes") {
    stages_tensor_list = list()
    weight_matrices_list = list()
    for (i in 1:length(vts_frames)) {
      stages_tensor_list = c(stages_tensor_list, list(r_corbit_stages_weights[[1]]))
      weight_matrices_list = c(weight_matrices_list, list(r_corbit_stages_weights[[2]]))
    }
  }
  number_levels = length(frame_names)
  r_corbit_plot_data = get_r_corbit_data(max_lag, max_stage, weight_matrices_list, stages_tensor_list, vts_frames, same_net=same_net, partial = partial)
  r_corbit_values = lapply(c(1:number_levels), function(x) {matrix(r_corbit_plot_data[seq(x + 1, nrow(r_corbit_plot_data) + x - number_levels, number_levels + 1), 3], nrow = max_lag, ncol = max_stage)})
  names(r_corbit_values) <- frame_names
  for (j in 1:number_levels) {
    colnames(r_corbit_values[[j]]) <- vapply(c(1:max_stage), function(x) {paste0("r-stage: ", as.character(x))}, "")
  }
  print(ggarrange(plot_corbit(r_corbit_plot_data, max_lag, max_stage,  viridis_color_option, size_option, r_corbit=r_corbit), 
            plot_r_corbit_legend(frame_names, partial=partial), widths = c(2, 0.5),
            ncol = 2, nrow = 1)
  )
  invisible(r_corbit_values)
}
