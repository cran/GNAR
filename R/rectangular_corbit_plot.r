# rectangular corbit plot

plot_rectangular_corbit <- function(corbit_data, max_lag, max_stage, viridis_color_option="viridis", size_option="absolute_val", r_corbit="no") {
  capt_text = "Corbit plot with max lag"
  #corbit_data = rbind(c(0, -1/3, -1, 0.1), corbit_data)
  aux_size = 4
  if (size_option == "absolute_val") {
    aux_size = 3
    colnames(corbit_data) <- c("x", "y", colnames(corbit_data)[3], paste0("abs(", colnames(corbit_data)[3], ")"))
  }
  x_rectangular = c(vapply(c(1:max_stage), function(x) {return (c(1:max_lag))}, rep(0, max_lag)))
  y_rectangular = c(vapply(c(1:max_stage), function(x) {return (rep(x, max_lag))}, rep(0, max_lag)))
  corbit_data[, 1] = x_rectangular
  corbit_data[, 2] = y_rectangular
  p <- ggplot(corbit_data, aes(x = .data$x, y= .data$y))
  print(p + geom_point(aes(color = .data[[colnames(corbit_data)[3]]], size = abs(.data[[colnames(corbit_data)[3]]]))) +
        scale_color_viridis(option=viridis_color_option) +
      theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 0, vjust = 0.5)
        ) +
      labs(title = paste(colnames(corbit_data)[3], capt_text, max_lag, "and max path length", max_stage, sep = " " ), 
             color = colnames(corbit_data)[3],
             size = colnames(corbit_data)[4]) +
      geom_vline(xintercept = c(1:max_lag), color = 'gray', linetype = 'dotted') + 
      scale_x_continuous("Lag", labels = as.character(c(1:max_lag)), breaks = c(1:max_lag)) +
      scale_y_continuous("Stage", labels = as.character(c(1:max_stage)), breaks = c(1:max_stage))
  )
}


corbit_plot_rectangular <- function(vts, net, max_lag, max_stage, weight_matrix, viridis_color_option="viridis", size_option="absolute_val", partial="no", r_corbit="no") {
  dummy_net = GNARtoigraph(net)
  adj_mat = as.matrix(dummy_net)
  # net_result = graph_from_adjacency_matrix(adj_mat, 'undirected')
  stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
  if (missing(weight_matrix)) {
    W_norm = weights_matrix(net)
  } else {
    W_norm = weight_matrix
  }
  missing_data_warning(missing_data = length(is.na(vts)) == 0)
  corbit_data = get_corbit_data(max_lag, max_stage, W_norm, stages_tensor, vts, size_option, partial = partial)
  corbit_values = matrix(corbit_data[, 3], nrow = max_lag, ncol = max_stage)
  colnames(corbit_values) <- vapply(c(1:max_stage), function(x) {paste0("r-stage: ", as.character(x))}, "")
  plot_rectangular_corbit(corbit_data, max_lag, max_stage,  viridis_color_option, size_option, r_corbit=r_corbit)
  invisible(corbit_values)
}


corbit_plot_rectangular_two <- function(vts, net, max_lag, max_stage, weight_matrix, viridis_color_option="viridis", size_option="absolute_val", partial="no", r_corbit="no") {
  dummy_net = GNARtoigraph(net)
  adj_mat = as.matrix(dummy_net)
  # net_result = graph_from_adjacency_matrix(adj_mat, 'undirected')
  stages_tensor = get_k_stages_adjacency_tensor(adj_mat, max_stage)
  if (missing(weight_matrix)) {
    W_norm = weights_matrix(net)
  } else {
    W_norm = weight_matrix
  }
  missing_data_warning(missing_data = length(is.na(vts)) == 0)
  if (partial == "yes") {
    corbit_data = get_pnacf_values(max_lag, max_stage, W_norm, stages_tensor, vts)
    plot_y_axis_name = "PNACF"
  } else {
    corbit_data = get_nacf_values(max_lag, max_stage, W_norm, stages_tensor, vts)
    plot_y_axis_name = "NACF"
  }
  r_stage_numbers = c(1:max_stage)
  colour_vector = palette(hcl.colors(max_stage, "viridis"))
  # more accurate pnacf values on y-axis ylim = c(min(corbit_data) - 0.25, max(corbit_data) + 0.25)
  plot(c(1:max_lag), corbit_data[, 1], type = "b", col=colour_vector[1], pch=as.character(r_stage_numbers[1]), ylim = c(min(corbit_data) - 0.25, max(corbit_data) + 0.25), xlab="Lag", ylab=plot_y_axis_name, main=paste0(plot_y_axis_name, " plot with max lag: ", as.character(max_lag),
        " and max r-stage: ", as.character(max_stage)))
  if (max_stage > 1) {
    for (i in 2:max_stage) {
      lines(c(1:max_lag), corbit_data[, i], type = "b", col=colour_vector[i], pch=as.character(r_stage_numbers[i]))
    }
  }
  legend("topright", legend=vapply(c(1:max_stage), function(x) {paste0("r-stage: ", as.character(x))}, ""), col=colour_vector, lty=rep(1, max_stage), cex=0.8)
  colnames(corbit_data) <- vapply(c(1:max_stage), function(x) {paste0("r-stage: ", as.character(x))}, "")
  invisible(corbit_data)
}