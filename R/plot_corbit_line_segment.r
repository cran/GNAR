plot_corbit_line_segment <-function(corbit_data, max_lag, max_stage, viridis_color_option="viridis", size_option="absolute_val", r_corbit="no") {
  corbit_data = rbind(c(0, 0, 0, 0), corbit_data)
  aux_line_start = sin((max_lag - 1) * 2 * pi / max_lag)
  if (r_corbit == "yes") {
    label_points = get_circle_points_polar(2 * max_stage + 1, max_lag)
    num_samples = (length(corbit_data[, 1]) - 1) / (max_lag * max_stage) - 1
    capt_text = paste0("R-Corbit plot with ", num_samples, " time frames or samples, max lag")
  } else {
    label_points = get_circle_points_polar(max_stage + 1, max_lag)
    capt_text = "Corbit plot with max lag"
  }
  aux_size = 4
  if (size_option == "absolute_val") {
    aux_size = 3
    colnames(corbit_data) <- c("x", "y", colnames(corbit_data)[3], paste0("abs(", colnames(corbit_data)[3], ")"))
  }
  plot_object <- ggplot(corbit_data, aes(x = .data$x, y = .data$y))
  print(plot_object + geom_point(aes(color = .data[[colnames(corbit_data)[3]]], size = abs(.data[[colnames(corbit_data)[3]]]))) +
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
      geom_segment(aes(x = corbit_data[2, 1], y = aux_line_start, xend = label_points[1, 1], yend = (label_points[1, 2] + label_points[max_lag, 2]) / 2), linetype = "dashed") +
      annotate("text", x=label_points[, 1], y=label_points[, 2], label=rownames(data.frame(label_points)), col="black", size = 5) +
      labs(caption = paste(colnames(corbit_data)[3], capt_text, max_lag, "and max path length", max_stage, sep = " " ), 
             color = colnames(corbit_data)[3],
             size = colnames(corbit_data)[4])
  )
}