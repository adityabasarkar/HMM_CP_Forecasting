library(dplyr)
library(grid)

data_dir <- "outputs/r_objects"
plot_dir <- "outputs/plots"

in_path <- file.path(data_dir, "mc_avg_set_sizes.rds")
avg_set_sizes <- readRDS(in_path)

draw_combined_table <- function(data) {
  horizons <- sort(unique(data$test_length))
  coverages <- sort(unique(data$target_coverage))
  
  n_rows <- length(horizons) + 1
  n_cols <- length(coverages) + 1
  
  value_color <- "black"
  header_fill <- "#D9E2F3"
  stripe_fill <- "#F7F7F7"
  border_color <- "#BDBDBD"
  
  grid.newpage()
  
  pushViewport(viewport(
    x = 0.5, y = 0.5,
    width = 0.96, height = 0.94,
    just = c("center", "center")
  ))
  
  grid.text(
    "Average Number of Sequences per Prediction Set",
    x = 0.5, y = 0.98,
    gp = gpar(fontsize = 20, fontface = "bold")
  )

  table_left <- 0.10
  table_right <- 0.98
  table_top <- 0.84
  table_bottom <- 0.14
  table_height <- table_top - table_bottom
  table_width <- table_right - table_left
  row_h <- table_height / n_rows
  col_w <- table_width / n_cols

  grid.text(
    "Target Coverage",
    x = (table_left + table_right) / 2,
    y = 0.90,
    gp = gpar(fontsize = 12, fontface = "bold")
  )

  grid.text(
    "Prediction Horizon",
    x = 0.03,
    y = (table_top + table_bottom) / 2,
    rot = 90,
    gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  for (row_idx in seq_len(n_rows)) {
    y_top <- table_top - (row_idx - 1) * row_h
    y_center <- y_top - row_h / 2
    
    for (col_idx in seq_len(n_cols)) {
      x_left <- table_left + (col_idx - 1) * col_w
      x_center <- x_left + col_w / 2
      
      fill <- "white"
      if (row_idx == 1 || col_idx == 1) {
        fill <- header_fill
      } else if ((row_idx - 1) %% 2 == 0) {
        fill <- stripe_fill
      }
      
      grid.rect(
        x = x_center,
        y = y_center,
        width = col_w,
        height = row_h,
        just = c("center", "center"),
        gp = gpar(fill = fill, col = border_color, lwd = 0.8)
      )
      
      if (row_idx == 1 && col_idx == 1) {
        grid.text(
          "Prediction\nHorizon",
          x = x_center,
          y = y_center,
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      } else if (row_idx == 1) {
        grid.text(
          sprintf("%.2f", coverages[col_idx - 1]),
          x = x_center,
          y = y_center,
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      } else if (col_idx == 1) {
        grid.text(
          as.character(horizons[row_idx - 1]),
          x = x_center,
          y = y_center,
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      } else {
        cell_data <- data %>%
          filter(
            test_length == horizons[row_idx - 1],
            target_coverage == coverages[col_idx - 1]
          )
        
        cp_label <- sprintf("%.2f", cell_data$avg_set_size_cp[[1]])
        lik_label <- sprintf("%.2f", cell_data$avg_set_size_lik[[1]])
        
        grid.text(
          cp_label,
          x = x_center,
          y = y_center + row_h * 0.16,
          gp = gpar(col = value_color, fontsize = 9.5, fontface = "bold")
        )
        grid.text(
          lik_label,
          x = x_center,
          y = y_center - row_h * 0.16,
          gp = gpar(col = value_color, fontsize = 9.5, fontface = "bold")
        )
      }
    }
  }

  legend_y <- 0.06
  grid.text(
    "Top number: CP set",
    x = 0.40,
    y = legend_y,
    gp = gpar(col = value_color, fontsize = 11, fontface = "bold")
  )
  grid.text(
    "|",
    x = 0.50,
    y = legend_y,
    gp = gpar(col = "gray40", fontsize = 11, fontface = "bold")
  )
  grid.text(
    "Bottom number: Likelihood set",
    x = 0.66,
    y = legend_y,
    gp = gpar(col = value_color, fontsize = 11, fontface = "bold")
  )
  
  popViewport()
}

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

pdf(file.path(plot_dir, "avg_set_size_mc.pdf"), width = 14, height = 8)
draw_combined_table(avg_set_sizes)
dev.off()

print("Saved: outputs/plots/avg_set_size_mc.pdf")
