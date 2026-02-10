# plot_R2_elbow_analysis.R

library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

load("analysis_R0.RData")   # objeto: analysis_R0$results
load("analysis_R30.RData")  # objeto: analysis_R30$results

# ----------------------------------------------------------
# [1] Função para preparar os dados com detecção de "elbow"
# ----------------------------------------------------------
prepare_r2_data <- function(r2_results, condition_label) {
  r2_values <- r2_results$r_squared
  
  # Segunda derivada para detectar "cotovelo"
  second_deriv <- diff(diff(r2_values))
  elbow_x <- which.max(abs(second_deriv)) + 2
  elbow_y <- r2_values[elbow_x]
  
  return(list(
    plot_data = data.frame(
      n_PCs = r2_results$n_pcs,
      R2 = r2_values,
      Condition = condition_label
    ),
    elbow_point = data.frame(
      x = elbow_x,
      y = elbow_y,
      label = sprintf("Elbow (n = %d)", elbow_x)
    )
  ))
}

# [2] Preparar dados
data_R0 <- prepare_r2_data(analysis_R0$results, "R0")
data_R30 <- prepare_r2_data(analysis_R30$results, "R30")

plot_data <- bind_rows(data_R0$plot_data, data_R30$plot_data)
elbow_points <- bind_rows(
  data_R0$elbow_point %>% mutate(Condition = "R0"),
  data_R30$elbow_point %>% mutate(Condition = "R30")
)

# ----------------------------------------------------------
# [3] Função de plotagem sem barras, com escala log
# ----------------------------------------------------------
create_r2_plot <- function(data, elbows, max_pcs, title_suffix) {
  filtered_data <- data %>% filter(n_PCs <= max_pcs)
  filtered_elbows <- elbows %>% filter(x <= max_pcs)
  
  ggplot(filtered_data, aes(x = n_PCs, y = R2, color = Condition)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_vline(data = filtered_elbows, aes(xintercept = x, color = Condition),
               linetype = "dashed", linewidth = 1) +
    geom_point(data = filtered_elbows, aes(x = x, y = y), 
               shape = 18, size = 6, inherit.aes = FALSE) +
    geom_text_repel(data = filtered_elbows,
                    aes(x = x, y = y, label = label, color = Condition),
                    size = 5, nudge_x = max_pcs/15) +
    labs(
      title = paste("Reconstruction R² with", title_suffix),
      x = "Number of Principal Components",
      y = expression(R^2 ~ " of reconstruction (log10 scale)"),
      color = "Condition"
    ) +
    scale_y_log10(
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_color_manual(values = c("R0" = "#1f78b4", "R30" = "#e6550d")) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# ----------------------------------------------------------
# [4] Criar e exportar os gráficos
# ----------------------------------------------------------
if (!dir.exists("images")) dir.create("images")

# Apenas 50 PCs e todos
p_50  <- create_r2_plot(plot_data, elbow_points, max_pcs = 50, title_suffix = "Up to 50 PCs")
p_all <- create_r2_plot(plot_data, elbow_points, max_pcs = max(plot_data$n_PCs), title_suffix = "All PCs")

# Combinar os dois gráficos verticalmente
combined_plot <- (p_50 / p_all) +
  plot_annotation(
    title = "Reconstruction R² by Number of Principal Components (Log Scale)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  ) &
  theme(legend.position = "bottom")

# Exportar em alta qualidade
ggsave(
  filename = "images/R2_reconstruction_elbowplot_combined.pdf",
  plot = combined_plot,
  device = cairo_pdf,
  width = 10, height = 12,
  dpi = 600,
  units = "in"
)
