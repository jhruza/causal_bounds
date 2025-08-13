library(ggplot2)
library(gridExtra)
library(grid)

results <- readRDS("results.rds")
reduced <- (results$true_estimand_values - results$bound_reduced_lower) / (results$bound_reduced_upper - results$bound_reduced_lower) # nolint: line_length_linter.
conditional <- (results$true_estimand_values - results$bound_conditional_lower) / (results$bound_conditional_upper - results$bound_conditional_lower) # nolint: line_length_linter.


# Create the first plot (your current plot)
p1 <- ggplot(data.frame(value = reduced), aes(x = value)) +
  geom_density(fill = "lightblue", color = NA,  alpha = 0.4) +
  geom_rug(sides = "b", length = unit(0.2, "npc"), alpha = 0.7) +
  geom_vline(xintercept = c(0, 1),  color = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    # aspect.ratio = 0.8
  ) +
  labs(x = "")

# Create the second plot with transformed data 
# (I'm using a different transformation as an example - adjust as needed)
p2 <- ggplot(data.frame(value = conditional), aes(x = value)) +
  geom_density(fill = "lightgreen", color = NA,  alpha = 0.4) +
  geom_rug(sides = "b", length = unit(0.2, "npc"), alpha = 0.7) +
  geom_vline(xintercept = c(0, 1), color = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    # aspect.ratio = 0.8
  ) +
  labs(x = "")

# Create empty plots with just titles
title1_plot <- ggplot() + 
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Reduced bounds", size = 6) +
  theme(plot.background = element_rect(fill = "white", color = NA))

title2_plot <- ggplot() + 
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = "Conditional bounds", size = 6) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Arrange the plots in a 2x2 grid
combined_plot <- grid.arrange(
  title1_plot, p1,
  title2_plot, p2,
  ncol = 2,
  widths = c(1, 3)  # Left column narrower, right column wider
)

# Print and save the combined plot
print(combined_plot)
ggsave("grid_plots.png", combined_plot, width = 10, height = 3, units = "in", dpi = 300)
