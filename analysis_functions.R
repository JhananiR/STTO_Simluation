library(ggrepel)
library(scales)
library(RColorBrewer)

extract_rda_scores <- function(rda_model, 
                              population_data = NULL,
                              pop_column = NULL,
                              color_column = NULL,
                              scaling = 2) {
  
  # Extract site scores (populations)
  site_scores <- as.data.frame(scores(rda_model, display = "sites", scaling = scaling))
  site_scores$Population <- rownames(site_scores)
  
  # Extract species scores (alleles) - optional for large datasets
  species_scores <- as.data.frame(scores(rda_model, display = "species", scaling = scaling))
  species_scores$Allele <- rownames(species_scores)
  
  # Extract environmental/spatial variable scores (biplot arrows)
  if (length(rda_model$CCA$eig) > 0) {  # Check if constrained axes exist
    env_scores <- as.data.frame(scores(rda_model, display = "bp", scaling = scaling))
    env_scores$Variable <- rownames(env_scores)
  } else {
    env_scores <- NULL
  }
  
  # Add population metadata if provided
  if (!is.null(population_data) && !is.null(pop_column)) {
    # Match population data
    pop_data_subset <- population_data[, c(pop_column, color_column)]
    names(pop_data_subset)[1] <- "Population"
    site_scores <- merge(site_scores, pop_data_subset, by = "Population", all.x = TRUE)
  }
  
  return(list(
    sites = site_scores,
    species = species_scores,
    env_vectors = env_scores,
    eigenvalues = rda_model$CCA$eig,
    variance_explained = summary(rda_model)
  ))
}

plot_rda_biplot <- function(rda_scores, 
                           axes = c(1, 2),
                           color_by = NULL,
                           show_species = FALSE,
                           species_alpha = 0.3,
                           vector_scale = 2,
                           point_size = 3,
                           title = "RDA Biplot") {
  
  # Prepare axis names
  axis_names <- paste0("RDA", axes)
  x_var <- axis_names[1]
  y_var <- axis_names[2]
  
  # Calculate variance explained for each axis
  if (!is.null(rda_scores$eigenvalues)) {
    total_var <- sum(rda_scores$eigenvalues)
    x_var_exp <- round((rda_scores$eigenvalues[axes[1]] / total_var) * 100, 1)
    y_var_exp <- round((rda_scores$eigenvalues[axes[2]] / total_var) * 100, 1)
    
    x_label <- paste0(x_var, " (", x_var_exp, "%)")
    y_label <- paste0(y_var, " (", y_var_exp, "%)")
  } else {
    x_label <- x_var
    y_label <- y_var
  }
  
  # Start building the plot
  p <- ggplot()
  
  # Add species scores (alleles) if requested and data exists
  if (show_species && !is.null(rda_scores$species) && nrow(rda_scores$species) < 1000) {
    p <- p + 
      geom_point(data = rda_scores$species, 
                aes_string(x = x_var, y = y_var), 
                color = "grey70", alpha = species_alpha, size = 0.5)
  }
  
  # Add site scores (populations)
  if (!is.null(color_by) && color_by %in% names(rda_scores$sites)) {
    p <- p + 
      geom_point(data = rda_scores$sites, 
                aes_string(x = x_var, y = y_var, color = color_by), 
                size = point_size, alpha = 0.8)
  } else {
    p <- p + 
      geom_point(data = rda_scores$sites, 
                aes_string(x = x_var, y = y_var), 
                color = "steelblue", size = point_size, alpha = 0.8)
  }
  
  # Add population labels
  p <- p + 
    geom_text_repel(data = rda_scores$sites,
                   aes_string(x = x_var, y = y_var, label = "Population"),
                   size = 6, max.overlaps = 10)
  
  # Add environmental vectors if they exist
  if (!is.null(rda_scores$env_vectors)) {
    # Scale vectors for better visualization
    env_data <- rda_scores$env_vectors
    env_data[, x_var] <- env_data[, x_var] * vector_scale
    env_data[, y_var] <- env_data[, y_var] * vector_scale
    
    p <- p + 
      geom_segment(data = env_data,
                  aes_string(x = 0, y = 0, xend = x_var, yend = y_var),
                  arrow = arrow(length = unit(0.3, "cm")), 
                  color = "red", size = 1, alpha = 0.8) +
      geom_text_repel(data = env_data,
                     aes_string(x = x_var, y = y_var, label = "Variable"),
                     color = "red", fontface = "plain", size = 6)
  }
  
  # Customize theme
  p <- p + 
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 14)
    ) +
    labs(
      x = x_label,
      y = y_label, 
      title = title
    ) +
    coord_fixed()
  
  return(p)
}
# Function to create variance partitioning pie chart
plot_variance_pie <- function(varpart_result, title = "Variance Partitioning: IBD vs IBE") {
  
  # Extract variance components (adjust these based on your varpart object structure)
  # Typically: [1] = pure environmental, [2] = pure spatial, [3] = shared, [4] = residual
  
  # Example data structure - replace with your actual values
  variance_data <- data.frame(
    Component = c("Pure Spatial\n(IBD)", "Pure Environmental\n(IBE)", 
                  "Shared\n(Spatial-Environmental)", "Unexplained"),
    Variance = c(varpart_result$part$indfract[2,3],  # Pure spatial
                 varpart_result$part$indfract[1,3],  # Pure environmental  
                 varpart_result$part$indfract[3,3],  # Shared
                 1 - sum(varpart_result$part$indfract[1:3,3])), # Residual
    stringsAsFactors = FALSE
  )
  
  # Remove negative values (can happen in variance partitioning)
  variance_data$Variance[variance_data$Variance < 0] <- 0
  
  # Calculate percentages
  variance_data$Percentage <- round(variance_data$Variance * 100, 1)
  
  # Create labels
  variance_data$Label <- paste0(variance_data$Component, "\n", 
                               variance_data$Percentage, "%")
  
  # Define colors
  colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#CCCCCC")  # Red, Blue, Green, Gray
  
  # Create pie chart
  p <- ggplot(variance_data, aes(x = "", y = Variance, fill = Component)) +
    geom_col(width = 1, color = "white", size = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = title) +
    geom_text(aes(label = paste0(Percentage, "%")), 
              position = position_stack(vjust = 0.5),
              size = 4, fontface = "bold", color = "white")
  
  return(p)
}

# Function to create variance partitioning bar chart
plot_variance_bars <- function(varpart_result, title = "Variance Components") {
  
  # Prepare data
  variance_data <- data.frame(
    Component = c("Pure Spatial\n(IBD)", "Pure Environmental\n(IBE)", 
                  "Shared", "Unexplained"),
    Variance = c(varpart_result$part$indfract[1,3],  # Pure spatial
                 varpart_result$part$indfract[2,3],  # Pure environmental
                 varpart_result$part$indfract[3,3],  # Shared
                 1 - sum(varpart_result$part$indfract[1:3,3])), # Residual
    Type = c("Spatial", "Environmental", "Shared", "Unexplained")
  )
  
  # Remove negative values
  variance_data$Variance[variance_data$Variance < 0] <- 0
  
  # Calculate percentages
  variance_data$Percentage <- round(variance_data$Variance * 100, 1)
  
  # Reorder by variance explained
  variance_data$Component <- factor(variance_data$Component, 
                                   levels = variance_data$Component[order(variance_data$Variance, decreasing = TRUE)])
  
  # Define colors
  colors <- c("Spatial" = "#E31A1C", "Environmental" = "#1F78B4", 
              "Shared" = "#33A02C", "Unexplained" = "#CCCCCC")
  
  p <- ggplot(variance_data, aes(x = Component, y = Percentage, fill = Type)) +
    geom_col(color = "white", size = 1, alpha = 0.8) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, max(variance_data$Percentage) + 5),
                      labels = function(x) paste0(x, "%")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      y = "Variance Explained (%)",
      title = title
    ) +
    geom_text(aes(label = paste0(Percentage, "%")), 
              vjust = -0.5, size = 4, fontface = "bold")
  
  return(p)
}

# Function to create Venn diagram-style overlap visualization
plot_variance_venn_style <- function(varpart_result, title = "IBD vs IBE Overlap") {
  
  # Extract values
  pure_spatial <- max(0, varpart_result$part$indfract[1,3])
  pure_env <- max(0, varpart_result$part$indfract[2,3])
  shared <- max(0, varpart_result$part$indfract[3,3])
  
  # Calculate percentages
  total_explained <- pure_spatial + pure_env + shared
  spatial_pct <- round((pure_spatial + shared) / total_explained * 100, 1)
  env_pct <- round((pure_env + shared) / total_explained * 100, 1)
  shared_pct <- round(shared / total_explained * 100, 1)
  
  # Create data for visualization
  overlap_data <- data.frame(
    x = c(-1, 1, 0),
    y = c(0, 0, 0),
    size = c(pure_spatial + shared, pure_env + shared, shared),
    label = c(paste0("Spatial\n", round(pure_spatial*100, 1), "%"),
              paste0("Environmental\n", round(pure_env*100, 1), "%"),
              paste0("Shared\n", round(shared*100, 1), "%")),
    color = c("Spatial", "Environmental", "Shared")
  )
  
  p <- ggplot(overlap_data, aes(x = x, y = y)) +
    geom_point(aes(size = size, color = color), alpha = 0.6) +
    scale_size_continuous(range = c(10, 30), guide = "none") +
    scale_color_manual(values = c("Spatial" = "#E31A1C", "Environmental" = "#1F78B4", "Shared" = "#33A02C")) +
    geom_text(aes(label = label), size = 4, fontface = "bold") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(title = title) +
    xlim(-2, 2) + ylim(-1, 1)
  
  return(p)
}

# Function to create stacked area chart showing the breakdown
plot_variance_stacked <- function(varpart_result, title = "Variance Decomposition") {
  
  # Prepare data
  variance_data <- data.frame(
    x = 1,
    Pure_Spatial = max(0, varpart_result$part$indfract[1,3]),
    Pure_Environmental = max(0, varpart_result$part$indfract[2,3]),
    Shared = max(0, varpart_result$part$indfract[3,3]),
    Unexplained = max(0, 1 - sum(varpart_result$part$indfract[1:3,3]))
  )
  
  # Reshape for plotting
  variance_long <- variance_data %>%
    tidyr::pivot_longer(cols = -x, names_to = "Component", values_to = "Variance") %>%
    mutate(
      Component = factor(Component, levels = c("Unexplained", "Shared", "Pure_Environmental", "Pure_Spatial")),
      Percentage = round(Variance * 100, 1),
      y_pos = cumsum(Variance) - Variance/2
    )
  
  # Define colors
  colors <- c("Pure_Spatial" = "#E31A1C", "Pure_Environmental" = "#1F78B4", 
              "Shared" = "#33A02C", "Unexplained" = "#CCCCCC")
  
  p <- ggplot(variance_long, aes(x = x, y = Variance, fill = Component)) +
    geom_col(width = 0.5, color = "white", size = 1) +
    scale_fill_manual(values = colors,
                     labels = c("Unexplained", "Shared\n(Spatial-Environmental)", 
                               "Pure Environmental\n(IBE)", "Pure Spatial\n(IBD)")) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "",
      y = "Proportion of Variance Explained",
      title = title
    ) +
    geom_text(aes(y = y_pos, label = paste0(Percentage, "%")), 
              size = 4, fontface = "bold", color = "white")
  
  return(p)
}