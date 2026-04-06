library(ggplot2)
library(dplyr)
library(scales)

# Generalized Manhattan Plot Function
plot_manhattan <- function(data, 
                          chr_col = "chr",           # Chromosome column name
                          pos_col = "midPos",        # Position column name
                          value_col = "PBS0",        # Value column to plot
                          transform = "none",        # "none", "log10", or "neglog10"
                          title = "Manhattan Plot",
                          significance_threshold = NULL,  # e.g., 0.05 for top 5%
                          threshold_value = NULL,    # Absolute threshold line (e.g., -log10(5e-8))
                          highlight_top_n = NULL,    # Highlight top N points
                          point_size = 2,
                          colors = c("#E69F00", "#56B4E9"),  # Alternating chr colors
                          highlight_color = "red",
                          ylim = NULL,
                          xlabel = "Chromosome",
                          ylabel = NULL) {           # Auto-generated if NULL
  
  # Rename columns for easier processing
  plot_data <- data %>%
    rename(chr = !!sym(chr_col),
           pos = !!sym(pos_col),
           value = !!sym(value_col))
  
  # Apply transformation
  if (transform == "log10") {
    plot_data <- plot_data %>%
      mutate(plot_value = log10(value))
    default_ylabel <- paste0("log10(", value_col, ")")
  } else if (transform == "neglog10") {
    plot_data <- plot_data %>%
      mutate(plot_value = -log10(value))
    default_ylabel <- paste0("-log10(", value_col, ")")
  } else {
    plot_data <- plot_data %>%
      mutate(plot_value = value)
    default_ylabel <- value_col
  }
  
  # Use provided ylabel or default
  if (is.null(ylabel)) {
    ylabel <- default_ylabel
  }
  
  # Remove NA values
  plot_data <- plot_data %>%
    filter(!is.na(plot_value), !is.infinite(plot_value))
  
  # Extract chromosome number for proper numeric sorting
  plot_data <- plot_data %>%
    mutate(chr_num = as.numeric(gsub("\\D", "", chr))) %>%
    # Handle any NA values from non-numeric chromosomes (X, Y, MT)
    mutate(chr_num = case_when(
      grepl("X", chr, ignore.case = TRUE) ~ max(chr_num, na.rm = TRUE) + 1,
      grepl("Y", chr, ignore.case = TRUE) ~ max(chr_num, na.rm = TRUE) + 2,
      grepl("M", chr, ignore.case = TRUE) ~ max(chr_num, na.rm = TRUE) + 3,
      TRUE ~ chr_num
    )) %>%
    arrange(chr_num, pos)
  
  # Calculate cumulative positions for x-axis
  chr_lengths <- plot_data %>%
    group_by(chr, chr_num) %>%
    summarise(chr_len = max(pos), .groups = "drop") %>%
    arrange(chr_num) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len)
  
  plot_data <- plot_data %>%
    left_join(chr_lengths %>% select(chr, tot), by = "chr") %>%
    mutate(BPcum = pos + tot)
  
  # Calculate axis labels position (maintaining proper order)
  axis_df <- plot_data %>%
    group_by(chr, chr_num) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop") %>%
    arrange(chr_num)
  
  # Determine which points to highlight
  if (!is.null(highlight_top_n)) {
    plot_data <- plot_data %>%
      mutate(is_highlight = rank(-plot_value) <= highlight_top_n)
  } else if (!is.null(significance_threshold)) {
    threshold_val <- quantile(plot_data$plot_value, 
                              probs = 1 - significance_threshold, 
                              na.rm = TRUE)
    plot_data <- plot_data %>%
      mutate(is_highlight = plot_value >= threshold_val)
  } else {
    plot_data <- plot_data %>%
      mutate(is_highlight = FALSE)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = BPcum, y = plot_value)) +
    # Alternating chromosome colors
    geom_point(aes(color = as.factor(chr_num %% 2)), 
               size = point_size, alpha = 0.7) +
    scale_color_manual(values = colors) +
    
    # Highlighted points
    geom_point(data = subset(plot_data, is_highlight), 
               color = highlight_color, 
               size = point_size * 1.5, 
               alpha = 1) +
    
    # Custom X axis
    scale_x_continuous(label = axis_df$chr, 
                      breaks = axis_df$center,
                      expand = c(0.01, 0.01)) +
    
    # Y axis
    scale_y_continuous(expand = c(0.02, 0.02)) +
    
    # Add threshold line if specified (for significance threshold)
    {if (!is.null(significance_threshold)) {
      threshold_val <- quantile(plot_data$plot_value, 
                               probs = 1 - significance_threshold, 
                               na.rm = TRUE)
      geom_hline(yintercept = threshold_val, 
                linetype = "dashed", 
                color = "gray40", 
                linewidth = 0.5)
    }} +
    
    # Add absolute threshold line if specified (for GWAS significance)
    {if (!is.null(threshold_value)) {
      geom_hline(yintercept = threshold_value, 
                linetype = "dashed", 
                color = "red", 
                linewidth = 0.7)
    }} +
    
    # Themes and labels
    labs(x = xlabel, 
         y = ylabel,
         title = title) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 0, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  # Apply y-axis limits if specified
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  return(p)
}
