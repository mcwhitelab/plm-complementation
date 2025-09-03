library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
print("loaded")

library(ggrastr)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot(font_size = 8, rel_small = 1, rel_large = 1, rel_tiny = 1, font_family = "Arial"))

# Read the data
data <- read_csv("data/compare_paralogs_output.txt")
print(colnames(data))
#print(data$num_embeddings %>% unique)
print(data$method %>% unique)

# Parse embedding filename into components using separate
parsed_data <- data %>%
  filter(model != "progen2-large") %>%
  # Extract the similarity filename from embedding_filename column
  mutate(similarity_filename = embedding_filename) %>%
  # Split the filename into meaningful parts
  separate(
    similarity_filename,
    into = c("input_file", "file_ext", "layer_config_parsed", "model_parsed", 
             "embedding_type", "idx", "similarities_method"),
    sep = "\\.",
    remove = FALSE,
    fill = "right"
  ) %>%
  # Fix the input file name which got split
  mutate(
    input_file = paste(input_file, file_ext, sep = "."),
    # Clean up the similarities_method column by removing "similarities_" prefix
    similarity_method_parsed = gsub("similarities_", "", similarities_method)
  ) %>%
  # Remove the temporary columns
  select(-file_ext, -idx) %>%
  mutate(model = case_when(model == "ESMplusplus_small" ~ "ESMC-300M",
                           model == "ESMplusplus_large" ~ "ESMC-600M",
                           TRUE ~ model))

parsed_data %>% 
write_csv("figures/Fig2_underlyingdata.csv")
# Summarize results by embedding type, model, and layer configuration
results <- parsed_data %>%
  group_by(similarity_method_parsed, embedding_type, model, layer_config) %>%
  summarize(
    total = n(),
    correct = sum(complement_higher == TRUE),
    accuracy = correct / total,
    .groups = "drop"
  ) %>%

  filter(!is.na(similarity_method_parsed))

print("Summary of results:")
print(results)
results %>%
write_csv("Fig2_underlying_data_summarized.csv")

# Using fixed baseline value of 19 for Sequence Alignment
alignment_baseline <- 19
total_cases <- 22

print(paste("Sequence Alignment baseline correct count:", alignment_baseline))
print(paste("Total number of cases:", total_cases))

print(results$model %>% unique)


# Define model ordering (t5 first, then ESM2 by size, then ESM++)
model_order <- c(
  "proST-esm1b",
  "progen2-small",
  "progen2-medium",
  "prot_bert_bfd",
  "prot_t5_xl_uniref50",    
  "esm2_t33_650M_UR50D",
  "esm2_t30_150M_UR50D",
  "esm2_t12_35M_UR50D",
  "esm2_t6_8M_UR50D",          # ESM2 models in order of increasing size
 


  "ESMC-300M",         # ESM++ models at the end
  "ESMC-600M",
  "Sequence Identity"
)
results %>% filter(!model %in% model_order)


# Create the embedding type comparison plot
print(results$model %>% unique)
print(results$layer_config %>% unique)
layer_config_order <- 
  c("Last layer only",
    "Last 4 layers",
    "Last 8 layers")

# Create plot functions to ensure consistent styling
create_complementation_plot <- function(data, facet = TRUE) {
  plot_data <- data %>%
    #filter(model != "prot_bert_bfd")
    mutate(
      embedding_type_label = case_when(
        embedding_type == "sequence_embeddings" ~ "Mean Averaging",
        embedding_type == "sequence_embeddings_swe" ~ "SWE",
        TRUE ~ embedding_type
      ),
      model = factor(model, levels = model_order),
      embedding_type_label = factor(embedding_type_label, levels = c("Mean Averaging", "SWE")),
      highlight = correct > alignment_baseline,
      layer_config = case_when(layer_config == "last_1" ~ "Last layer only",
                               layer_config == "last_4" ~ "Last 4 layers",
                               layer_config == "last_8" ~ "Last 8 layers"),
      layer_config = factor(layer_config, levels = layer_config_order),
    )
  
  
  p <- ggplot(plot_data, aes(x = model, y = correct, fill = embedding_type_label)) +

    
    geom_hline(yintercept = total_cases, linetype = "solid", color = "#0e8e42", linewidth = 1) +
    geom_hline(yintercept = alignment_baseline, linetype = "solid", color = "#0e8e42", linewidth = 0.5) +
    geom_bar(stat = "identity", position = "dodge") +
    #geom_bar(aes(linetype = highlight, color = highlight), stat = "identity", position = "dodge") +
    scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = NULL), guide = "none") +


    annotate("text", x = -Inf, y = total_cases, label = "Total cases (22)", 
             hjust = -0.1, vjust = -0.5, color = "#0e8e42", fontface = "italic", size = 8/2.8) +
    annotate("text", x = -Inf, y = alignment_baseline, label = "Alignment baseline (19)", 
             hjust = -0.1, vjust = -0.5, color = "#0e8e42", fontface = "italic", size = 8/2.8) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = "Number Correct (out of 22)",
      fill = "Embedding Type"
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      strip.text = element_text(size = 8),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank()
    ) +
    scale_fill_manual(values = c("grey80", "grey60")) +
    ylim(0, 25) +
    geom_text(
      aes(label = correct), 
      position = position_dodge(width = 0.9), 
      vjust = -0.5,
      size = 8/2.8
    ) 
  
  if (facet) {
    p <- p + facet_wrap(~layer_config)
  }
  
  return(p)
}


# 1. Create Supplemental Figure 1 (all layer configs)
embperf_layerchoice<- results %>%
  filter(similarity_method_parsed == "cosine") %>%
  filter(layer_config != "all_layers") %>%
  create_complementation_plot(facet = TRUE)

embperf_layerchoice

# Save Supplemental Figure 1
ggsave("figures/SupplementalFigure1.pdf", embperf_layerchoice, width = 12, height = 6)
ggsave("figures/SupplementalFigure1.png", embperf_layerchoice, width = 12, height = 6, bg = "white")
print("Plot saved as SupplementalFigure1.pdf")

# 2. Create Figure 2 (just last_8 layer config)
embperf <- results %>%
  filter(similarity_method_parsed == "cosine") %>%
  filter(layer_config == "last_8") %>%
  create_complementation_plot(facet = FALSE)

embperf

# Save Figure 2
ggsave("data/embperf.pdf", embperf, width = 6, height = 6)
print("Plot saved as embperf.pdf")



# 4. Create a paired dot plot showing complement vs non-complement similarities
# First, extract ESMplusplus_large with last_4 and SWE Embeddings data
esm_plus_plus_data <- parsed_data %>%
  filter(model == "ESMC-600M",
         layer_config == "last_1",
         embedding_type == "sequence_embeddings_swe",
         similarity_method_parsed == "cosine") %>%
  select(yeast_id, complementing_human_id, noncomplementing_human_id, 
         complement_sim, noncomplement_sim, complement_higher)

# Also get alignment baseline data
alignment_data <- data %>%
  filter(method == "Sequence Alignment") %>%
  select(yeast_id, complementing_human_id, noncomplementing_human_id, 
         complement_sim, noncomplement_sim, complement_higher) %>%
  unique()

print(esm_plus_plus_data)
# Prepare data for paired dot plot by converting to long format
esm_plus_plus_long <- esm_plus_plus_data %>%
  mutate(case_id = row_number()) %>%
  pivot_longer(
    cols = c(complement_sim, noncomplement_sim),
    names_to = "similarity_type",
    values_to = "similarity_value"
  ) %>%
  mutate(
    similarity_type = factor(
      similarity_type, 
      levels = c("noncomplement_sim", "complement_sim"),
      labels = c("Non-complementing", "Complementing")
    ),
    # Add label for connecting lines
    correct_prediction = ifelse(complement_higher, "Correct", "Incorrect")
  )

print(esm_plus_plus_long)
print(esm_plus_plus_long$complement_higher)

alignment_long <- alignment_data %>%
  mutate(case_id = row_number()) %>%
  pivot_longer(
    cols = c(complement_sim, noncomplement_sim),
    names_to = "similarity_type",
    values_to = "similarity_value"
  ) %>%
  mutate(
    similarity_type = factor(
      similarity_type, 
      levels = c("noncomplement_sim", "complement_sim"),
      labels = c("Non-complementing", "Complementing")
    ),
    # Add label for connecting lines
    correct_prediction = ifelse(complement_higher, "Correct", "Incorrect")
  )

# Create the paired dot plot function
create_paired_dot_plot <- function(data, title, y_label = "Similarity") {
  # Order data to ensure "Incorrect" (blue) points are plotted last (on top)
  ordered_data <- data %>%
    arrange(correct_prediction == "Incorrect")
  
  p <- ggplot(ordered_data, aes(x = similarity_type, y = similarity_value, color = correct_prediction, group = case_id)) +
    # Add connecting lines colored by whether prediction was correct
    geom_line(alpha = 0.7) +
    # Add dots at each end of the lines
    geom_point(size = 3, alpha = 0.7, aes(shape = similarity_type)) +
   
    scale_color_manual(values = c("Correct" = "#1F77B4", "Incorrect" = "red")) +
    # Labels
    scale_shape_manual(values = c(1,16)) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      color = "Prediction"
    ) +
    theme(
      
      legend.position = "bottom",
      panel.grid.major.x = element_blank()  # Remove vertical gridlines
    )
  
  return(p)
}

# Create two paired dot plots side by side with new colors
paired_plot_esm <- create_paired_dot_plot(
  esm_plus_plus_long, 
  "ESMC 600M\nSWE", 
  "Embedding Similarity"
)

paired_plot_alignment <- create_paired_dot_plot(
  alignment_long, 
  "Sequence Alignment", 
  "Alignment Similarity"
)

# Combine the plots
library(patchwork)
paired_plots <- paired_plot_alignment + paired_plot_esm + theme(legend.position = "none") + plot_layout(ncol = 2)

paired_plots


# Create Figure 4b - ESM++ Large similarities vs sequence identity
esm_plus_plus_with_alignment <- parsed_data %>%
  filter(model == "ESMC-600M",
         layer_config == "last_8",
         embedding_type == "sequence_embeddings_swe",
         similarity_method_parsed == "cosine") %>%
  # Join with alignment data to get sequence identities
  left_join(
    alignment_data,
    by = c("yeast_id", "complementing_human_id", "noncomplementing_human_id"),
    suffix = c("", "_alignment")
  ) %>%
  select(yeast_id, complementing_human_id, noncomplementing_human_id,
         complement_sim, noncomplement_sim, complement_higher,
         complement_sim_alignment, noncomplement_sim_alignment)

# Convert to long format for plotting
esm_plus_plus_vs_alignment_long <- esm_plus_plus_with_alignment %>%
  mutate(case_id = row_number()) %>%
  pivot_longer(
    cols = c(complement_sim, noncomplement_sim),
    names_to = "similarity_type",
    values_to = "embedding_similarity"
  ) %>%
  # Add corresponding alignment similarities
  mutate(
    alignment_similarity = ifelse(
      similarity_type == "complement_sim",
      complement_sim_alignment,
      noncomplement_sim_alignment
    ),
    similarity_type = factor(
      similarity_type,
      levels = c("noncomplement_sim", "complement_sim"),
      labels = c("Non-complementing", "Complementing")
    ),
    correct_prediction = ifelse(complement_higher, "Correct", "Incorrect")
  )

# Create the scatter plot
embVident <- ggplot(esm_plus_plus_vs_alignment_long, 
       aes(x = alignment_similarity, y = embedding_similarity, 
           color = correct_prediction, shape = similarity_type)) +
  geom_point(size = 3, alpha = 0.7) +
  # Add lines connecting pairs
  geom_line(aes(group = case_id), alpha = 0.2) +
  scale_color_manual(values = c("Correct" = "#1F77B4", "Incorrect" = "red")) +
  scale_shape_manual(values = c(1,16)) +
  labs(
    title = "ESM++ Large Embedding\nvs.\nSequence Alignment Similarities",
    x = "Sequence Identity",
    y = "Embedding Similarity",
    color = "Prediction",
    shape = "Pair Type"
  ) +
  theme_cowplot(font_size = 8) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

embVident

# Panel A: keep legend, put it on the right
embperf_labeled <- embperf +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    legend.position = "right"
  )

# Panel B: no legend
paired_plot_alignment_labeled <- paired_plot_alignment +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 8, face = "bold"), legend.position = "none")

# Panel (center): no legend, no label
paired_plot_esm_nolegend <- paired_plot_esm + theme(legend.position = "none")

# Panel C: legend on right, label
embVident_labeled <- embVident +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    legend.position = "right"
  )

# Arrange the panels: 2 rows, 3 columns
fig2 <- (embperf_labeled) /
        (paired_plot_alignment_labeled + paired_plot_esm_nolegend + embVident_labeled +
         plot_layout(ncol = 3, widths = c(1, 1, 1.2))) +
        plot_layout(heights = c(1, 1.2))

# Save the combined plot with more vertical space
ggsave("figures/Fig2_qfo.pdf", fig2, width = 8, height = 6, units = "in")
ggsave("figures/Fig2_qfo.png", fig2, width = 8, height = 6, units = "in", bg = "white")
print("Plot saved as Fig2.pdf")




# Create paralog ranking heatmap
create_paralog_ranking_heatmap <- function() {
  # Read the paralog groups data
  paralog_data <- read_csv("data/paralog_groups_output.csv")%>%
    mutate(model = case_when(model == "ESMplusplus_small" ~ "ESMC-300M",
                             model == "ESMplusplus_large" ~ "ESMC-600M",
                             TRUE ~ model))
  head(paralog_data)
  
  # Filter for last_8 layer configuration, cosine similarity, and SWE embeddings
  filtered_data_align <- paralog_data %>%
    filter(
      method == "Sequence Alignment"
    ) %>%
    select(-embedding_filename, -model, -layer_config, -similarity_method) %>%
    mutate(model = "Sequence Identity") %>%
    unique() 

  filtered_data <- paralog_data %>%
    filter(
      layer_config == "last_8",
      similarity_method == "cosine",
      grepl("swe", embedding_filename, ignore.case = TRUE)
    )
  
  # Count paralogs per yeast gene and create ordering
  paralog_counts <- filtered_data %>%
    group_by(yeast_id) %>%
    summarize(
      num_paralogs = n(),
      .groups = "drop"
    ) %>%
    arrange(num_paralogs) %>%  # Sort by number of paralogs
    pull(yeast_id)  # Get ordered yeast IDs
  
  # Process data for the heatmap
  heatmap_data <- filtered_data %>%
    bind_rows(filtered_data_align) %>%
    # Ensure yeast_id is ordered by number of paralogs
    mutate(
      yeast_id = factor(yeast_id, levels = paralog_counts),
      # Convert complements to a factor for color mapping
      complement_status = ifelse(complements, "Complementing", "Non-complementing"),
      model = factor(model, levels = rev(model_order))
      
    )

  # Create the heatmap using ggplot2
  heatmap_plot <- ggplot(heatmap_data, 
         aes(x = yeast_id, 
             y = similarity_rank, 
             fill = complement_status)) +
    # Fill cells based on complementation status
    geom_tile(color = "white", size = 0.2) +
    # Use colorblind-friendly colors
    scale_fill_manual(values = c("Complementing" = "skyblue3", "Non-complementing" = "red")) +
    # Reverse y-axis so rank 1 is at the top
    scale_y_reverse(breaks = 1:10) +
    # Use minimal theme for cleaner look
    theme_minimal() +
    # Customize appearance
    theme(
      axis.text.x = element_blank(),  # Hide yeast gene IDs
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),   # Remove grid lines
      panel.border = element_rect(fill = NA, color = "grey80"),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    ) +
    # Add labels
    labs(
      x = "Yeast Genes (ordered by number of paralogs)",
      y = "Similarity Rank (1 = Highest Similarity)",
      fill = "Status"
    ) +
    facet_wrap(~model)
  
  # Save the plot
  ggsave("figures/Fig5_Paralog_Ranking_Heatmap.pdf", heatmap_plot, width = 8, height = 6)
  ggsave("figures/Fig5_Paralog_Ranking_Heatmap.png", heatmap_plot, width = 8, height = 6, bg = "white")
  print("Plot saved as Fig5_Paralog_Ranking_Heatmap.pdf")
  
  # Also create a summary plot showing the distribution of complementing paralogs by rank
  summary_data <- heatmap_data %>%
    filter(complements == TRUE) %>%
    group_by(similarity_rank) %>%
    summarize(count = n(), .groups = "drop") %>%
    # Calculate percentage
    mutate(
      total = sum(count),
      percentage = count / total * 100
    )
  
  # Plot the distribution of complementing paralog ranks
  rank_distribution_plot <- ggplot(summary_data, 
       aes(x = similarity_rank, y = percentage)) +
    geom_col(fill = "#4DAF4A") +
    scale_x_continuous(breaks = 1:10) +
    theme_minimal() +
    labs(
      title = "Distribution of Complementing Paralog Ranks",
      x = "Similarity Rank",
      y = "Percentage of Complementing Paralogs"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  
  # Save the summary plot
  ggsave("data/Figure6_Complementing_Paralog_Rank_Distribution.pdf", rank_distribution_plot, width = 6, height = 4)
  print("Plot saved as Figure6_Complementing_Paralog_Rank_Distribution.pdf")
  
  return(list(heatmap_plot = heatmap_plot, rank_distribution_plot = rank_distribution_plot))
}

# Call the function at the end of your script
paralog_ranking_plots <- create_paralog_ranking_heatmap()

paralog_ranking_plots$heatmap_plot
ggsave("figures/Supplemental_Fig2.pdf" , paralog_ranking_plots$heatmap_plot, width = 8, height = 6)
ggsave("figures/Supplemental_Fig2.png" , paralog_ranking_plots$heatmap_plot, width = 8, height = 6, dpi = 600, bg = "white")




