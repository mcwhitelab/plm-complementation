#!/usr/bin/env Rscript

# Libraries
library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(scales)
library(patchwork)
library(ggrastr)

# Set theme
theme_set(theme_cowplot(font_size = 8, rel_small = 1, rel_large = 1, rel_tiny = 1))

# Function to get file paths
get_file_paths <- function() {
  list(
    # Mouse-Human (10090-9606)
    mouse_human = list(
      pairs = "static/orthologies/inp9_10090_9606_prot.pairs",
      emb_sim = "data/similarities/inp9_10090_9606_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_10090_9606_prot.pairs.align.similarities"
    ),
    # Mouse-Yeast (10090-559292)
    mouse_yeast = list(
      pairs = "static/orthologies/inp9_10090_559292_prot.pairs",
      emb_sim = "data/similarities/inp9_10090_559292_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_10090_559292_prot.pairs.align.similarities"
    ),
    # Mouse-Pombe (10090-284812)
    mouse_pombe = list(
      pairs = "static/orthologies/inp9_10090_284812_prot.pairs",
      emb_sim = "data/similarities/inp9_10090_284812_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_10090_284812_prot.pairs.align.similarities"
    ),
    # Human-Yeast (9606-559292)
    yeast_human = list(
      pairs = "static/orthologies/inp9_559292_9606_prot.pairs",
      emb_sim = "data/similarities/inp9_559292_9606_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_559292_9606_prot.pairs.align.similarities"
    ),
    # Human-Pombe (9606-284812)
    pombe_human = list(
      pairs = "static/orthologies/inp9_284812_9606_prot.pairs",
      emb_sim = "data/similarities/inp9_284812_9606_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_284812_9606_prot.pairs.align.similarities"
    ),
    # Pombe-Yeast (284812-559292)
    pombe_yeast = list(
      pairs = "static/orthologies/inp9_284812_559292_prot.pairs",
      emb_sim = "data/similarities/inp9_284812_559292_prot.pairs.sequence_embeddings_swe.esmppl.1l.similarities",
      seq_sim = "data/similarities/inp9_284812_559292_prot.pairs.align.similarities"
    )
  )
}

# Function to process pairs and similarities files
process_similarity_data <- function(pairs_file, emb_sim_file, seq_sim_file, species_label) {
  # Read input files
  pairs <- read_csv(pairs_file)
  emb_sims <- read_csv(emb_sim_file)
  seq_sims <- read_csv(seq_sim_file)
  
  # Identify many-to-many group_ids in this pairs file
  many_many_paralog_group_ids <- pairs %>%
    filter(label == "many-to-many") %>%
    pull(group_id) %>%
    unique()



one_two_paralog_group_ids <- pairs %>%
    filter(label == "one-to-many" | label == "many-to-one") %>%
    group_by(group_id) %>%
        mutate(group_size = n()) %>%
    filter(group_size == 2) %>%
    pull(group_id) %>%
    unique()    

 one_many_paralog_group_ids <- pairs %>%
    filter(label == "one-to-many" | label == "many-to-one") %>%
    filter(!group_id %in% one_two_paralog_group_ids) %>%
    pull(group_id) %>%
    unique()

one_one_ortholog_group_ids <- pairs %>%
    filter(label == "one-to-one") %>%
    pull(group_id) %>%
    unique()

pairs_w_info <- pairs %>%
   mutate(group_type = case_when(
    label != "inparalogs" & group_id %in% one_one_ortholog_group_ids ~ "one-to-one",
    label != "inparalogs" & group_id %in% one_two_paralog_group_ids ~ "one-to-two",
    label != "inparalogs" & group_id %in% one_many_paralog_group_ids ~ "one-to-many",
    label != "inparalogs" & group_id %in% many_many_paralog_group_ids ~ "many-to-many",
    label == "inparalogs" & group_id %in% one_two_paralog_group_ids ~ "two-copy",
    label == "inparalogs" & group_id %in% one_many_paralog_group_ids ~ "many-copy (1-many)",
    label == "inparalogs" & group_id %in% many_many_paralog_group_ids ~ "many-copy (many-many)",
    TRUE ~ "other"
   ),
    color_category = case_when(
        label == "one-to-one" ~ "LDO",
        is_ldo == TRUE  & label != "inparalogs"~ "LDO",
        is_ldo == FALSE & label != "inparalogs" ~ "non-LDO",
        TRUE ~ "other")
    ) %>%
        ungroup()
 
  pairs_w_info %>% write_csv("/groups/clairemcwhite/claire_workspace/github/seqsim_project/test3.csv")
  # Combine all data
  combined_data <- pairs_w_info %>%
    left_join(emb_sims, by = c("protein1_id", "protein2_id")) %>%
    left_join(seq_sims %>% select(protein1_id, protein2_id, seq_identity = similarity),
              by = c("protein1_id", "protein2_id"))
  
  combined_data$species_label <- species_label
  
  return(combined_data)
}

# Function to load or create all_data
get_all_data <- function(force_reprocess = TRUE) {
  if (!force_reprocess && file.exists("data/all_data.csv")) {
    # Load existing data
    all_data <- read_csv("data/all_data.csv")
    return(all_data)
  }
  
  # Get file paths
  file_paths <- get_file_paths()

# Process all datasets
all_data <- bind_rows(
    process_similarity_data(file_paths$mouse_human$pairs, 
                          file_paths$mouse_human$emb_sim, 
                          file_paths$mouse_human$seq_sim, 
                          "Mouse-Human"),
    process_similarity_data(file_paths$mouse_yeast$pairs, 
                          file_paths$mouse_yeast$emb_sim, 
                          file_paths$mouse_yeast$seq_sim, 
                          "Mouse-S. cerevisiae"),
    process_similarity_data(file_paths$mouse_pombe$pairs, 
                          file_paths$mouse_pombe$emb_sim, 
                          file_paths$mouse_pombe$seq_sim, 
                          "Mouse-S. pombe"),
    process_similarity_data(file_paths$yeast_human$pairs, 
                          file_paths$yeast_human$emb_sim, 
                          file_paths$yeast_human$seq_sim, 
                          "S. cerevisiae-Human"),
    process_similarity_data(file_paths$pombe_human$pairs, 
                          file_paths$pombe_human$emb_sim, 
                          file_paths$pombe_human$seq_sim, 
                          "S. pombe-Human"),
    process_similarity_data(file_paths$pombe_yeast$pairs, 
                          file_paths$pombe_yeast$emb_sim, 
                          file_paths$pombe_yeast$seq_sim, 
                          "S. pombe-S. cerevisiae")
  )
  
  # Create data directory if it doesn't exist
  dir.create("data", showWarnings = FALSE)
  dir.create("figures", showWarnings = FALSE)
  # Save the processed data
  write_csv(all_data, "data/all_data.csv")
  
  return(all_data)
}

# After loading all_data, create the annotated version immediately
all_data <- get_all_data(force_reprocess = FALSE)
print(unique(all_data$species_label))



all_data %>% write_csv("test7.csv")

# Function to create LDO analysis data for a species pair
create_ldo_data <- function(data, species_pair) {
  data %>%
    filter(
      species_label == species_pair,
      label != "inparalogs",
      group_type %in% c("one-to-one", "one-to-two", "one-to-many")
    ) %>%
    group_by(group_id) %>%
    mutate(
      max_seq_identity = max(seq_identity),
      color_category = ifelse(seq_identity == max_seq_identity, "LDO", "non-LDO"),
      max_similarity = max(similarity),
      ldo_similarity = first(similarity[seq_identity == max_seq_identity]),
      expected = ifelse(any(seq_identity == max_seq_identity & similarity == max_similarity),
        "Expected (LDO highest)",
        "Unexpected (LDO not highest)"
      )
    ) %>%
    ungroup()
}

print("HEY")

# Then modify all downstream analyses to use all_data_annot instead of all_data
# For example, the plot function would become:
plot_ldo_analysis <- function(data) {
  # Now we can use the pre-computed annotations
  plot_data <- data %>%
    filter(color_category %in% c("LDO", "non-LDO")) %>%
    mutate(group_type = factor(group_type, 
                              levels = c("one-to-one", "one-to-two", "one-to-many")))
  
  ggplot(plot_data, aes(x = seq_identity, y = similarity, 
                         group = group_id)) +
    geom_point(alpha = 0.5, size = 0.5, aes(color = color_category)) +
    geom_line(alpha = 0.1) +
    facet_grid(expected ~ group_type) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "LDO" = "#D7191C",      # red
      "non-LDO" = "#2C7BB6"   # blue
    )) +
    labs(x = "Sequence Identity",
         y = "Embedding Similarity",
         color = "Category") +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())  # Remove all gridlines
}


# Create and plot LDO analysis for both species pairs
species_pairs <- c("Mouse-Human", "S. cerevisiae-Human")
for (species in species_pairs) {
  ortholog_data <- create_ldo_data(all_data, species)
  p <- plot_ldo_analysis(ortholog_data)
  filename <- paste0("figures/ldo_analysis_", tolower(gsub("[. ]", "_", species)), ".pdf")
  ggsave(filename, p, width = 6, height = 4)
}

# Analyze unexpected groups with significant differences
unexpected_groups <- ortholog_data %>%
  filter(expected == "Unexpected (LDO not highest)") %>%
  group_by(group_id) %>%
  summarize(
    group_type = first(group_type),
    protein_pairs = paste(paste(protein1_id, protein2_id, sep="-"), collapse = ";"),
    all_seq_identities = paste(seq_identity, collapse = ";"),
    all_similarities = paste(similarity, collapse = ";"),
    all_categories = paste(color_category, collapse = ";"),
    similarity_diff = max(similarity) - similarity[color_category == "LDO"]
  ) %>%
  filter(similarity_diff > 0.01) %>%
  arrange(desc(similarity_diff))

print("Unexpected groups with >0.01 embedding similarity difference:")
print(unexpected_groups)

# Save to file
write.table(unexpected_groups, "data/unexpected_groups.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# And the summary statistics would use all_data_annot directly
expectation_stats <- ortholog_data %>%
  group_by(group_type, expected) %>%
  summarize(
    count = n_distinct(group_id),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_emb_similarity = mean(similarity, na.rm = TRUE)
  ) %>%
  arrange(group_type, expected)

print("Expectation Analysis Summary:")
print(expectation_stats)


# Print summary statistics for group analysis
summary_stats <- all_data %>%
  group_by(group_type, color_category) %>%
  summarize(
    n = n(),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_emb_similarity = mean(similarity, na.rm = TRUE)
  ) %>%
  arrange(group_type, color_category)

print("Group Analysis Summary:")
print(summary_stats)

# Print summary statistics for LDO analysis
summary_stats <- all_data %>%
  filter(
    label != "inparalogs",
    group_type %in% c("one-to-one", "one-to-two", "one-to-many")
  ) %>%
  group_by(group_type, color_category) %>%
  summarize(
    n = n(),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_emb_similarity = mean(similarity, na.rm = TRUE)
  ) %>%
  arrange(group_type, color_category)

print("LDO Analysis Summary:")
print(summary_stats)


# Now, update paralog_ortholog_data to include these two types
paralog_ortholog_data <- all_data %>%
  filter(
    (group_type == "one-to-one") |  # 1-1 orthologs
    (label == "inparalogs" & group_type == "two-copy") |    # two-copy inparalogs
    (label == "inparalogs" & group_type == "many-copy (1-many)") |
    (label == "inparalogs" & group_type == "many-copy (many-many)") # only annotated many-copy
  ) %>%
  mutate(
    comparison_type = factor(case_when(
      group_type == "one-to-one" & label != "inparalogs" ~ "one-to-one orthologs",
      group_type == "two-copy" & label == "inparalogs" ~ "two-copy inparalogs",
      group_type == "many-copy (1-many)" & label == "inparalogs" ~ "many-copy inparalogs\n(1-to-many)",
      group_type == "many-copy (many-many)" & label == "inparalogs" ~ "many-copy inparalogs\n(many-to-many)"
    ), levels = c(
      "one-to-one orthologs",
      "two-copy inparalogs",
      "many-copy inparalogs\n(1-to-many)",
      "many-copy inparalogs\n(many-to-many)"
    ))
  )

# Function to plot ortholog vs paralog comparison
plot_ortholog_paralog_comparison <- function(data) {
  # Define the color palette
  custom_colors <- c(
    "#2f4f4f", "#556b2f", "#8b4513", "#228b22", "#708090", "#8b0000", "#808000",
    "#483d8b", "#3cb371", "#bc8f8f", "#008080", "#cd853f", "#4682b4", "#9acd32",
    "#cd5c5c", "#4b0082", "#32cd32", "#daa520", "#8fbc8f", "#8b008b", "#b03060",
    "#9932cc", "#ff4500", "#00ced1", "#ff8c00", "#ffd700", "#6a5acd", "#ffff00",
    "#0000cd", "#00ff00", "#00ff7f", "#e9967a", "#dc143c", "#00bfff", "#0000ff",
    "#a020f0", "#adff2f", "#ff6347", "#da70d6", "#d8bfd8", "#ff00ff", "#f0e68c",
    "#6495ed", "#dda0dd", "#add8e6", "#ff1493", "#98fb98", "#7fffd4", "#ff69b4",
    "#ffe4c4"
  )
  
  # Get unique group IDs for paralogs only and shuffle them randomly
  set.seed(42)  # for reproducibility
  paralog_data <- data %>% filter(comparison_type != "one-to-one orthologs")
  group_ids <- sample(unique(paralog_data$group_id))
  n_groups <- length(group_ids)
  
  # Repeat the color palette as needed to match the number of groups
  color_palette <- rep(custom_colors, ceiling(n_groups/length(custom_colors)))[1:n_groups]
  
  # Create named vector of colors with shuffled group IDs
  group_colors <- setNames(color_palette, group_ids)
  
  # Split the data into orthologs and paralogs
  orthologs <- data %>% filter(comparison_type == "one-to-one orthologs")
  paralogs <- data %>% filter(comparison_type != "one-to-one orthologs")
  
  # Create base plot
  p <- ggplot() +
    # Add grey points for orthologs
    geom_point(data = orthologs, 
               aes(x = seq_identity, y = similarity),
               color = "grey70",
               alpha = 0.2, 
               size = 0.1) +
    # Add colored points for paralogs
    geom_point(data = paralogs,
               #aes(x = seq_identity, y = similarity, color = as.character(group_id)),
               aes(x = seq_identity, y = similarity),
               color = "grey70",
               alpha = 0.2,
               size = 0.1) +
    facet_wrap(species_label~comparison_type, ncol = 4, drop = FALSE) +
    scale_y_continuous(limits = c(0.9, 1)) +  # Set y-axis minimum to 0.9
    scale_color_manual(values = group_colors) +
    labs(x = "Sequence Identity",
         y = "Embedding Similarity") +
    theme(legend.position = "none",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())  # Remove all gridlines
  p_raster <- rasterize(p)
  return(p_raster)
}

# Generate and save the ortholog-paralog comparison plot
paralog_ortholog_data <- paralog_ortholog_data %>%
  filter(species_label == "Mouse-Human" | species_label == "S. cerevisiae-Human") %>%
  mutate(species_label = factor(species_label, 
                              levels = c("S. cerevisiae-Human", "Mouse-Human")))
p_comparison <- plot_ortholog_paralog_comparison(paralog_ortholog_data)
ggsave("figures/ortholog_paralog_comparison.pdf", 
       p_comparison, width = 4, height = 3)

paralog_ortholog_data %>% 
group_by(group_id) %>% 
  mutate(n = n()) %>%
  arrange(desc(n)) %>% 
  write_csv("ortholog_paralog_comparison.csv") 

# Calculate summary statistics for this comparison
comparison_stats <- paralog_ortholog_data %>%
  group_by(comparison_type) %>%
  summarize(
    n = n(),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_emb_similarity = mean(similarity, na.rm = TRUE),
    median_seq_identity = median(seq_identity, na.rm = TRUE),
    median_emb_similarity = median(similarity, na.rm = TRUE)
  )

print("Ortholog-Paralog Comparison Summary:")
print(comparison_stats)

# Save comparison data
write.table(paralog_ortholog_data, 
            "data/ortholog_paralog_comparison.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Find cases with high sequence identity but low embedding similarity
unusual_similarity_cases <- paralog_ortholog_data %>%
  filter(
    species_label == "Mouse-Human",
    similarity < 0.95,
    seq_identity > 0.6
  ) %>%
  select(
    protein1_id,
    protein2_id,
    comparison_type,
    seq_identity,
    similarity
  ) %>%
  arrange(similarity)

print("Cases with embedding similarity < 0.95 and sequence identity > 0.6 in Mouse-Human:")
print(unusual_similarity_cases)

# Save to file
write.table(unusual_similarity_cases,
            "data/unusual_similarity_cases.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Find many-copy orthologs with low embedding similarity in Human-Yeast
human_yeast_many_copy <- paralog_ortholog_data %>%
  filter(
    species_label == "Human-S. cerevisiae",
    comparison_type == "many-copy inparalogs",
    similarity < 0.925
  ) %>%
  select(
    protein1_id,
    protein2_id,
    comparison_type,
    seq_identity,
    similarity
  ) %>%
  arrange(similarity)

print("Human-Yeast many-copy cases with embedding similarity < 0.925:")
print(human_yeast_many_copy)

# Save to file
write.table(human_yeast_many_copy,
            "data/human_yeast_many_copy_low_sim.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Find one-to-one orthologs with high sequence identity but low embedding similarity in Human-Yeast
human_yeast_one_to_one <- paralog_ortholog_data %>%
  filter(
    species_label == "Human-S. cerevisiae",
    comparison_type == "one-to-one orthologs",
    similarity < 0.93,
    seq_identity > 0.5
  ) %>%
  select(
    protein1_id,
    protein2_id,
    comparison_type,
    seq_identity,
    similarity
  ) %>%
  arrange(similarity)

print("Human-Yeast one-to-one cases with embedding similarity < 0.93 and sequence identity > 0.5:")
print(human_yeast_one_to_one)

# Save to file
write.table(human_yeast_one_to_one,
            "data/human_yeast_one_to_one_unusual.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Find many-copy cases with low sequence similarity in Mouse-Human
mouse_human_many_copy_low_seq <- paralog_ortholog_data %>%
  filter(
    species_label == "Mouse-Human",
    comparison_type == "many-copy inparalogs",
    seq_identity < 0.3
  ) %>%
  select(
    protein1_id,
    protein2_id,
    comparison_type,
    seq_identity,
    similarity
  ) %>%
  arrange(seq_identity)

print("Mouse-Human many-copy cases with sequence identity < 0.3:")
print(mouse_human_many_copy_low_seq)

# Save to file
write.table(mouse_human_many_copy_low_seq,
            "data/mouse_human_many_copy_low_seq.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Create data for ortholog conjecture plot
ortholog_conjecture_data <- all_data %>%
  filter(
    (group_type == "one-to-one") |  # 1-1 orthologs
    (label == "inparalogs" & group_type %in% c("two-copy", "many-copy (1-many)", "many-copy (many-many)"))  # paralogs
  ) %>%
  mutate(
    comparison_type = factor(case_when(
      group_type == "one-to-one" & label != "inparalogs" ~ "one-to-one\northologs\n",
      group_type == "two-copy" & label == "inparalogs" ~ "two-copy\ninparalogs\n",
      group_type == "many-copy (1-many)" & label == "inparalogs" ~ "many-copy\ninparalogs\n(1-to-many)",
      group_type == "many-copy (many-many)" & label == "inparalogs" ~ "many-copy\ninparalogs\n(many-to-many)"
    ), levels = c("one-to-one\northologs\n", "two-copy\ninparalogs\n", "many-copy\ninparalogs\n(1-to-many)", "many-copy\ninparalogs\n(many-to-many)")),
    # Create sequence identity bins
    seq_bin = factor(cut(seq_identity, 
                  breaks = seq(0, 1, by = 0.1),
                  labels = paste0("(", seq(0, 0.9, by = 0.1), ",", seq(0.1, 1, by = 0.1), "]"),
                  include.lowest = TRUE),
                  levels = paste0("(", seq(0, 0.9, by = 0.1), ",", seq(0.1, 1, by = 0.1), "]"))
  )

# Calculate mean and standard error for each bin
ortholog_conjecture_summary <- ortholog_conjecture_data %>%
  group_by(comparison_type, seq_bin) %>%
  summarise(
    mean_similarity = mean(similarity, na.rm = TRUE),
    se_similarity = sd(similarity, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(seq_bin))

# Calculate species-specific summary
ortholog_conjecture_species_summary <- ortholog_conjecture_data %>%
  group_by(species_label, comparison_type, seq_bin) %>%
  summarise(
    mean_similarity = mean(similarity, na.rm = TRUE),
    se_similarity = sd(similarity, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(seq_bin))

# Create the ortholog conjecture plot
plot_ortholog_conjecture <- function(data) {
  ggplot(data, aes(x = seq_bin, y = mean_similarity, 
                   color = comparison_type, 
                   group = comparison_type,
                   shape = comparison_type)) +  
    geom_line(alpha = 0.7) +
    geom_point(size = 2) +  
    geom_errorbar(aes(ymin = mean_similarity - se_similarity, 
                      ymax = mean_similarity + se_similarity),
                  width = 0.2) +
    scale_color_manual(values = c(
      "one-to-one\northologs\n" = "green",
      "two-copy\ninparalogs\n" = "blue",
      "many-copy\ninparalogs\n(1-to-many)" = "orange",
      "many-copy\ninparalogs\n(many-to-many)" = "purple"
    )) +
    scale_shape_manual(values = c(
      "one-to-one\northologs\n" = 16,  # filled circle
      "two-copy\ninparalogs\n" = 17,   # filled triangle
      "many-copy\ninparalogs\n(1-to-many)" = 15,    # filled square
      "many-copy\ninparalogs\n(many-to-many)" = 18  # filled diamond
    )) +
    labs(x = "Sequence Identity",
         y = "Embedding Similarity",
         color = NULL,
         shape = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.spacing.y = unit(0.1, "cm"),
      legend.box.just = "left",
      legend.justification = "left",
      legend.box.margin = margin(0, 0, 0, -0.5, "cm"),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank())
}

# Create species-specific ortholog conjecture plot
plot_ortholog_conjecture_by_species <- function(data) {
  ggplot(data, aes(x = seq_bin, y = mean_similarity, 
                   color = comparison_type, group = comparison_type)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_similarity - se_similarity, 
                      ymax = mean_similarity + se_similarity),
                  width = 0.2) +
    facet_wrap(~species_label, ncol = 2) +
    scale_color_manual(values = c(
      "one-to-one\northologs\n" = "red",
      "two-copy\ninparalogs\n" = "blue",
      "many-copy\ninparalogs\n(1-to-many)" = "orange",
      "many-copy\ninparalogs\n(many-to-many)" = "purple"
    )) +
    labs(x = "Sequence Identity",
         y = "Embedding Similarity",
         color = "Comparison Type") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(size = 8),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank())  # Remove all gridlines
}

# Create a blank panel for the diagram
blank_panel <- ggplot() + 
  theme_void() +
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 8, face = "bold"))


print(ortholog_conjecture_summary)
ortholog_conjecture_summary %>% write_csv("ortholog_conjecture_summary.csv")
# Modify the other plots with new panel labels
p_conjecture <- plot_ortholog_conjecture(ortholog_conjecture_summary) +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 8, face = "bold"))

p_comparison <- plot_ortholog_paralog_comparison(paralog_ortholog_data) +
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 8, face = "bold"))

# Combine plots using patchwork
# First combine A and B side by side, then C below
Fig4 <- (blank_panel + p_conjecture + plot_layout(widths = c(1.5, 1))) / p_comparison +
  plot_layout(heights = c(1, 2.2)) 

ggsave("figures/Fig4_qfo.pdf", Fig4, width = 8, height = 6)
ggsave("figures/Fig4_qfo.png", Fig4, width = 8, height = 6, dpi = 300, bg = "white")


# Still save individual plots as before
ggsave("figures/ortholog_conjecture.pdf", p_conjecture, width = 4, height = 3)
ggsave("figures/ortholog_paralog_comparison.pdf", p_comparison, width = 4, height = 3)


# Print summary statistics for each comparison type
conjecture_stats <- ortholog_conjecture_data %>%
  group_by(comparison_type) %>%
  summarise(
    n = n(),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_similarity = mean(similarity, na.rm = TRUE)
  )

print("Ortholog Conjecture Analysis Summary:")
print(conjecture_stats)

# Print summary statistics for each species pair and comparison type
conjecture_species_stats <- ortholog_conjecture_data %>%
  group_by(species_label, comparison_type) %>%
  summarise(
    n = n(),
    mean_seq_identity = mean(seq_identity, na.rm = TRUE),
    mean_similarity = mean(similarity, na.rm = TRUE)
  )

print("Ortholog Conjecture Analysis Summary by Species:")
print(conjecture_species_stats)

print("here")
# Create combined LDO analysis plot for both species pairs
plot_combined_ldo_analysis <- function(mouse_human_data, yeast_human_data) {
  # Filter and prepare mouse-human data
  mh_data <- mouse_human_data %>%
    filter(group_type %in% c("one-to-one", "one-to-two")) %>%
    mutate(species = "Mouse-Human",
           group_type = factor(group_type, levels = c("one-to-one", "one-to-two")))
  
  # Filter and prepare yeast-human data
  yh_data <- yeast_human_data %>%
    filter(group_type %in% c("one-to-one", "one-to-two")) %>%
    mutate(species = "S. cerevisiae-Human",
           group_type = factor(group_type, levels = c("one-to-one", "one-to-two")))
  
  # Combine the data
  combined_data <- bind_rows(mh_data, yh_data) %>%
    mutate(
      species = factor(species, levels = c("Mouse-Human", "S. cerevisiae-Human")),
      # Create a new column that combines group type and expected status
      plot_group = case_when(
        group_type == "one-to-one" ~ "one-to-one",
        group_type == "one-to-two" & expected == "Expected (LDO highest)" ~ "one-to-two\nConcordant",
        group_type == "one-to-two" & expected == "Unexpected (LDO not highest)" ~ "one-to-two\nDiscordant",
        TRUE ~ NA_character_
      ),
      plot_group = factor(plot_group, 
                         levels = c("one-to-one",
                                  "one-to-two\nConcordant",
                                  "one-to-two\nDiscordant")),
      # Create a new color category that combines group type and LDO status
      plot_color = case_when(
        group_type == "one-to-one" ~ "one-to-one",
        color_category == "LDO" ~ "LDO",
        color_category == "non-LDO" ~ "non-LDO",
        TRUE ~ "other"
      )
    ) %>%
    filter(!is.na(plot_group))  # Remove any NA values
  
  # Calculate counts for each facet
  facet_counts <- combined_data %>%
    group_by(species, plot_group) %>%
    summarise(n_groups = n_distinct(group_id),
             .groups = 'drop') %>%
    # Add dummy x and y coordinates for the text
    mutate(x = 0.15,  # This will be in log10 scale
           y = 1)     # This will be in log10 scale
  print(facet_counts)
  
  # Create the plot
  ggplot(combined_data, 
         aes(x = seq_identity, y = similarity)) +
    geom_point(aes(color = plot_color, group = group_id), alpha = 0.5, size = 0.5) +
    geom_line(aes(group = group_id), alpha = 0.1) +
    # Add count labels
    geom_text(data = facet_counts,
              aes(x = x, y = y, label = paste0("n=", n_groups)),
              hjust = 0, vjust = 1,
              size = 2.5) +
    facet_grid(species ~ plot_group, drop = TRUE) +  # Add drop=TRUE to remove empty facets
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(
      "one-to-one" = "grey80",
      "LDO" = "#D7191C",      # red
      "non-LDO" = "#2C7BB6"   # blue
    )) +
    labs(x = "Sequence Identity",
         y = "Embedding Similarity",
         color = "Category") +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8),
          strip.background = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())
}

# Create LDO analysis data for both species pairs
mouse_human_ldo <- create_ldo_data(all_data, "Mouse-Human")
yeast_human_ldo <- create_ldo_data(all_data, "S. cerevisiae-Human")


# Create LDO analysis plot with panel label
p_ldo_both <- plot_combined_ldo_analysis(mouse_human_ldo, yeast_human_ldo) +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 8, face = "bold"))

p_ldo_both


# Create blank panel for new Figure 3
blank_panel_fig3 <- ggplot() + 
  theme_void() +
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 8, face = "bold"))


# Combine plots for Figure 3
Fig3 <- blank_panel_fig3 / rasterize(p_ldo_both) +
  plot_layout(heights = c(1, 2))

ggsave("figures/Fig3_qfo.pdf", Fig3, width = 8, height = 6)
ggsave("figures/Fig3_qfo.png", Fig3, width = 8, height = 6, dpi = 300, bg = "white")
