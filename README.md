# Three_Circular_Plots_Readable
# R script for Three Circular Plots
library(tidyverse)
library(ggalluvial)
library(patchwork)

kegg <- read.delim("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/kegg_enrichment.txt",
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)

print(head(kegg))
print(colnames(kegg))

top_n <- 20  
kegg_top <- kegg %>%
  arrange(Enrichment.FDR) %>%
  slice_head(n = top_n)

kegg_long <- kegg_top %>%
  separate_rows(Genes, sep = "\\s+") %>%      
  rename(Gene = Genes) %>%
  mutate(Gene = str_trim(Gene)) %>%          
  filter(Gene != "" & !is.na(Gene))           

p_alluvial <- ggplot(kegg_long,
                     aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
  geom_alluvium(aes(fill = Pathway), alpha = 0.7, width = 1/6) +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(0.15, 0.15)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11, face = "bold")
  ) +
  labs(x = NULL, y = NULL)

p_dot <- kegg %>%
  mutate(
    GeneRatio = nGenes / Pathway.Genes,
    logFDR = -log10(Enrichment.FDR)
  ) %>%
  arrange(desc(GeneRatio)) %>%
  ggplot(aes(x = GeneRatio,
             y = reorder(Pathway, GeneRatio),
             size = nGenes,
             color = logFDR)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#2166ac", high = "#b2182b", 
                       name = expression(-log[10](FDR))) +
  scale_size_continuous(range = c(4, 9), 
                        name = "Gene Count",
                        breaks = unique(sort(kegg$nGenes)),  
                        labels = unique(sort(kegg$nGenes)))   
theme_minimal(base_size = 12) +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

final_plot <- p_alluvial + p_dot +
  plot_layout(widths = c(2.8, 1)) +
  plot_annotation(
    title = "KEGG Pathway Enrichment Analysis",
    subtitle = "Left: Alluvial plot of genes to pathways | Right: Dotplot of enrichment significance",
    caption = "Produced with ggplot2, ggalluvial & patchwork"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5))

print(final_plot)

ggsave("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/KEGG_Alluvial_Dotplot.pdf",
       plot = final_plot,
       width = 18, height = 11, dpi = 300, device = cairo_pdf)

ggsave("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/KEGG_Alluvial_Dotplot.png",
       plot = final_plot,
       width = 18, height = 11, dpi = 300)

#####################################################
library(tidyverse)
library(ggalluvial)
library(patchwork)

# Read file
kegg <- read.delim("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/GO biological process.txt",
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)

# Rename column to Pathway for compatibility
colnames(kegg)[colnames(kegg) == "GO.biological.process"] <- "Pathway"

# Select top 20 pathways based on FDR
top_n <- 20  
kegg_top <- kegg %>%
  arrange(Enrichment.FDR) %>%
  slice_head(n = top_n)

# Prepare data for alluvial plot (split genes)
kegg_long <- kegg_top %>%
  separate_rows(Genes, sep = "\\s+") %>%      
  rename(Gene = Genes) %>%
  mutate(Gene = str_trim(Gene)) %>%          
  filter(Gene != "" & !is.na(Gene))

# Alluvial plot (left)
p_alluvial <- ggplot(kegg_long,
                     aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
  geom_alluvium(aes(fill = Pathway), alpha = 0.7, width = 1/6) +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(0.15, 0.15)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11, face = "bold")
  ) +
  labs(x = NULL, y = NULL)

# Dot plot (right)
p_dot <- kegg_top %>%
  mutate(
    GeneRatio = nGenes / Pathway.Genes,
    logFDR = -log10(Enrichment.FDR)
  ) %>%
  arrange(desc(GeneRatio)) %>%
  ggplot(aes(x = GeneRatio,
             y = reorder(Pathway, GeneRatio),
             size = nGenes,
             color = logFDR)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#2166ac", high = "#b2182b", 
                       name = expression(-log[10](FDR))) +
  scale_size_continuous(range = c(4, 9), 
                        name = "Gene Count") +
  theme_minimal(base_size = 12) +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Combine two plots
final_plot <- p_alluvial + p_dot +
  plot_layout(widths = c(2.8, 1)) +
  plot_annotation(
    title = "GO Biological Process Enrichment Analysis",
    subtitle = "Left: Alluvial plot of genes to pathways | Right: Dotplot of enrichment significance",
    caption = "Produced with ggplot2, ggalluvial & patchwork"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5))

# Display plot
print(final_plot)

# Save plot
ggsave("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/GO biological process_Dotplot.png",
       plot = final_plot,
       width = 18, height = 11, dpi = 300)

################################################################

library(tidyverse)
library(ggalluvial)
library(patchwork)

# Read file
reactome <- read.delim("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/Reactome pathway.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

# Rename columns for convenience
reactome <- reactome %>%
  rename(
    Pathway = term.description,                  # Pathway name
    nGenes = observed.gene.count,                # Observed gene count
    Background = background.gene.count,          # Background gene count
    FDR = false.discovery.rate,                   # FDR
    Genes = matching.proteins.in.your.network..labels.  # Genes
  )

# Sort by FDR and select all (since only 3 pathways)
reactome_sorted <- reactome %>%
  arrange(FDR)

# Prepare data for alluvial plot (split genes)
reactome_long <- reactome_sorted %>%
  separate_rows(Genes, sep = ",") %>%
  rename(Gene = Genes) %>%
  mutate(Gene = str_trim(Gene)) %>%
  filter(Gene != "" & !is.na(Gene))

# Alluvial plot (left)
p_alluvial <- ggplot(reactome_long,
                     aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
  geom_alluvium(aes(fill = Pathway), alpha = 0.7, width = 1/6) +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(0.15, 0.15)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 11, face = "bold")
  ) +
  labs(x = NULL, y = NULL)

# Dot plot (right)
p_dot <- reactome_sorted %>%
  mutate(
    GeneRatio = nGenes / Background,
    logFDR = -log10(FDR)
  ) %>%
  arrange(desc(GeneRatio)) %>%
  ggplot(aes(x = GeneRatio,
             y = reorder(Pathway, GeneRatio),
             size = nGenes,
             color = logFDR)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#2166ac", high = "#b2182b",
                       name = expression(-log[10](FDR))) +
  scale_size_continuous(range = c(4, 9),
                        name = "Gene Count") +
  theme_minimal(base_size = 12) +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Combine two plots
final_plot <- p_alluvial + p_dot +
  plot_layout(widths = c(2.8, 1)) +
  plot_annotation(
    title = "Reactome Pathway Enrichment Analysis",
    subtitle = "Left: Alluvial plot of genes to pathways | Right: Dotplot of enrichment significance",
    caption = "Produced with ggplot2, ggalluvial & patchwork"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5))

# Display plot
print(final_plot)

# Save plot
ggsave("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis/Reactome_Pathway_Dotplot_Alluvial.png",
       plot = final_plot,
       width = 18, height = 8, dpi = 300)

#=======================================Improved Circular Chord Diagram=====================================
library(tidyverse)
library(ggalluvial)
library(patchwork)

create_alluvial_dotplot <- function(file_path, 
                                    output_name, 
                                    title_text,
                                    pathway_col_name = NULL,
                                    genes_sep = "\\s+",
                                    top_n = 20,
                                    width = 18, height = 11) {
  
  # ---- read file ----
  if (grepl("\\.csv$", file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # ---- manually rename columns based on known patterns ----
  col_mapping <- list(
    Pathway = c("Pathway", "GO.biological.process", "term.description", "term description"),
    Enrichment.FDR = c("Enrichment.FDR", "Enrichment FDR", "FDR", "false.discovery.rate", "false discovery rate"),
    nGenes = c("nGenes", "observed.gene.count", "observed gene count"),
    Pathway.Genes = c("Pathway.Genes", "Pathway Genes", "Background", "background.gene.count", "background gene count"),
    Genes = c("Genes", "matching.proteins.in.your.network..labels.", "matching proteins in your network (labels)")
  )
  
  # If user specified a custom pathway column name, add it to the list
  if (!is.null(pathway_col_name)) {
    col_mapping$Pathway <- c(col_mapping$Pathway, pathway_col_name)
  }
  
  # Apply renaming
  for (std_name in names(col_mapping)) {
    for (orig_name in col_mapping[[std_name]]) {
      if (orig_name %in% names(df)) {
        names(df)[names(df) == orig_name] <- std_name
        break
      }
    }
  }
  
  # After renaming, check essential columns
  required <- c("Pathway", "Enrichment.FDR", "nGenes", "Pathway.Genes", "Genes")
  missing <- required[!required %in% names(df)]
  if (length(missing) > 0) {
    message("Available columns: ", paste(names(df), collapse=", "))
    stop("Missing required columns: ", paste(missing, collapse=", "))
  }
  
  # ---- data cleaning ----
  df <- df %>% filter(!is.na(Genes), Genes != "", trimws(Genes) != "")
  df$Enrichment.FDR <- as.numeric(df$Enrichment.FDR)
  df$nGenes <- as.numeric(df$nGenes)
  df$Pathway.Genes <- as.numeric(df$Pathway.Genes)
  
  # Select top_n by lowest FDR
  df_top <- df %>% arrange(Enrichment.FDR) %>% slice_head(n = top_n)
  if (nrow(df_top) == 0) {
    message("⚠️ No valid data in ", file_path)
    return(NULL)
  }
  
  # Split genes
  df_long <- df_top %>%
    separate_rows(Genes, sep = genes_sep) %>%
    rename(Gene = Genes) %>%
    mutate(Gene = trimws(Gene)) %>%
    filter(Gene != "" & !is.na(Gene))
  
  if (nrow(df_long) == 0) {
    message("⚠️ No gene-pathway pairs after splitting")
    return(NULL)
  }
  
  # ---- alluvial plot ----
  p_alluvial <- ggplot(df_long, aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
    geom_alluvium(aes(fill = Pathway), alpha = 0.7, width = 1/6) +
    geom_stratum(width = 1/6, fill = "grey90", color = "grey50") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(0.15, 0.15)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 10, face = "bold")
    ) +
    labs(x = NULL, y = NULL)
  
  # ---- dot plot ----
  df_dot <- df_top %>%
    mutate(
      GeneRatio = nGenes / Pathway.Genes,
      logFDR = -log10(Enrichment.FDR)
    ) %>%
    arrange(desc(GeneRatio))
  
  p_dot <- ggplot(df_dot, aes(x = GeneRatio, y = reorder(Pathway, GeneRatio),
                              size = nGenes, color = logFDR)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#2166ac", high = "#b2182b", 
                         name = expression(-log[10](FDR))) +
    scale_size_continuous(range = c(4, 9), name = "Gene Count") +
    theme_minimal(base_size = 12) +
    labs(x = "Gene Ratio", y = NULL) +
    theme(axis.text.y = element_text(size = 9),
          legend.position = "right")
  
  # ---- combine ----
  final_plot <- p_alluvial + p_dot +
    plot_layout(widths = c(2.5, 1.2)) +
    plot_annotation(
      title = title_text,
      subtitle = "Left: Alluvial plot of genes to pathways | Right: Dotplot",
      caption = "ggplot2, ggalluvial & patchwork"
    ) &
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5))
  
  # Save at 300 dpi
  output_png <- paste0(output_name, "_Alluvial_Dotplot.png")
  ggsave(filename = output_png, plot = final_plot,
         width = width, height = height, dpi = 300)
  
  message("✅ Saved: ", output_png)
  return(final_plot)
}

setwd("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis")

# 1) KEGG
create_alluvial_dotplot("kegg_enrichment.txt", "KEGG", "KEGG Pathway Enrichment Analysis", top_n = 10)

# 2) GO biological process
create_alluvial_dotplot("GO biological process.txt", "GO_BP", "GO Biological Process Enrichment Analysis", top_n = 10)

# 3) Reactome
create_alluvial_dotplot("Reactome pathway.txt", "Reactome", "Reactome Pathway Enrichment Analysis", genes_sep = ",", top_n = 10)

#==============================Circ=================================================

create_circular_plot <- function(file_path, output_name, title_text, 
                                 top_n = 12, 
                                 width = 12, height = 12,
                                 genes_sep = "\\s+") {
  
  # ---- read file ----
  if (grepl("\\.csv$", file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  # ---- standardise column names ----
  col_mapping <- list(
    Pathway = c("Pathway", "GO.biological.process", "term.description", "term description"),
    Enrichment.FDR = c("Enrichment.FDR", "Enrichment FDR", "FDR", "false.discovery.rate", "false discovery rate"),
    nGenes = c("nGenes", "observed.gene.count", "observed gene count"),
    Pathway.Genes = c("Pathway.Genes", "Pathway Genes", "Background", "background.gene.count", "background gene count"),
    Genes = c("Genes", "matching.proteins.in.your.network..labels.", "matching proteins in your network (labels)")
  )
  
  for (std_name in names(col_mapping)) {
    for (orig_name in col_mapping[[std_name]]) {
      if (orig_name %in% names(df)) {
        names(df)[names(df) == orig_name] <- std_name
        break
      }
    }
  }
  
  required <- c("Pathway", "Enrichment.FDR", "nGenes", "Pathway.Genes", "Genes")
  missing <- required[!required %in% names(df)]
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse=", "))
  }
  
  df <- df %>% filter(!is.na(Genes), Genes != "", trimws(Genes) != "")
  df$Enrichment.FDR <- as.numeric(df$Enrichment.FDR)
  df$nGenes <- as.numeric(df$nGenes)
  df$Pathway.Genes <- as.numeric(df$Pathway.Genes)
  
  df_top <- df %>% arrange(Enrichment.FDR) %>% slice_head(n = top_n)
  if (nrow(df_top) == 0) return(NULL)
  
  df_long <- df_top %>%
    separate_rows(Genes, sep = genes_sep) %>%
    rename(Gene = Genes) %>%
    mutate(Gene = trimws(Gene)) %>%
    filter(Gene != "" & !is.na(Gene))
  
  if (nrow(df_long) == 0) return(NULL)
  
  chord_data <- df_long %>%
    select(Pathway, Gene) %>%
    group_by(Pathway, Gene) %>%
    summarise(Count = n(), .groups = "drop")
  
  gene_ratio_df <- df_top %>%
    mutate(GeneRatio = nGenes / Pathway.Genes) %>%
    select(Pathway, GeneRatio)
  
  pathways <- unique(chord_data$Pathway)
  n_path <- length(pathways)
  if (n_path <= 8) {
    pathway_colors <- setNames(brewer.pal(n_path, "Set2"), pathways)
  } else {
    pathway_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(n_path), pathways)
  }
  chord_data$col <- pathway_colors[chord_data$Pathway]
  
  all_sectors <- c(pathways, unique(chord_data$Gene))
  max_ratio <- max(gene_ratio_df$GeneRatio, na.rm = TRUE)
  
  output_png <- paste0(output_name, "_Circular_GeneRatio.png")
  png(filename = output_png, width = width, height = height, units = "in", res = 300)
  
  # Set global font to Times (serif)
  par(family = "serif")
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2.5, clock.wise = TRUE)
  
  chordDiagram(chord_data,
               order = all_sectors,
               transparency = 0.3,
               col = chord_data$col,
               link.border = NA,
               link.lwd = 0.3,
               directional = 0,
               annotationTrack = c("grid", "axis"),
               annotationTrackHeight = c(0.04, 0.02),
               preAllocateTracks = list(track.height = 0.15))
  
  # Track 2: labels with larger serif font
  circos.track(track.index = 2, panel.fun = function(x, y) {
    label <- CELL_META$sector.index
    if (nchar(label) > 25) label <- paste0(substr(label, 1, 22), "...")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1,
                label, 
                facing = "clockwise", 
                niceFacing = TRUE,
                adj = c(0, 0.5), 
                cex = 0.7,
                font = 2,
                family = "serif")
  }, bg.border = NA)
  
  # Track 3: Gene Ratio bars
  circos.track(track.index = 3, 
               factors = all_sectors,
               ylim = c(0, max_ratio * 1.05),
               panel.fun = function(x, y) {
                 sector_name <- CELL_META$sector.index
                 if (sector_name %in% gene_ratio_df$Pathway) {
                   ratio <- gene_ratio_df$GeneRatio[gene_ratio_df$Pathway == sector_name]
                   circos.rect(xleft = CELL_META$xlim[1], ybottom = 0,
                               xright = CELL_META$xlim[2], ytop = ratio,
                               col = pathway_colors[sector_name], border = NA)
                   circos.text(CELL_META$xcenter, ratio + 0.02 * max_ratio,
                               labels = round(ratio, 2), cex = 0.5, adj = c(0.5, 0),
                               family = "serif")
                 } else {
                   circos.rect(xleft = CELL_META$xlim[1], ybottom = 0,
                               xright = CELL_META$xlim[2], ytop = 0.02 * max_ratio,
                               col = "grey80", border = NA)
                 }
               })
  
  title(title_text, cex.main = 1.4, font.main = 2, line = -1, family = "serif")
  
  # Legends
  legend("bottomleft", legend = pathways, fill = pathway_colors[pathways],
         title = "Pathways", cex = 0.65, bty = "o", box.lwd = 0.5, bg = "white",
         text.font = 2, title.font = 2)
  
  legend("bottomright",
         legend = c(paste("Gene Ratio =", round(min(gene_ratio_df$GeneRatio),2)),
                    paste("           →", round(max_ratio,2))),
         fill = c(pathway_colors[pathways[1]], pathway_colors[pathways[1]]),
         title = "Bar height = Gene Ratio",
         cex = 0.65, bty = "o", box.lwd = 0.5, bg = "white", border = NA,
         text.font = 1, title.font = 2)
  
  dev.off()
  message("✅ Saved: ", output_png)
}

setwd("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis")

create_circular_plot("kegg_enrichment.txt", "KEGG", "KEGG Enrichment", top_n = 10)
create_circular_plot("GO biological process.txt", "GO_BP", "GO Biological Process", top_n = 10)
create_circular_plot("Reactome pathway.txt", "Reactome", "Reactome Enrichment", top_n = 10, genes_sep = ",")

#============================All Picture====================================

library(cowplot)
library(magick)
library(ggplot2)

# Read images
img1 <- image_read("KEGG_Circular_GeneRatio.png")
img2 <- image_read("GO_BP_Circular_GeneRatio.png")
img3 <- image_read("Reactome_Circular_GeneRatio.png")

# Convert to grob for cowplot
g1 <- ggdraw() + draw_image(img1)
g2 <- ggdraw() + draw_image(img2)
g3 <- ggdraw() + draw_image(img3)

# Arrange in one row (three plots) with small labels a, b, c
combined <- plot_grid(g1, g2, g3, 
                      ncol = 3, 
                      labels = c("a", "b", "c"),
                      label_size = 20,
                      label_fontface = "bold",
                      label_colour = "white",
                      hjust = -0.5, vjust = 1.5)

# Save with 300 dpi quality
ggsave("Three_Circular_Plots.png", combined, 
       width = 20, height = 8, dpi = 300)

# Display
print(combined)

#=====================================Test For Final=========================

create_circular_plot <- function(file_path, output_name, title_text, 
                                 top_n = 10, 
                                 width = 9, height = 9,
                                 genes_sep = "\\s+") {
  
  # (same column mapping and data preparation as before)
  if (grepl("\\.csv$", file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  col_mapping <- list(
    Pathway = c("Pathway", "GO.biological.process", "term.description", "term description"),
    Enrichment.FDR = c("Enrichment.FDR", "Enrichment FDR", "FDR", "false.discovery.rate", "false discovery rate"),
    nGenes = c("nGenes", "observed.gene.count", "observed gene count"),
    Pathway.Genes = c("Pathway.Genes", "Pathway Genes", "Background", "background.gene.count", "background gene count"),
    Genes = c("Genes", "matching.proteins.in.your.network..labels.", "matching proteins in your network (labels)")
  )
  
  for (std_name in names(col_mapping)) {
    for (orig_name in col_mapping[[std_name]]) {
      if (orig_name %in% names(df)) {
        names(df)[names(df) == orig_name] <- std_name
        break
      }
    }
  }
  
  required <- c("Pathway", "Enrichment.FDR", "nGenes", "Pathway.Genes", "Genes")
  missing <- required[!required %in% names(df)]
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse=", "))
  
  df <- df %>% filter(!is.na(Genes), Genes != "", trimws(Genes) != "")
  df$Enrichment.FDR <- as.numeric(df$Enrichment.FDR)
  df$nGenes <- as.numeric(df$nGenes)
  df$Pathway.Genes <- as.numeric(df$Pathway.Genes)
  
  df_top <- df %>% arrange(Enrichment.FDR) %>% slice_head(n = top_n)
  if (nrow(df_top) == 0) return(NULL)
  
  df_long <- df_top %>%
    separate_rows(Genes, sep = genes_sep) %>%
    rename(Gene = Genes) %>%
    mutate(Gene = trimws(Gene)) %>%
    filter(Gene != "" & !is.na(Gene))
  if (nrow(df_long) == 0) return(NULL)
  
  chord_data <- df_long %>%
    select(Pathway, Gene) %>%
    group_by(Pathway, Gene) %>%
    summarise(Count = n(), .groups = "drop")
  
  gene_ratio_df <- df_top %>%
    mutate(GeneRatio = nGenes / Pathway.Genes) %>%
    select(Pathway, GeneRatio)
  
  pathways <- unique(chord_data$Pathway)
  n_path <- length(pathways)
  if (n_path <= 8) {
    pathway_colors <- setNames(brewer.pal(n_path, "Set2"), pathways)
  } else {
    pathway_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(n_path), pathways)
  }
  chord_data$col <- pathway_colors[chord_data$Pathway]
  
  all_sectors <- c(pathways, unique(chord_data$Gene))
  max_ratio <- max(gene_ratio_df$GeneRatio, na.rm = TRUE)
  
  output_png <- paste0(output_name, "_Circular_GeneRatio.png")
  png(filename = output_png, width = width, height = height, units = "in", res = 300)
  
  par(family = "serif", mar = c(0.5, 0.5, 0.5, 0.5))
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2.5, clock.wise = TRUE,
             cell.padding = c(0.02, 0.02, 0.02, 0.02))
  
  chordDiagram(chord_data,
               order = all_sectors,
               transparency = 0.3,
               col = chord_data$col,
               link.border = NA,
               link.lwd = 0.3,
               directional = 0,
               annotationTrack = c("grid", "axis"),
               annotationTrackHeight = c(0.04, 0.02),
               preAllocateTracks = list(track.height = 0.18))
  
  circos.track(track.index = 2, panel.fun = function(x, y) {
    label <- CELL_META$sector.index
    if (nchar(label) > 25) label <- paste0(substr(label, 1, 22), "...")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.12,
                label, 
                facing = "clockwise", 
                niceFacing = TRUE,
                adj = c(0, 0.5), 
                cex = 0.85,
                font = 2,
                family = "serif")
  }, bg.border = NA)
  
  circos.track(track.index = 3, 
               factors = all_sectors,
               ylim = c(0, max_ratio * 1.05),
               panel.fun = function(x, y) {
                 sector_name <- CELL_META$sector.index
                 if (sector_name %in% gene_ratio_df$Pathway) {
                   ratio <- gene_ratio_df$GeneRatio[gene_ratio_df$Pathway == sector_name]
                   circos.rect(xleft = CELL_META$xlim[1], ybottom = 0,
                               xright = CELL_META$xlim[2], ytop = ratio,
                               col = pathway_colors[sector_name], border = NA)
                   circos.text(CELL_META$xcenter, ratio + 0.02 * max_ratio,
                               labels = round(ratio, 2), cex = 0.75, adj = c(0.5, 0),
                               font = 2, family = "serif")
                 } else {
                   circos.rect(xleft = CELL_META$xlim[1], ybottom = 0,
                               xright = CELL_META$xlim[2], ytop = 0.02 * max_ratio,
                               col = "grey80", border = NA)
                 }
               })
  
  title(title_text, cex.main = 1.3, font.main = 2, line = -1, family = "serif")
  
  legend("bottomleft", legend = pathways, fill = pathway_colors[pathways],
         title = "Pathways", cex = 0.8, bty = "o", box.lwd = 0.5, bg = "white",
         text.font = 2, title.font = 2)
  
  legend("bottomright",
         legend = c(paste("Gene Ratio =", round(min(gene_ratio_df$GeneRatio),3)),
                    paste("           →", round(max_ratio,3))),
         fill = c(pathway_colors[pathways[1]], pathway_colors[pathways[1]]),
         title = "Bar height = Gene Ratio",
         cex = 0.8, bty = "o", box.lwd = 0.5, bg = "white", border = NA,
         text.font = 1, title.font = 2)
  
  dev.off()
  message("✅ Saved: ", output_png)
}

setwd("G:/My Paper/Full Paper/3 Breed haplo/Data/ShinyGO software for functional analysis")

create_circular_plot("kegg_enrichment.txt", "KEGG", "KEGG Enrichment", top_n = 10)
create_circular_plot("GO biological process.txt", "GO_BP", "GO Biological Process", top_n = 10)
create_circular_plot("Reactome pathway.txt", "Reactome", "Reactome Enrichment", top
