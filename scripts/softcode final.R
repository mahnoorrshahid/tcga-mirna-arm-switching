## ---- Project paths -----------------------------------------------------
project_root <- normalizePath(getwd())

data_root    <- file.path(project_root, "data")
results_root <- file.path(project_root, "results")
annot_root   <- file.path(project_root, "annotation")

dir.create(results_root, showWarnings = FALSE, recursive = TRUE)

gff_path <- file.path(annot_root, "hsa.gff3.txt")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tibble)

prepare_sample_metadata <- function(sample_data, expression_columns) {
  # Rename Sample.ID to sample_id if needed (dot notation, not space)
  if ("Sample.ID" %in% colnames(sample_data)) {
    sample_data <- dplyr::rename(sample_data, sample_id = Sample.ID)
  }
  
  # If sample_id column already exists, just proceed
  if (!"sample_id" %in% colnames(sample_data)) {
    stop("❌ No 'sample_id' column found in sample_data.")
  }
  
  # Filter to samples that match expression matrix
  sample_data <- sample_data %>%
    filter(sample_id %in% expression_columns)
  
  # Assign tumor/normal labels
  sample_data <- sample_data %>%
    mutate(
      condition = ifelse(substr(sample_id, 14, 15) == "01", "Tumor",
                         ifelse(substr(sample_id, 14, 15) == "11", "Normal", NA)),
      patient_id = substr(sample_id, 1, 12)
    )
  
  # Identify patients with both Tumor + Normal
  paired_patients <- sample_data %>%
    filter(!is.na(condition)) %>%
    group_by(patient_id) %>%
    summarise(n = n_distinct(condition), .groups = "drop") %>%
    filter(n == 2) %>%
    pull(patient_id)
  
  # Keep only paired
  paired_sample_data <- sample_data %>%
    filter(patient_id %in% paired_patients) %>%
    arrange(patient_id, desc(condition))
  
  paired_sample_count <- length(paired_patients)
  
  return(list(
    metadata = paired_sample_data,
    paired_count = paired_sample_count
  ))
}


# Load GFF3 and extract miRNA annotation

load_gff_annotation <- function(gff_path) {
  gff <- read.delim(gff_path, header = FALSE, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
  colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  # Filter mature miRNAs
  
  mature <- gff[gff$type == "miRNA", ]
  
  # Extract IDs and names
  
  mature$MIMAT_ID <- sub(".*ID=([^;]+);.*", "\\1", mature$attributes)
  mature$miRNA_Name <- sub(".*Name=([^;]+);.*", "\\1", mature$attributes)
  mature$Precursor <- sub(".*Derives_from=([^;]+).*", "\\1", mature$attributes)
  
  # Extract relevant fields for annotation
  
  mirna_lookup <- unique(mature[, c("MIMAT_ID", "miRNA_Name", "Precursor", "seqid", "start", "end")])
  return(mirna_lookup)
}

# Differential expression 
run_deseq <- function(count_data, metadata) {
  metadata <- metadata %>% filter(!is.na(condition))
  
  # Remove miRNAs with low counts across all samples
  count_data_filtered <- count_data[rowSums(count_data[, metadata$sample_id]) >= 10, ]
  
  
  # Create DESeq2 dataset using patient pairing in the design
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_data_filtered[, metadata$sample_id],
    colData = metadata,
    design = ~ patient_id + condition
  )
  dds <- DESeq(dds)
  res <- results(dds)
  return(as.data.frame(res))
}

# Generate DESeq2 volcano plot

plot_deseq_volcano <- function(res_df, output_path, padj_thresh = 0.05, lfc_thresh = 1) {
  res_df <- res_df %>%
    rownames_to_column("Gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      significance = ifelse(padj < padj_thresh & abs(log2FoldChange) > lfc_thresh, "Significant", "Not Significant"),
      neg_log10_padj = -log10(padj)
    )
  
  # Get top 15 most significant genes
  
  top_genes <- res_df %>%
    filter(significance == "Significant") %>%
    arrange(padj) %>%
    slice_head(n = 15)
  
  # Plot
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = "DESeq2 Differential Expression", x = "log2 Fold Change", y = "-log10(padj)") +
    theme_minimal()
  
  # Save to file
  
  ggsave(output_path, plot = p, width = 7, height = 6, dpi = 600)
  return(top_genes$Gene)
}

# arm switching calculation 

calculate_arm_switching <- function(annotated_counts, mirna_lookup) {
  
  # Clean up MIMAT_ID formats to match
  
  annotated_counts$MIMAT_ID <- sub("\\..*", "", annotated_counts$MIMAT_ID)
  mirna_lookup$MIMAT_ID <- sub("\\..*", "", mirna_lookup$MIMAT_ID)
  mirna_lookup$MIMAT_ID <- gsub("_\\d+$", "", mirna_lookup$MIMAT_ID)
  
  # Join annotations
  
  annotated_counts <- annotated_counts %>%
    left_join(mirna_lookup, by = "MIMAT_ID")
  
  # Sanity check
  
  if (!"miRNA_Name" %in% colnames(annotated_counts) || all(is.na(annotated_counts$miRNA_Name))) {
    stop("\u274c ERROR: 'miRNA_Name' is missing or all NA after join. Check MIMAT_ID format.")
  }
  
  # Annotate strand arm (5p/3p)
  
  annotated_counts <- annotated_counts %>%
    mutate(strand_arm = case_when(
      grepl("-5p$", miRNA_Name) ~ "5p",
      grepl("-3p$", miRNA_Name) ~ "3p",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(strand_arm))
  
  # Summarise counts by precursor and condition/strand
  
  arm_summary <- annotated_counts %>%
    group_by(Precursor, condition, strand_arm) %>%
    summarise(total_reads = sum(count), .groups = "drop")
  
  # Reshape and calculate log-odds ratio, standard deviation, z-score, and FDR
  
  arm_wide <- arm_summary %>%
    pivot_wider(names_from = c(condition, strand_arm), values_from = total_reads, values_fill = 0)
  
  arm_wide <- arm_wide %>%
    mutate(
      LOR = log((Tumor_5p + 1) / (Tumor_3p + 1)) - log((Normal_5p + 1) / (Normal_3p + 1)),
      SD = sqrt(1 / (Tumor_5p + 1) + 1 / (Tumor_3p + 1) + 1 / (Normal_5p + 1) + 1 / (Normal_3p + 1)),
      Z_score = LOR / SD,
      p_value = 2 * pnorm(-abs(Z_score)),
      qvals = p.adjust(p_value, method = "fdr")
    )
  
  return(arm_wide)
}

# Volcano plot for arm switching results

plot_volcano <- function(arm_data, output_path, q_threshold = 0.01, lor_threshold = 1) {
  arm_data <- arm_data %>%
    mutate(
      significance = ifelse(qvals < q_threshold & abs(LOR) > lor_threshold, "Significant", "Not Significant"),
      neg_log10_q = -log10(qvals)
    )
  
  # Top hits for annotation
  
  top_mirnas <- arm_data %>%
    filter(qvals < q_threshold & abs(LOR) > lor_threshold) %>%
    arrange(qvals) %>%
    slice_head(n = 15)
  
  # Plot
  
  p <- ggplot(arm_data, aes(x = LOR, y = neg_log10_q, color = significance)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_vline(xintercept = c(-lor_threshold, lor_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(q_threshold), linetype = "dashed") +
    geom_text_repel(data = top_mirnas, aes(label = Precursor), size = 3.0) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = "miRNA Arm Usage Bias", x = "Log Odds Ratio (LOR)", y = "-log10(q-value)") +
    theme_minimal()
  
  ggsave(output_path, plot = p, width = 7, height = 6, dpi = 600)
  
  return(top_mirnas$Precursor)
  
}

## Heatmap for log2(5p/3p) ratios across matched tumor and normal samples
# Draws side by side heatmaps for tumor and normal tissues
# Filters top miRNAs and restricts to significant hits
# Saves output as pdf


generate_dual_heatmap <- function(raw_counts_annotated, sample_data, sig_mirnas = NULL, output_pdf = NULL) {
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  
  # Ensure sample_data has unique sample_id entries
  sample_data_unique <- sample_data %>%
    distinct(sample_id, .keep_all = TRUE)
  
  # Calculate log2(5p/3p) per sample per precursor
  arm_counts <- raw_counts_annotated %>%
    group_by(sample_id, Precursor, strand_arm) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = strand_arm, values_from = total, values_fill = 0) %>%
    mutate(log2_5p_3p = log2((`5p` + 1) / (`3p` + 1)))
  
  # Annotate with metadata
  annotated <- arm_counts %>%
    left_join(sample_data_unique[, c("sample_id", "patient_id", "condition")], by = "sample_id")
  
  # Create a wide matrix of log2 ratios
  log_data_matrix <- annotated %>%
    filter(!is.na(Precursor), !is.na(log2_5p_3p)) %>%
    mutate(column_label = paste0(tolower(condition), "_", patient_id)) %>%
    select(Precursor, column_label, log2_5p_3p) %>%
    group_by(Precursor, column_label) %>%
    summarise(value = mean(log2_5p_3p), .groups = "drop") %>%
    pivot_wider(names_from = column_label, values_from = value) %>%
    column_to_rownames("Precursor") %>%
    as.matrix()
  
  # Get matched patients
  matched_ids <- sample_data_unique %>%
    group_by(patient_id) %>%
    filter(n_distinct(condition) == 2) %>%
    distinct(patient_id) %>%
    pull(patient_id)
  
  # Order columns: normal, tumor per patient
  ordered_cols <- as.vector(t(outer(c("normal", "tumor"), matched_ids, paste, sep = "_")))
  ordered_cols <- intersect(ordered_cols, colnames(log_data_matrix))
  matrix_matched <- log_data_matrix[, ordered_cols, drop = FALSE]
  
  # Filter for variance
  matrix_matched <- matrix_matched[
    apply(matrix_matched, 1, function(x) sd(x, na.rm = TRUE) > 0.5),
  ]
  
  # Filter to significant miRNAs (if provided)
  if (!is.null(sig_mirnas)) {
    matrix_matched <- matrix_matched[rownames(matrix_matched) %in% sig_mirnas, , drop = FALSE]
  }
  
  # Split tumor and normal
  tumor_cols <- grep("^tumor_", colnames(matrix_matched), value = TRUE)
  normal_cols <- grep("^normal_", colnames(matrix_matched), value = TRUE)
  tumor_matrix <- matrix_matched[, tumor_cols, drop = FALSE]
  normal_matrix <- matrix_matched[, normal_cols, drop = FALSE]
  
  # Clean column names
  colnames(tumor_matrix) <- gsub("^tumor_", "", colnames(tumor_matrix))
  colnames(normal_matrix) <- gsub("^normal_", "", colnames(normal_matrix))
  
  # Define color scale
  max_val <- max(abs(matrix_matched), na.rm = TRUE)
  col_fun <- colorRamp2(c(-max_val, 0, max_val), c("skyblue", "black", "salmon"))
  
  
  # 11. Draw heatmaps
  normal_heatmap <- Heatmap(normal_matrix,
                            name = "Log2Ratio",
                            col = col_fun,
                            cluster_rows = TRUE,
                            cluster_columns = FALSE,
                            column_title = "Non-Tumor",
                            row_names_gp = gpar(fontsize = 9),
                            column_names_gp = gpar(fontsize = 7),
                            heatmap_legend_param = list(title = "Log2(5p/3p)"))
  
  tumor_heatmap <- Heatmap(tumor_matrix,
                           name = "Log2Ratio",
                           col = col_fun,
                           cluster_rows = TRUE,
                           cluster_columns = FALSE,
                           column_title = "Tumor",
                           row_names_gp = gpar(fontsize = 9),
                           column_names_gp = gpar(fontsize = 7),
                           heatmap_legend_param = list(title = "Log2(5p/3p)"))
  
  ht_list <- draw(normal_heatmap + tumor_heatmap, merge_legend = TRUE,
                  column_title_gp = gpar(fontsize = 5, fontface = "bold"))
  
  # save as PDF
  if (!is.null(output_pdf)) {
    pdf(output_pdf, width = 10, height = 6)
    draw(ht_list)
    dev.off()
    message("Heatmap saved to: ", output_pdf)
  } else {
    message("Heatmap rendered in R session.")
  }
}

# Modify pipeline to pass expression column names
run_arm_bias_pipeline <- function(cancer_type) {
  message("Running analysis for: ", cancer_type)
  
  data_dir <- file.path("~/miRNA_TCGA_project", "data", cancer_type)
  results_dir <- file.path("~/miRNA_TCGA_project", "results", cancer_type)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Detect files based on known patterns
  expr_file <- list.files(data_dir, pattern = "miRs_counts", full.names = TRUE)
  sample_file <- list.files(data_dir, pattern = "miRs_samples", full.names = TRUE)
  
  # Validate they exist
  if (length(expr_file) == 0 || length(sample_file) == 0) {
    message("⚠️ Skipping ", cancer_type, ": Expression or sample file not found.")
    return(invisible(NULL))
  }
  
  # Load expression and sample data
  raw_counts <- read.delim(expr_file, row.names = 1, check.names = FALSE)
  sample_data <- read.delim(sample_file, header = TRUE, stringsAsFactors = FALSE)
  sample_data <- if (file.exists(sample_file)) {
    read.delim(sample_file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    data.frame(sample_id = colnames(raw_counts), stringsAsFactors = FALSE)
  }
  
  mirna_lookup <- load_gff_annotation("~/miRNA_TCGA_project/hsa.gff3.txt")
  meta_result <- prepare_sample_metadata(sample_data, colnames(raw_counts))
  metadata <- meta_result$metadata
  paired_count <- meta_result$paired_count
  
  if (paired_count == 0 || nrow(metadata) == 0) {
    writeLines("0", file.path(results_dir, "paired_sample_count.txt"))
    message("⛔ Skipping ", cancer_type, ": No paired tumor-normal samples.")
    return(invisible(NULL))
  }
  
  writeLines(as.character(paired_count), file.path(results_dir, "paired_sample_count.txt"))
  write.csv(metadata, file.path(results_dir, "sample_metadata.csv"), row.names = FALSE)
  
  res_df <- run_deseq(raw_counts, metadata)
  write.csv(res_df, file.path(results_dir, "deseq_results.csv"), row.names = TRUE)
  
  top_deseq_genes <- plot_deseq_volcano(
    res_df,
    file.path(results_dir, "deseq_volcano_plot.tiff")
  )
  write.csv(data.frame(Gene = top_deseq_genes), file.path(results_dir, "top_deseq_genes.csv"), row.names = FALSE)
  
  long_counts <- raw_counts %>%
    rownames_to_column("MIMAT_ID") %>%
    pivot_longer(-MIMAT_ID, names_to = "sample_id", values_to = "count") %>%
    left_join(metadata, by = "sample_id")
  write.csv(long_counts, file.path(results_dir, "long_counts_annotated.csv"), row.names = FALSE)
  
  arm_wide <- calculate_arm_switching(long_counts, mirna_lookup)
  write.csv(arm_wide, file.path(results_dir, "arm_switching_stats.csv"), row.names = FALSE)
  
  top_precursors <- plot_volcano(arm_wide, file.path(results_dir, "volcano_plot.tiff"))
  write.csv(data.frame(Precursor = top_precursors), file.path(results_dir, "top_volcano_miRNAs.csv"), row.names = FALSE)
  
  if (length(top_precursors) > 0) {
    message(" Generating side-by-side tumor/normal heatmap for top miRNAs..")
    
    annotated_long_counts <- long_counts %>%
      left_join(mirna_lookup, by = "MIMAT_ID") %>%
      mutate(strand_arm = case_when(
        grepl("-5p$", miRNA_Name) ~ "5p",
        grepl("-3p$", miRNA_Name) ~ "3p",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(strand_arm))
    
    write.csv(annotated_long_counts, file.path(results_dir, "long_counts_with_strand.csv"), row.names = FALSE)
    generate_dual_heatmap(
      raw_counts_annotated = annotated_long_counts,
      sample_data = metadata,
      sig_mirnas = top_precursors,
      output_pdf = file.path(results_dir, paste0("log2ratio_", cancer_type, "_heatmap.pdf"))
    )
  } else {
    message(" No significant miRNAs found for heatmap.")
  }
}

# Batch run over all cancer types from folder
run_all_cancers <- function(base_data_dir = "~/miRNA_TCGA_project/data") {
  cancer_types <- list.dirs(base_data_dir, full.names = FALSE, recursive = FALSE)
  for (cancer in cancer_types) {
    tryCatch({
      run_arm_bias_pipeline(cancer)
    }, error = function(e) {
      message("❌ Error in ", cancer, ": ", e$message)
    })
  }
}

setwd ("~/miRNA_TCGA_project")

run_all_cancers()
