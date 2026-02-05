
## ---- Project paths -----------------------------------------------------
project_root <- normalizePath(getwd())

data_root    <- file.path(project_root, "data")
results_root <- file.path(project_root, "results")
annot_root   <- file.path(project_root, "annotation")

dir.create(results_root, showWarnings = FALSE, recursive = TRUE)

gff_path <- file.path(annot_root, "hsa.gff3.txt")

if (!file.exists(gff_path)) {
  stop(
    "❌ GFF annotation not found at: ", gff_path,
    "\nPut your file here: ", annot_root,
    "\nExpected filename: hsa.gff3.txt"
  )
}

## ---- Libraries ---------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tibble)

## ---- Metadata ----------------------------------------------------------
prepare_sample_metadata <- function(sample_data, expression_columns) {
  # Rename Sample.ID to sample_id if needed
  if ("Sample.ID" %in% colnames(sample_data)) {
    sample_data <- dplyr::rename(sample_data, sample_id = Sample.ID)
  }
  
  if (!"sample_id" %in% colnames(sample_data)) {
    stop("❌ No 'sample_id' column found in sample_data.")
  }
  
  # Keep only samples present in expression matrix
  sample_data <- sample_data %>%
    filter(sample_id %in% expression_columns)
  
  # Assign tumor/normal labels from TCGA barcode
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
  
  paired_sample_data <- sample_data %>%
    filter(patient_id %in% paired_patients) %>%
    arrange(patient_id, desc(condition))
  
  paired_sample_count <- length(paired_patients)
  
  return(list(
    metadata = paired_sample_data,
    paired_count = paired_sample_count
  ))
}

## ---- Annotation --------------------------------------------------------
load_gff_annotation <- function(gff_path) {
  gff <- read.delim(
    gff_path, header = FALSE, comment.char = "#", sep = "\t",
    stringsAsFactors = FALSE
  )
  colnames(gff) <- c("seqid", "source", "type", "start", "end", "score",
                     "strand", "phase", "attributes")
  
  mature <- gff[gff$type == "miRNA", ]
  
  mature$MIMAT_ID   <- sub(".*ID=([^;]+);.*", "\\1", mature$attributes)
  mature$miRNA_Name <- sub(".*Name=([^;]+);.*", "\\1", mature$attributes)
  mature$Precursor  <- sub(".*Derives_from=([^;]+).*", "\\1", mature$attributes)
  
  mirna_lookup <- unique(mature[, c("MIMAT_ID", "miRNA_Name", "Precursor", "seqid", "start", "end")])
  return(mirna_lookup)
}

## ---- DESeq2 ------------------------------------------------------------
run_deseq <- function(count_data, metadata) {
  metadata <- metadata %>% filter(!is.na(condition))
  
  # Filter low-count miRNAs (sum across kept samples)
  count_data_filtered <- count_data[rowSums(count_data[, metadata$sample_id, drop = FALSE]) >= 10, , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_data_filtered[, metadata$sample_id, drop = FALSE],
    colData = metadata,
    design = ~ patient_id + condition
  )
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  return(as.data.frame(res))
}

plot_deseq_volcano <- function(res_df, output_path, padj_thresh = 0.05, lfc_thresh = 1) {
  res_df <- res_df %>%
    rownames_to_column("Gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      significance = ifelse(padj < padj_thresh & abs(log2FoldChange) > lfc_thresh,
                            "Significant", "Not Significant"),
      neg_log10_padj = -log10(padj)
    )
  
  top_genes <- res_df %>%
    filter(significance == "Significant") %>%
    arrange(padj) %>%
    slice_head(n = 15)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = "DESeq2 Differential Expression", x = "log2 Fold Change", y = "-log10(padj)") +
    theme_minimal()
  
  ggsave(output_path, plot = p, width = 7, height = 6, dpi = 600)
  return(top_genes$Gene)
}

## ---- Arm switching -----------------------------------------------------
calculate_arm_switching <- function(annotated_counts, mirna_lookup) {
  annotated_counts$MIMAT_ID <- sub("\\..*", "", annotated_counts$MIMAT_ID)
  mirna_lookup$MIMAT_ID <- sub("\\..*", "", mirna_lookup$MIMAT_ID)
  mirna_lookup$MIMAT_ID <- gsub("_\\d+$", "", mirna_lookup$MIMAT_ID)
  
  annotated_counts <- annotated_counts %>%
    left_join(mirna_lookup, by = "MIMAT_ID")
  
  if (!"miRNA_Name" %in% colnames(annotated_counts) || all(is.na(annotated_counts$miRNA_Name))) {
    stop("❌ ERROR: 'miRNA_Name' is missing or all NA after join. Check MIMAT_ID format.")
  }
  
  annotated_counts <- annotated_counts %>%
    mutate(strand_arm = case_when(
      grepl("-5p$", miRNA_Name) ~ "5p",
      grepl("-3p$", miRNA_Name) ~ "3p",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(strand_arm))
  
  arm_summary <- annotated_counts %>%
    group_by(Precursor, condition, strand_arm) %>%
    summarise(total_reads = sum(count), .groups = "drop")
  
  arm_wide <- arm_summary %>%
    pivot_wider(names_from = c(condition, strand_arm), values_from = total_reads, values_fill = 0) %>%
    mutate(
      LOR = log((Tumor_5p + 1) / (Tumor_3p + 1)) - log((Normal_5p + 1) / (Normal_3p + 1)),
      SD = sqrt(1 / (Tumor_5p + 1) + 1 / (Tumor_3p + 1) + 1 / (Normal_5p + 1) + 1 / (Normal_3p + 1)),
      Z_score = LOR / SD,
      p_value = 2 * pnorm(-abs(Z_score)),
      qvals = p.adjust(p_value, method = "fdr")
    )
  
  return(arm_wide)
}

plot_volcano <- function(arm_data, output_path, q_threshold = 0.01, lor_threshold = 1) {
  arm_data <- arm_data %>%
    mutate(
      significance = ifelse(qvals < q_threshold & abs(LOR) > lor_threshold, "Significant", "Not Significant"),
      neg_log10_q = -log10(qvals)
    )
  
  top_mirnas <- arm_data %>%
    filter(qvals < q_threshold & abs(LOR) > lor_threshold) %>%
    arrange(qvals) %>%
    slice_head(n = 15)
  
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

## ---- Heatmap -----------------------------------------------------------
generate_dual_heatmap <- function(raw_counts_annotated, sample_data, sig_mirnas = NULL, output_pdf = NULL) {
  library(ComplexHeatmap)
  library(circlize)
  
  sample_data_unique <- sample_data %>%
    distinct(sample_id, .keep_all = TRUE)
  
  arm_counts <- raw_counts_annotated %>%
    group_by(sample_id, Precursor, strand_arm) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = strand_arm, values_from = total, values_fill = 0) %>%
    mutate(log2_5p_3p = log2((`5p` + 1) / (`3p` + 1)))
  
  annotated <- arm_counts %>%
    left_join(sample_data_unique[, c("sample_id", "patient_id", "condition")], by = "sample_id")
  
  log_data_matrix <- annotated %>%
    filter(!is.na(Precursor), !is.na(log2_5p_3p)) %>%
    mutate(column_label = paste0(tolower(condition), "_", patient_id)) %>%
    select(Precursor, column_label, log2_5p_3p) %>%
    group_by(Precursor, column_label) %>%
    summarise(value = mean(log2_5p_3p), .groups = "drop") %>%
    pivot_wider(names_from = column_label, values_from = value) %>%
    column_to_rownames("Precursor") %>%
    as.matrix()
  
  matched_ids <- sample_data_unique %>%
    group_by(patient_id) %>%
    filter(n_distinct(condition) == 2) %>%
    distinct(patient_id) %>%
    pull(patient_id)
  
  ordered_cols <- as.vector(t(outer(c("normal", "tumor"), matched_ids, paste, sep = "_")))
  ordered_cols <- intersect(ordered_cols, colnames(log_data_matrix))
  matrix_matched <- log_data_matrix[, ordered_cols, drop = FALSE]
  
  matrix_matched <- matrix_matched[
    apply(matrix_matched, 1, function(x) sd(x, na.rm = TRUE) > 0.5),
  ]
  
  if (!is.null(sig_mirnas)) {
    matrix_matched <- matrix_matched[rownames(matrix_matched) %in% sig_mirnas, , drop = FALSE]
  }
  
  tumor_cols  <- grep("^tumor_", colnames(matrix_matched), value = TRUE)
  normal_cols <- grep("^normal_", colnames(matrix_matched), value = TRUE)
  
  tumor_matrix  <- matrix_matched[, tumor_cols, drop = FALSE]
  normal_matrix <- matrix_matched[, normal_cols, drop = FALSE]
  
  colnames(tumor_matrix)  <- gsub("^tumor_", "", colnames(tumor_matrix))
  colnames(normal_matrix) <- gsub("^normal_", "", colnames(normal_matrix))
  
  max_val <- max(abs(matrix_matched), na.rm = TRUE)
  col_fun <- colorRamp2(c(-max_val, 0, max_val), c("skyblue", "black", "salmon"))
  
  normal_heatmap <- Heatmap(
    normal_matrix,
    name = "Log2Ratio",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_title = "Non-Tumor",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 7),
    heatmap_legend_param = list(title = "Log2(5p/3p)")
  )
  
  tumor_heatmap <- Heatmap(
    tumor_matrix,
    name = "Log2Ratio",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_title = "Tumor",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 7),
    heatmap_legend_param = list(title = "Log2(5p/3p)")
  )
  
  if (!is.null(output_pdf)) {
    pdf(output_pdf, width = 10, height = 6)
    draw(normal_heatmap + tumor_heatmap, merge_legend = TRUE,
         column_title_gp = gpar(fontsize = 10, fontface = "bold"))
    dev.off()
    message("✅ Heatmap saved to: ", output_pdf)
  } else {
    draw(normal_heatmap + tumor_heatmap, merge_legend = TRUE,
         column_title_gp = gpar(fontsize = 10, fontface = "bold"))
    message("ℹ️ Heatmap rendered in R session.")
  }
}

## ---- Pipeline ----------------------------------------------------------
run_arm_bias_pipeline <- function(cancer_type) {
  message("Running analysis for: ", cancer_type)
  
  data_dir    <- file.path(data_root, cancer_type)
  results_dir <- file.path(results_root, cancer_type)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(data_dir)) {
    message("⚠️ Skipping ", cancer_type, ": data folder not found: ", data_dir)
    return(invisible(NULL))
  }
  
  expr_file   <- list.files(data_dir, pattern = "miRs_counts", full.names = TRUE)
  sample_file <- list.files(data_dir, pattern = "miRs_samples", full.names = TRUE)
  
  if (length(expr_file) == 0 || length(sample_file) == 0) {
    message("⚠️ Skipping ", cancer_type, ": Expression or sample file not found.")
    return(invisible(NULL))
  }
  
  raw_counts  <- read.delim(expr_file[1], row.names = 1, check.names = FALSE)
  sample_data <- read.delim(sample_file[1], header = TRUE, stringsAsFactors = FALSE)
  
  mirna_lookup <- load_gff_annotation(gff_path)
  
  meta_result <- prepare_sample_metadata(sample_data, colnames(raw_counts))
  metadata    <- meta_result$metadata
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
  write.csv(data.frame(Gene = top_deseq_genes),
            file.path(results_dir, "top_deseq_genes.csv"),
            row.names = FALSE)
  
  long_counts <- raw_counts %>%
    rownames_to_column("MIMAT_ID") %>%
    pivot_longer(-MIMAT_ID, names_to = "sample_id", values_to = "count") %>%
    left_join(metadata, by = "sample_id")
  
  write.csv(long_counts, file.path(results_dir, "long_counts_annotated.csv"), row.names = FALSE)
  
  arm_wide <- calculate_arm_switching(long_counts, mirna_lookup)
  write.csv(arm_wide, file.path(results_dir, "arm_switching_stats.csv"), row.names = FALSE)
  
  top_precursors <- plot_volcano(arm_wide, file.path(results_dir, "volcano_plot.tiff"))
  write.csv(data.frame(Precursor = top_precursors),
            file.path(results_dir, "top_volcano_miRNAs.csv"),
            row.names = FALSE)
  
  if (length(top_precursors) > 0) {
    message("✅ Generating side-by-side tumor/normal heatmap for top miRNAs...")
    
    annotated_long_counts <- long_counts %>%
      left_join(mirna_lookup, by = "MIMAT_ID") %>%
      mutate(strand_arm = case_when(
        grepl("-5p$", miRNA_Name) ~ "5p",
        grepl("-3p$", miRNA_Name) ~ "3p",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(strand_arm))
    
    write.csv(annotated_long_counts,
              file.path(results_dir, "long_counts_with_strand.csv"),
              row.names = FALSE)
    
    generate_dual_heatmap(
      raw_counts_annotated = annotated_long_counts,
      sample_data = metadata,
      sig_mirnas = top_precursors,
      output_pdf = file.path(results_dir, paste0("log2ratio_", cancer_type, "_heatmap.pdf"))
    )
  } else {
    message("ℹ️ No significant miRNAs found for heatmap.")
  }
  
  invisible(TRUE)
}

run_all_cancers <- function(base_data_dir = data_root) {
  if (!dir.exists(base_data_dir)) {
    stop("❌ data folder not found: ", base_data_dir)
  }
  
  cancer_types <- list.dirs(base_data_dir, full.names = FALSE, recursive = FALSE)
  
  for (cancer in cancer_types) {
    tryCatch(
      run_arm_bias_pipeline(cancer),
      error = function(e) message("❌ Error in ", cancer, ": ", e$message)
    )
  }
  
  invisible(TRUE)
}