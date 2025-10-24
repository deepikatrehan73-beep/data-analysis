# data-analysis
A script for microarray normalization in R
Combine sample information
sample_ids <- c(blca_sample_ids, normal_sample_ids)
conditions <- c(blca_conditions, normal_conditions)
sample_info <- data.frame(sample_id = sample_ids, condition = conditions)


blca_files <- list.files("C:/Users/Dell/Documents/BLCA", pattern = ".tsv", full.names = TRUE)
> blca_sample_ids <- sub("^.*/BLCA_(\\d+)\\.tsv$", "\\1", blca_files)
> blca_conditions <- rep("BLCA", length(blca_sample_ids))
> normal_files <- list.files("C:/Users/Dell/Documents/Normal Bladder", pattern = ".tsv", full.names = TRUE)
> normal_sample_ids <- sub("^.*/Normal_(\\d+)\\.tsv$", "\\1", normal_files)
> normal_conditions <- rep("Normal", length(normal_sample_ids))

10-06-2023



blca_files <- list.files("E:/BLCA", pattern = ".tsv", full.names = TRUE)
> blca_sample_ids <- sub("^.*/BLCA_(\\d+)\\.tsv$", "\\1", blca_files)
> blca_conditions <- rep("BLCA", length(blca_sample_ids))
> normal_files <- list.files("E:/Normal Bladder", pattern = ".tsv", full.names = TRUE)
> normal_sample_ids <- sub("^.*/Normal_(\\d+)\\.tsv$", "\\1", normal_files)
blca_counts <- lapply(blca_files, function(file) {
+     data <- read.table(file, header = TRUE, sep = "\t")
+     data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
+     colnames(data) <- c("gene_id", "gene_name", "fpkm_unstranded")
+     data <- data[complete.cases(data), ]
+     data <- data.frame(data)
+     data$fpkm_unstranded <- as.numeric(data$fpkm_unstranded)
+     data$condition <- "BLCA"
+ data
+ })
> count_matrix_blca <- do.call(rbind, blca_counts)
> normal_counts <- lapply(normal_files, function(file) {
+     data <- read.table(file, header = TRUE, sep = "\t")
+     data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
+     colnames(data) <- c("gene_id", "gene_name", "fpkm_unstranded")
+     data <- data[complete.cases(data), ]
+     data <- data.frame(data)
+     data$fpkm_unstranded <- as.numeric(data$fpkm_unstranded)
+     data$condition <- "Normal"
+ data
+ })
> count_matrix_normal <- do.call(rbind, normal_counts)
> count_matrix <- rbind(count_matrix_blca, count_matrix_normal)
sample_info <- data.frame(
+     sample_id = colnames(count_matrix),
+     gene_id = count_matrix$gene_id,
+     gene_name = count_matrix$gene_name,
+     fpkm_unstranded = count_matrix$fpkm_unstranded,
+     condition = count_matrix$condition,
+     stringsAsFactors = FALSE
+ )

sample_info <- data.frame(sample_id = character(),
+                           gene_id = character(),
+                           gene_name = character(),
+                           fpkm_unstranded = numeric(),
+                           condition = character()) 
> for (file in blca_files) {
+     data <- read.table(file, header = TRUE, sep = "\t")
+     sample_id <- sub("^.*/BLCA_(\\d+)\\.tsv$", "\\1", file)
+     sample_data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
+     sample_data$condition <- "BLCA"  # Add 'condition' column and assign 'BLCA'
+     sample_data <- data.frame(sample_id = sample_id, sample_data)
+     sample_info <- rbind(sample_info, sample_data)
+ }
> for (file in normal_files) {
+     data <- read.table(file, header = TRUE, sep = "\t")
+     sample_id <- sub("^.*/Normal_(\\d+)\\.tsv$", "\\1", file)
+     sample_data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
+     sample_data$condition <- "Normal" 
+ sample_data <- data.frame(sample_id = sample_id, sample_data)
+ sample_info <- rbind(sample_info, sample_data)
+ }

11/06/2023
blca_files <- list.files("E:/BLCA", pattern = ".tsv", full.names = TRUE)
blca_counts <- lapply(blca_files, function(file) {
+     data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
+     data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
+     colnames(data) <- c("gene_id", "gene_name", "fpkm_unstranded")
+     data <- data[complete.cases(data), ]
+     data$gene_id <- as.character(data$gene_id)
+     data$gene_name <- as.character(data$gene_name)
+     data$fpkm_unstranded <- as.numeric(data$fpkm_unstranded)
+     data$condition <- "BLCA"
+     data
+ })
> count_matrix_blca <- do.call(rbind, blca_counts)
normal_counts <- lapply(normal_files, function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data <- data[, c("gene_id", "gene_name", "fpkm_unstranded")]
  colnames(data) <- c("gene_id", "gene_name", "fpkm_unstranded")
  data <- data[complete.cases(data), ]
  data$gene_id <- as.character(data$gene_id)
  data$gene_name <- as.character(data$gene_name)
  data$fpkm_unstranded <- as.numeric(data$fpkm_unstranded)
  data$condition <- "Normal Bladder"
  data
})
count_matrix_normal <- do.call(rbind, normal_counts)

count_matrix <- rbind(count_matrix_blca, count_matrix_normal)

13/06/2023
DEGs_df <- as.data.frame(DEGs)
> 
> upregulated_genes <- DEGs_df[DEGs_df$log2FoldChange > 0 & DEGs_df$padj < 0.05, ]
> downregulated_genes <- DEGs_df[DEGs_df$log2FoldChange < 0 & DEGs_df$padj < 0.05, ]
> 
> significant_upregulated_genes <- upregulated_genes[upregulated_genes$log2FoldChange > 0 & upregulated_genes$pvalue < 0.05, ]
> significant_downregulated_genes <- downregulated_genes[downregulated_genes$log2FoldChange < 0 & downregulated_genes$pvalue < 0.05, ]
> 
> enhanced_volcano <- EnhancedVolcano(DEGs_df, x = "log2FoldChange", y = "pvalue", lab = rownames(DEGs_df), title = "Volcano Plot",
+                                     xlim = c(-10, 10), ylim = c(0, 15),
+                                     col = ifelse(DEGs_df$log2FoldChange > 0, "red", "blue"))
> 
> enhanced_volcano + 
+     geom_point(data = significant_upregulated_genes, aes(x = log2FoldChange, y = pvalue), color = "red") +
+     geom_point(data = significant_downregulated_genes, aes(x = log2FoldChange, y = pvalue), color = "blue") +
+     geom_point(data = DEGs_df[!(DEGs_df %in% c(significant_upregulated_genes, significant_downregulated_genes)), ],
+                aes(x = log2FoldChange, y = pvalue), color = "black", show.legend = FALSE) +
+     theme_classic()
