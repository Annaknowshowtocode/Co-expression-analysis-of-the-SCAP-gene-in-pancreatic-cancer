args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[1] else "expression.txt"
output_dir <- if (length(args) >= 2) args[2] else "codfinder_output"
gene_of_interest <- if (length(args) >= 3) args[3] else "SCAP"

id_column <- "Entrez_Gene_Id"
symbol_column <- "Hugo_Symbol"
txdb_package <- "TxDb.Hsapiens.UCSC.hg38.knownGene"

window_sizes <- 3:20
expression_quantile <- 0.30
log_transform <- TRUE
zoom_flank <- 5
n_permutations <- 1000
random_seed <- 42

genome_chromosomes <- paste0("chr", c(1:22, "X"))
chromosome_order <- paste0("chr", c(1:22, "X", "Y", "M"))

font_family <- "serif"
colors <- colorRampPalette(c("#2b56a6", "white", "#A0352B"))(101)
color_limit <- 0.6

matrix_dir <- file.path(output_dir, "matrix")
window_dir <- file.path(output_dir, "window_test")
figures_dir <- file.path(output_dir, "figures")
for (d in c(matrix_dir, window_dir, figures_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_packages <- c("GenomicFeatures", txdb_package)
missing <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) BiocManager::install(missing, ask = FALSE, update = FALSE)

library(GenomicFeatures)
library(txdb_package, character.only = TRUE)
txdb <- get(txdb_package)


message("[1/15] Loading expression data: ", input_file)

genes_exp_data <- read.delim(input_file, check.names = FALSE)
stopifnot(all(c(id_column, symbol_column) %in% names(genes_exp_data)))

names(genes_exp_data)[names(genes_exp_data) == id_column] <- "Entrez_Gene_Id"
names(genes_exp_data)[names(genes_exp_data) == symbol_column] <- "Hugo_Symbol"

genes_exp_data <- genes_exp_data[!is.na(genes_exp_data$Entrez_Gene_Id), ]

duplicates <- duplicated(genes_exp_data$Entrez_Gene_Id)
message("Duplicated Entrez IDs (kept): ", sum(duplicates))


message("[2/15] Adding genomic coordinates")

gene_annot <- genes_exp_data[, c("Entrez_Gene_Id", "Hugo_Symbol")]

gene_ranges <- GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)

ids <- as.character(gene_annot$Entrez_Gene_Id)
idx <- match(ids, names(gene_ranges))

gene_annot$chromosome <- as.character(seqnames(gene_ranges))[idx]
gene_annot$start <- start(gene_ranges)[idx]
gene_annot$end <- end(gene_ranges)[idx]
gene_annot$strand <- as.character(strand(gene_ranges))[idx]

keep <- complete.cases(gene_annot[, c("chromosome", "start", "end")])
gene_annot <- gene_annot[keep, , drop = FALSE]
genes_exp_data <- genes_exp_data[keep, , drop = FALSE]


message("[3/15] Sorting genes by genomic position")

chr_number <- match(gene_annot$chromosome, chromosome_order)
ord <- order(chr_number, gene_annot$chromosome, gene_annot$start, gene_annot$end)

gene_annot <- gene_annot[ord, , drop = FALSE]
genes_exp_data <- genes_exp_data[ord, , drop = FALSE]

stopifnot(identical(gene_annot$Entrez_Gene_Id, genes_exp_data$Entrez_Gene_Id))


message("[4/15] Building and filtering the expression matrix")

sample_cols <- setdiff(names(genes_exp_data), c("Entrez_Gene_Id", "Hugo_Symbol"))

stopifnot(all(sapply(genes_exp_data[sample_cols], is.numeric)))

if (length(sample_cols) < 3) stop("At least 3 samples are required for correlation")
exp_matrix <- as.matrix(genes_exp_data[, sample_cols, drop = FALSE])

if (log_transform) exp_matrix <- log2(exp_matrix + 1)

mean_expression <- rowMeans(exp_matrix)

expression_threshold <- quantile(
  mean_expression,
  probs = expression_quantile,
  na.rm = TRUE
)

svg(file.path(figures_dir, "expression_before_filtering.svg"),
    width = 8, height = 5, family = font_family)

hist(
  mean_expression,
  breaks = 60,
  col = "#2b56a6",
  border = "white",
  main = "Gene expression before filtering",
  xlab = "Mean expression across samples",
  ylab = "Number of genes"
)

abline(
  v = expression_threshold,
  col = "#A0352B",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = paste0(round(100 * expression_quantile), "th percentile: ",
                  round(expression_threshold, 2)),
  col = "#A0352B",
  lwd = 2,
  lty = 2,
  bty = "n"
)

dev.off()

keep_expr <- mean_expression > expression_threshold

gene_variance <- apply(exp_matrix, 1, var)

keep <- keep_expr & !is.na(gene_variance) & gene_variance > 0

write.table(gene_annot[!keep, ], file.path(matrix_dir, "removed_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

exp_matrix <- exp_matrix[keep, , drop = FALSE]
gene_annot <- gene_annot[keep, , drop = FALSE]

gene_names <- as.character(gene_annot$Hugo_Symbol)
no_name <- is.na(gene_names) | gene_names == ""
gene_names[no_name] <- paste0("ENTREZ_", gene_annot$Entrez_Gene_Id[no_name])
gene_annot$Hugo_Symbol <- gene_names

gene_annot$matrix_name <- make.unique(gene_names)
rownames(exp_matrix) <- gene_annot$matrix_name

exp_transposed <- t(exp_matrix)

use_gene <- !is.null(gene_of_interest) && !is.na(gene_of_interest) &&
  nzchar(gene_of_interest)
if (use_gene && !(gene_of_interest %in% gene_annot$Hugo_Symbol)) {
  message(gene_of_interest, " not found after filtering: gene-specific steps will be skipped")
  use_gene <- FALSE
}


message("[5/15] Computing per-chromosome correlations")

chromosomes <- unique(gene_annot$chromosome)
correl_list <- list()
annot_list <- list()

for (chromosome in chromosomes) {

  take <- gene_annot$chromosome == chromosome
  annot <- gene_annot[take, , drop = FALSE]
  annot$gene_index <- seq_len(nrow(annot))

  correl <- cor(exp_transposed[, take, drop = FALSE], method = "pearson")

  correl_list[[chromosome]] <- correl
  annot_list[[chromosome]] <- annot

  prefix <- file.path(matrix_dir, chromosome)
  write.table(correl, paste0(prefix, ".correlation.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(annot, paste0(prefix, ".annotation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message(chromosome, " genes: ", nrow(annot))
}

correlation_values <- unlist(
  lapply(
    correl_list[intersect(genome_chromosomes, names(correl_list))],
    function(x) x[upper.tri(x)]
  ),
  use.names = FALSE
)

correlation_values <- correlation_values[is.finite(correlation_values)]

svg(file.path(figures_dir, "gene_pair_correlation_distribution.svg"),
    width = 8, height = 5, family = font_family)

hist(
  correlation_values,
  breaks = 100,
  col = "#2b56a6",
  border = "white",
  xlim = c(-1, 1),
  main = "Distribution of gene-pair correlations",
  xlab = "Pearson correlation",
  ylab = "Number of gene pairs"
)

abline(
  v = c(-color_limit, color_limit),
  col = "#A0352B",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = paste("Color limits: ±", color_limit),
  col = "#A0352B",
  lwd = 2,
  lty = 2,
  bty = "n"
)

dev.off()

message("Correlations outside the color scale: ",
        round(100 * mean(abs(correlation_values) > color_limit), 1), "%")
print(quantile(abs(correlation_values), probs = c(0.90, 0.95, 0.99)))


message("[6/15] Scanning window sizes")

binsignal_by_window <- list()
thresholds <- data.frame(window = integer(), threshold = numeric())
cod_counts <- data.frame(window = integer(), chromosome = character(),
                         n_genes_chr = integer(), n_COD = integer())
all_cods <- data.frame(window = integer(), chromosome = character(),
                       START = integer(), END = integer(), n_genes = integer(),
                       first_gene_name = character(), last_gene_name = character(),
                       start_bp = numeric(), end_bp = numeric(), size_kb = numeric(),
                       intra_coexpression = numeric())

for (window.size in window_sizes) {

  w_name <- paste0("w", window.size)
  binsignal_by_window[[w_name]] <- list()
  genome_values <- c()

  for (chromosome in chromosomes) {

    correl <- correl_list[[chromosome]]
    n_bins <- nrow(correl)

    mean.cf <- rep(NA_real_, n_bins)
    pvalue <- rep(NA_real_, n_bins)

    if (n_bins >= 2) {
      for (i in 1:(n_bins - 1)) {
        lowerbound <- max(1, i - window.size + 1)
        upperbound <- min(i + window.size, n_bins)
        diamond <- correl[lowerbound:i, (i + 1):upperbound, drop = FALSE]
        mean.cf[i] <- mean(diamond)
      }
    }

    if (n_bins >= 2 * window.size + 1) {
      for (i in window.size:(n_bins - window.size - 1)) {
        upstream <- mean.cf[(i - window.size + 1):i]
        downstream <- mean.cf[(i + 1):(i + window.size)]
        ttest <- try(t.test(upstream, downstream, var.equal = TRUE), silent = TRUE)
        if (!inherits(ttest, "try-error")) {
          pvalue[i] <- ttest$p.value
        }
      }
    }

    binsignal_by_window[[w_name]][[chromosome]] <- data.frame(
      mean.cf = mean.cf, pvalue = pvalue,
      matrix_name = annot_list[[chromosome]]$matrix_name)

    if (chromosome %in% genome_chromosomes) {
      genome_values <- c(genome_values, mean.cf)
    }
  }

  threshold <- mean(genome_values, na.rm = TRUE)
  thresholds <- rbind(thresholds, data.frame(window = window.size,
                                             threshold = threshold))

  edge <- max(window.size, 4)

  for (chromosome in chromosomes) {

    correl <- correl_list[[chromosome]]
    annot <- annot_list[[chromosome]]
    mean.cf <- binsignal_by_window[[w_name]][[chromosome]]$mean.cf
    pvalue <- binsignal_by_window[[w_name]][[chromosome]]$pvalue
    n_bins <- length(mean.cf)

    pvalue[is.na(pvalue)] <- 1

    domain <- rep("", n_bins)
    i <- 2

    while (i < n_bins - edge) {

      if (mean.cf[i] > threshold &&
          sum(mean.cf[i:(i + 3)] > threshold) > 2 &&
          pvalue[i - 1] < 0.05) {

        domain[i - 1] <- "START"

        while (i < n_bins - edge &&
               any(mean.cf[(i + 1):(i + 4)] > threshold) &&
               (any(mean.cf[(i + 1):(i + 2)] > threshold) || pvalue[i] > 0.05)) {

          i <- i + 1
          domain[i - 1] <- "DOMAIN"
          domain[i] <- "DOMAIN"

          if (mean.cf[i + 1] < threshold && pvalue[i + 1] < 0.05) {
            domain[i + 1] <- "END"
          } else {
            domain[i] <- "END"
          }
        }
      }
      i <- i + 1
    }

    starts <- which(domain == "START")
    ends <- which(domain == "END")
    n_cod <- 0

    for (first_gene in starts) {

      possible_ends <- ends[ends > first_gene]
      if (length(possible_ends) == 0) next
      last_gene <- possible_ends[1]

      near_end <- max(1, last_gene - 2):min(n_bins, last_gene + 3)
      significant_end <- any(pvalue[near_end] < 0.05)
      number_of_genes <- last_gene - first_gene + 1

      if (significant_end && number_of_genes >= 4) {

        block <- correl[first_gene:last_gene, first_gene:last_gene]
        intra <- mean(block[upper.tri(block)])

        start_bp <- min(annot$start[first_gene:last_gene])
        end_bp <- max(annot$end[first_gene:last_gene])

        all_cods <- rbind(all_cods, data.frame(
          window = window.size, chromosome = chromosome,
          START = first_gene, END = last_gene, n_genes = number_of_genes,
          first_gene_name = annot$matrix_name[first_gene],
          last_gene_name = annot$matrix_name[last_gene],
          start_bp = start_bp, end_bp = end_bp,
          size_kb = (end_bp - start_bp + 1) / 1000,
          intra_coexpression = intra))
        n_cod <- n_cod + 1
      }
    }

    cod_counts <- rbind(cod_counts, data.frame(window = window.size,
                                               chromosome = chromosome,
                                               n_genes_chr = n_bins,
                                               n_COD = n_cod))
  }

  message("window: ", window.size, "  threshold: ", round(threshold, 4),
          "  COD total: ", sum(cod_counts$n_COD[cod_counts$window == window.size]))
}

cod_counts$chromosome <- factor(cod_counts$chromosome, levels = chromosomes)
cod_wide <- as.data.frame.matrix(xtabs(n_COD ~ chromosome + window, data = cod_counts))
colnames(cod_wide) <- paste0("w", colnames(cod_wide))
cod_wide <- rbind(cod_wide, total = colSums(cod_wide))
print(cod_wide)

write.table(cod_wide, file.path(window_dir, "COD_count_by_chromosome.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(all_cods, file.path(window_dir, "all_COD_all_windows.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(thresholds, file.path(window_dir, "thresholds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

window_summary <- data.frame()
for (window.size in window_sizes) {
  for (set in c("chr1", "genome")) {
    if (set == "chr1") {
      take <- all_cods$window == window.size & all_cods$chromosome == "chr1"
    } else {
      take <- all_cods$window == window.size & all_cods$chromosome %in% genome_chromosomes
    }
    cods <- all_cods[take, ]
    window_summary <- rbind(window_summary, data.frame(
      window = window.size,
      set = set,
      threshold = thresholds$threshold[thresholds$window == window.size],
      n_COD = nrow(cods),
      mean_n_genes = mean(cods$n_genes),
      median_n_genes = median(cods$n_genes),
      mean_intra_coexpression = mean(cods$intra_coexpression)))
  }
}
print(window_summary)
write.table(window_summary, file.path(window_dir, "window_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

chr1_summary <- window_summary[window_summary$set == "chr1", ]
if (!any(is.finite(chr1_summary$mean_intra_coexpression))) {
  stop("No COD found on chr1 for any window: cannot select the best window")
}
best_window <- chr1_summary$window[which.max(chr1_summary$mean_intra_coexpression)]

window.size <- best_window
w_name <- paste0("w", best_window)
threshold <- thresholds$threshold[thresholds$window == best_window]

message("Best window (chr1): ", best_window, "  threshold: ", round(threshold, 4))


message("[7/15] Writing final COD tables")

final_cods <- all_cods[all_cods$window == best_window, ]

for (chromosome in chromosomes) {
  prefix <- file.path(matrix_dir, chromosome)
  write.table(binsignal_by_window[[w_name]][[chromosome]],
              paste0(prefix, ".binsignal.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(final_cods[final_cods$chromosome == chromosome, ],
              paste0(prefix, ".COD.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
write.table(final_cods, file.path(matrix_dir, "all_COD.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

chr_summary <- data.frame(chromosome = character(), mean_binsignal = numeric())
for (chromosome in chromosomes) {
  current_mean <- mean(binsignal_by_window[[w_name]][[chromosome]]$mean.cf, na.rm = TRUE)
  if (is.nan(current_mean)) current_mean <- NA_real_
  chr_summary <- rbind(chr_summary,
                       data.frame(chromosome = chromosome, mean_binsignal = current_mean))
}
write.table(chr_summary, file.path(matrix_dir, "mean_binsignal_by_chromosome.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(window = best_window, mean_binsignal = threshold),
            file.path(matrix_dir, "mean_binsignal_genome.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_chromosomes <- chromosome_order[chromosome_order %in% chromosomes &
                                       chromosome_order != "chrM"]
chr_stats <- data.frame()
for (chromosome in plot_chromosomes) {
  cods <- final_cods[final_cods$chromosome == chromosome, ]
  n_genes_chr <- nrow(annot_list[[chromosome]])
  n_cod <- nrow(cods)
  chr_stats <- rbind(chr_stats, data.frame(
    chromosome = chromosome,
    n_genes_chr = n_genes_chr,
    n_COD = n_cod,
    COD_per_100_genes = 100 * n_cod / n_genes_chr,
    mean_genes_per_COD = if (n_cod > 0) mean(cods$n_genes) else NA,
    sem_genes_per_COD = if (n_cod > 1) sd(cods$n_genes) / sqrt(n_cod) else NA))
}
print(chr_stats)
write.table(chr_stats, file.path(matrix_dir, "COD_stats_by_chromosome.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


message("[8/15] Drawing heatmaps with binsignal")

for (chromosome in chromosomes) {

  correl <- correl_list[[chromosome]]
  n_bins <- nrow(correl)

  if (n_bins < 2) next

  mean.cf <- binsignal_by_window[[w_name]][[chromosome]]$mean.cf
  COD <- final_cods[final_cods$chromosome == chromosome, ]

  correl_plot <- pmin(pmax(correl, -color_limit), color_limit)

  svg(file.path(figures_dir, paste0("heatmap_", chromosome, ".svg")),
      width = 10, height = 11, family = font_family)

  layout(matrix(c(1, 2), ncol = 1), heights = c(4, 1))
  par(mar = c(0.5, 4, 4, 1))

  image(1:n_bins, 1:n_bins, t(correl_plot[n_bins:1, , drop = FALSE]),
        col = colors, zlim = c(-color_limit, color_limit), useRaster = TRUE,
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(0.5, n_bins + 0.5), ylim = c(0.5, n_bins + 0.5),
        xaxs = "i", yaxs = "i")
  box()

  plot_title <- paste0(chromosome,
                       "   |   genes: ", n_bins,
                       "   |   window: ", window.size,
                       "   |   COD: ", nrow(COD))
  mtext(plot_title, side = 3, line = 1.2, adj = 0, cex = 1.2, font = 2)

  par(xpd = NA)
  x0 <- grconvertX(0.86, "npc", "user")
  x1 <- grconvertX(1.00, "npc", "user")
  y0 <- grconvertY(1.015, "npc", "user")
  y1 <- grconvertY(1.030, "npc", "user")
  x_steps <- seq(x0, x1, length.out = length(colors) + 1)
  rect(x_steps[-length(x_steps)], y0, x_steps[-1], y1, col = colors, border = NA)
  rect(x0, y0, x1, y1, border = "black", lwd = 0.5)
  text(c(x0, (x0 + x1) / 2, x1), y1, labels = c(-color_limit, 0, color_limit),
       pos = 3, cex = 0.6, offset = 0.2)
  text(x0, (y0 + y1) / 2, labels = "Pearson r", pos = 2, cex = 0.6, offset = 0.3)
  par(xpd = FALSE)

  bar_colors <- rep("grey60", n_bins)
  for (k in seq_len(nrow(COD))) {
    first_gene <- COD$START[k]
    last_gene <- COD$END[k]
    rect(first_gene - 0.5, n_bins - last_gene + 0.5,
         last_gene + 0.5, n_bins - first_gene + 1.5, border = "black")
    bar_colors[first_gene:last_gene] <- "firebrick"
  }

  par(mar = c(4, 4, 0.5, 1))
  plot(1:n_bins, mean.cf, type = "h", col = bar_colors,
       xlim = c(0.5, n_bins + 0.5), ylim = c(-1, 1), xaxs = "i",
       xlab = "Gene index (genomic order)", ylab = "Binsignal")
  abline(h = 0, col = "grey80")
  abline(h = threshold, col = "#4fd138", lty = 2)

  dev.off()
}


if (use_gene) {
  message("[9/15] Zoomed COD around ", gene_of_interest)

  gene_chr <- gene_annot$chromosome[gene_annot$Hugo_Symbol == gene_of_interest]
  gene_chr <- gene_chr[1]

  annot <- annot_list[[gene_chr]]
  gene_index <- which(annot$Hugo_Symbol == gene_of_interest)[1]

  chr_cods <- final_cods[final_cods$chromosome == gene_chr, ]
  in_cod <- which(chr_cods$START <= gene_index & chr_cods$END >= gene_index)

  if (length(in_cod) > 0) {
    cod_start <- chr_cods$START[in_cod[1]]
    cod_end <- chr_cods$END[in_cod[1]]
    zoom_title <- paste0(gene_of_interest, " COD: ", gene_chr, ":",
                         format(chr_cods$start_bp[in_cod[1]], big.mark = ","), "-",
                         format(chr_cods$end_bp[in_cod[1]], big.mark = ","),
                         "  |  ", chr_cods$n_genes[in_cod[1]], " genes",
                         "  |  ", round(chr_cods$size_kb[in_cod[1]]), " kb",
                         "  |  mean r = ", round(chr_cods$intra_coexpression[in_cod[1]], 2),
                         "  |  window: ", window.size)
  } else {
    message(gene_of_interest, " is not in any COD; drawing its neighbourhood")
    cod_start <- gene_index
    cod_end <- gene_index
    zoom_title <- paste0(gene_of_interest, " (", gene_chr, "): no COD",
                         "  |  window: ", window.size)
  }

  zoom_from <- max(1, cod_start - zoom_flank)
  zoom_to <- min(nrow(annot), cod_end + zoom_flank)
  zoom_genes <- zoom_from:zoom_to
  n_zoom <- length(zoom_genes)

  correl_zoom <- correl_list[[gene_chr]][zoom_genes, zoom_genes]
  correl_zoom <- pmin(pmax(correl_zoom, -color_limit), color_limit)

  color_index <- round((correl_zoom + color_limit) / (2 * color_limit) *
                         (length(colors) - 1)) + 1

  zoom_table <- annot[zoom_genes, c("gene_index", "Hugo_Symbol", "chromosome",
                                    "start", "end")]
  zoom_table$in_COD <- zoom_genes >= cod_start & zoom_genes <= cod_end & length(in_cod) > 0
  zoom_table$r_with_gene <- correl_list[[gene_chr]][zoom_genes, gene_index]
  write.table(zoom_table,
              file.path(figures_dir, paste0(gene_of_interest, "_COD_genes.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  zoom_width <- max(7, n_zoom * 0.3 + 1)
  svg(file.path(figures_dir, paste0(gene_of_interest, "_COD_zoom.svg")),
      width = zoom_width, height = zoom_width * 0.5 + 2.2, family = font_family)
  par(mar = c(7, 0.5, 3, 0.5))
  plot.new()
  plot.window(xlim = c(0.5, n_zoom + 0.5), ylim = c(0, n_zoom / 2 + 0.5), asp = 1)

  for (a in 1:n_zoom) {
    for (b in a:n_zoom) {
      cx <- (a + b) / 2
      cy <- (b - a) / 2
      if (a == b) {
        polygon(c(cx - 0.5, cx, cx + 0.5), c(0, 0.5, 0),
                col = colors[color_index[a, b]], border = "white", lwd = 0.6)
      } else {
        polygon(c(cx - 0.5, cx, cx + 0.5, cx), c(cy, cy + 0.5, cy, cy - 0.5),
                col = colors[color_index[a, b]], border = "white", lwd = 0.6)
      }
    }
  }

  if (length(in_cod) > 0) {
    s <- cod_start - zoom_from + 1
    e <- cod_end - zoom_from + 1
    polygon(c(s - 0.5, (s + e) / 2, e + 0.5), c(0, (e - s) / 2 + 0.5, 0),
            border = "black", lwd = 1.8, col = NA)
  }

  label_colors <- ifelse(annot$Hugo_Symbol[zoom_genes] == gene_of_interest,
                         "firebrick", "black")
  label_font <- ifelse(annot$Hugo_Symbol[zoom_genes] == gene_of_interest, 2, 1)
  text(1:n_zoom, -0.2, labels = annot$Hugo_Symbol[zoom_genes],
       srt = 90, adj = 1, xpd = NA, cex = 0.8,
       col = label_colors, font = label_font)

  mtext(zoom_title, side = 3, line = 1, adj = 0, cex = 1, font = 2)

  x1 <- n_zoom + 0.5
  x0 <- x1 - min(5, n_zoom / 4)
  y1 <- n_zoom / 2 + 0.5
  y0 <- y1 - 0.3
  x_steps <- seq(x0, x1, length.out = length(colors) + 1)
  rect(x_steps[-length(x_steps)], y0, x_steps[-1], y1, col = colors, border = NA)
  rect(x0, y0, x1, y1, border = "black", lwd = 0.5)
  text(c(x0, (x0 + x1) / 2, x1), y0, labels = c(-color_limit, 0, color_limit),
       pos = 1, cex = 0.6, offset = 0.2)
  text(x0, (y0 + y1) / 2, labels = "Pearson r", pos = 2, cex = 0.6, offset = 0.3)

  dev.off()
}


message("[10/15] Barplots by chromosome")

chr_labels <- sub("chr", "", chr_stats$chromosome)
bar_color <- "#3b5ba5"

svg(file.path(figures_dir, "COD_by_chromosome.svg"),
    width = 10, height = 10, family = font_family)
par(mfrow = c(3, 1), mar = c(4, 5, 2.5, 1), las = 1)

mids <- barplot(chr_stats$n_COD, names.arg = chr_labels, col = bar_color,
                border = NA, ylim = c(0, max(chr_stats$n_COD, 1) * 1.15),
                xlab = "Chromosome", ylab = "Number of CODs")
text(mids, chr_stats$n_COD, labels = chr_stats$n_COD, pos = 3, cex = 0.8, offset = 0.2)
title(main = "A", adj = 0)

density_max <- max(chr_stats$COD_per_100_genes, 0.1, na.rm = TRUE)
mids <- barplot(chr_stats$COD_per_100_genes, names.arg = chr_labels, col = bar_color,
                border = NA, ylim = c(0, density_max * 1.15),
                xlab = "Chromosome", ylab = "CODs per 100 genes")
abline(h = 100 * nrow(final_cods[final_cods$chromosome %in% genome_chromosomes, ]) /
         sum(chr_stats$n_genes_chr[chr_stats$chromosome %in% genome_chromosomes]),
       col = "darkorange3", lty = 2)
title(main = "B", adj = 0)

upper <- chr_stats$mean_genes_per_COD + ifelse(is.na(chr_stats$sem_genes_per_COD), 0,
                                               chr_stats$sem_genes_per_COD)
mids <- barplot(chr_stats$mean_genes_per_COD, names.arg = chr_labels, col = bar_color,
                border = NA, ylim = c(0, max(upper, 1, na.rm = TRUE) * 1.15),
                xlab = "Chromosome", ylab = "Genes per COD (mean ± SEM)")
has_sem <- !is.na(chr_stats$sem_genes_per_COD) & chr_stats$sem_genes_per_COD > 0
arrows(mids[has_sem], chr_stats$mean_genes_per_COD[has_sem] - chr_stats$sem_genes_per_COD[has_sem],
       mids[has_sem], chr_stats$mean_genes_per_COD[has_sem] + chr_stats$sem_genes_per_COD[has_sem],
       angle = 90, code = 3, length = 0.03)
title(main = "C", adj = 0)

dev.off()


message("[11/15] COD distributions")

genome_cods <- final_cods[final_cods$chromosome %in% genome_chromosomes, ]
fill_color <- "#7ff5f5"
line_color <- "#1f2bff"

svg(file.path(figures_dir, "COD_distributions.svg"),
    width = 13, height = 4.3, family = font_family)
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1), las = 1)

max_genes <- max(genome_cods$n_genes, 4)
gene_counts <- as.numeric(table(factor(genome_cods$n_genes, levels = 0:max_genes)))
plot(0:max_genes, gene_counts, type = "n", bty = "l",
     xlab = "Number of genes per COD", ylab = "Number of CODs",
     ylim = c(0, max(gene_counts, 1) * 1.1))
polygon(c(0, 0:max_genes, max_genes), c(0, gene_counts, 0),
        col = fill_color, border = NA)
lines(0:max_genes, gene_counts, col = line_color, lwd = 1.5)
legend("topright", bty = "n",
       legend = c(paste0("Average = ", round(mean(genome_cods$n_genes), 1), " genes"),
                  paste0("Median = ", median(genome_cods$n_genes), " genes")))
title(main = "A", adj = 0, cex.main = 1.6)

size_breaks <- seq(0, ceiling(max(genome_cods$size_kb, 100) / 100) * 100 + 100, by = 100)
size_hist <- hist(genome_cods$size_kb, breaks = size_breaks, plot = FALSE)
plot(size_hist$mids, size_hist$counts, type = "n", bty = "l",
     xlab = "kb per COD", ylab = "Number of CODs",
     xlim = c(0, max(size_breaks)), ylim = c(0, max(size_hist$counts, 1) * 1.1))
polygon(c(0, size_hist$mids, max(size_hist$mids)), c(0, size_hist$counts, 0),
        col = fill_color, border = NA)
lines(size_hist$mids, size_hist$counts, col = line_color, lwd = 1.5)
legend("topright", bty = "n",
       legend = c(paste0("Average = ", round(mean(genome_cods$size_kb) / 1000, 2), " Mb"),
                  paste0("Median = ", round(median(genome_cods$size_kb) / 1000, 2), " Mb")))
title(main = "B", adj = 0, cex.main = 1.6)

coex_breaks <- seq(floor(min(genome_cods$intra_coexpression, 0) * 50) / 50,
                   ceiling(max(genome_cods$intra_coexpression, 0.1) * 50) / 50 + 0.02,
                   by = 0.02)
hist(genome_cods$intra_coexpression, breaks = coex_breaks,
     col = "#e8432e", border = "black", main = "",
     xlab = "Average intra-COD coexpression", ylab = "Number of CODs")
legend("topright", bty = "n",
       legend = c(paste0("Average = ", round(mean(genome_cods$intra_coexpression), 2)),
                  paste0("Median = ", round(median(genome_cods$intra_coexpression), 2))))
title(main = "C", adj = 0, cex.main = 1.6)

dev.off()


message("[12/15] Window selection figures")

point_color <- "#2b3bd6"

for (set in c("chr1", "genome")) {

  if (set == "chr1") {
    set_cods <- all_cods[all_cods$chromosome == "chr1", ]
  } else {
    set_cods <- all_cods[all_cods$chromosome %in% genome_chromosomes, ]
  }
  set_summary <- window_summary[window_summary$set == set, ]
  window_factor <- factor(set_cods$window, levels = window_sizes)

  svg(file.path(window_dir, paste0("window_selection_", set, ".svg")),
      width = 5, height = 11, family = font_family)
  par(mfrow = c(3, 1), mar = c(4.5, 5.5, 2.5, 1), las = 1, font.lab = 2)

  stripchart(set_cods$intra_coexpression ~ window_factor, vertical = TRUE,
             method = "jitter", jitter = 0.2, pch = 16, cex = 0.5, col = point_color,
             xlab = "Window size", ylab = "Average intra-COD\ncoexpression")
  for (k in seq_along(window_sizes)) {
    values <- set_cods$intra_coexpression[set_cods$window == window_sizes[k]]
    if (length(values) == 0) next
    segments(k - 0.3, median(values), k + 0.3, median(values), lwd = 3,
             col = ifelse(window_sizes[k] == best_window, "firebrick", "black"))
  }
  title(main = "A", adj = 0, cex.main = 1.5)

  plot(set_summary$window, set_summary$n_COD, pch = 19, cex = 1.5, bty = "l",
       xlim = c(min(window_sizes) - 0.5, max(window_sizes) + 0.5),
       ylim = c(0, max(set_summary$n_COD, 1) * 1.2), xaxt = "n",
       xlab = "Window size", ylab = "Number of CODs detected")
  axis(1, at = window_sizes)
  title(main = "B", adj = 0, cex.main = 1.5)

  stripchart(set_cods$n_genes ~ window_factor, vertical = TRUE,
             method = "jitter", jitter = 0.2, pch = 16, cex = 0.5, col = point_color,
             xlab = "Window size", ylab = "COD size (number of genes)")
  for (k in seq_along(window_sizes)) {
    values <- set_cods$n_genes[set_cods$window == window_sizes[k]]
    if (length(values) == 0) next
    m <- mean(values)
    segments(k - 0.3, m, k + 0.3, m, lwd = 3)
    if (length(values) > 1) {
      sem <- sd(values) / sqrt(length(values))
      if (sem > 0) arrows(k, m - sem, k, m + sem, angle = 90, code = 3,
                          length = 0.06, lwd = 2)
    }
  }
  title(main = "C", adj = 0, cex.main = 1.5)

  dev.off()
}


message("[13/15] Mean correlation between CODs")

cod_x_cod_matrix_dir <- file.path(matrix_dir, "COD_x_COD")
cod_x_cod_figures_dir <- file.path(figures_dir, "COD_x_COD")
dir.create(cod_x_cod_matrix_dir, showWarnings = FALSE)
dir.create(cod_x_cod_figures_dir, showWarnings = FALSE)

inter_cods <- final_cods[final_cods$chromosome %in% genome_chromosomes, ]
inter_cods <- inter_cods[order(match(inter_cods$chromosome, chromosome_order),
                               inter_cods$START), ]
inter_cods <- as.data.frame(inter_cods)
rownames(inter_cods) <- NULL
n_inter <- nrow(inter_cods)

if (n_inter == 0) stop("No CODs available for the coexpression analysis")

inter_cods$COD_ID <- paste0(inter_cods$chromosome, "_COD_",
                            ave(seq_len(n_inter), inter_cods$chromosome,
                                FUN = seq_along))
cod_ids <- inter_cods$COD_ID

write.table(inter_cods, file.path(cod_x_cod_matrix_dir, "COD_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cod_matrix <- matrix(NA_real_, nrow = n_inter, ncol = n_inter,
                     dimnames = list(cod_ids, cod_ids))

chromosomes_with_cod <- unique(inter_cods$chromosome)

for (chr_a_index in seq_along(chromosomes_with_cod)) {
  chr_a <- chromosomes_with_cod[chr_a_index]
  cod_a_indices <- which(inter_cods$chromosome == chr_a)
  annot_a <- annot_list[[chr_a]]

  for (chr_b_index in chr_a_index:length(chromosomes_with_cod)) {
    chr_b <- chromosomes_with_cod[chr_b_index]
    cod_b_indices <- which(inter_cods$chromosome == chr_b)
    annot_b <- annot_list[[chr_b]]

    if (chr_a == chr_b) {
      gene_cor <- correl_list[[chr_a]]
    } else {
      genes_a <- unique(unlist(lapply(cod_a_indices, function(k) {
        annot_a$matrix_name[inter_cods$START[k]:inter_cods$END[k]]
      })))
      genes_b <- unique(unlist(lapply(cod_b_indices, function(k) {
        annot_b$matrix_name[inter_cods$START[k]:inter_cods$END[k]]
      })))
      gene_cor <- cor(exp_transposed[, genes_a, drop = FALSE],
                      exp_transposed[, genes_b, drop = FALSE],
                      method = "pearson")
    }

    for (a in cod_a_indices) {
      genes_a <- annot_a$matrix_name[inter_cods$START[a]:inter_cods$END[a]]

      for (b in cod_b_indices) {
        if (b < a) next
        genes_b <- annot_b$matrix_name[inter_cods$START[b]:inter_cods$END[b]]
        block <- gene_cor[genes_a, genes_b, drop = FALSE]

        if (a == b) {
          value <- mean(block[upper.tri(block)])
        } else {
          value <- mean(block)
        }
        cod_matrix[a, b] <- value
        cod_matrix[b, a] <- value
      }
    }
    message("Done: ", chr_a, " and ", chr_b)
  }
}

stopifnot(all(is.finite(cod_matrix)))
for (chr in chromosomes_with_cod) {
  x <- inter_cods[inter_cods$chromosome == chr, ]
  if (nrow(x) > 1 && any(x$START[-1] <= x$END[-nrow(x)])) {
    stop("Overlapping CODs on ", chr, ": check the boundary calls")
  }
}

write.table(cod_matrix, file.path(cod_x_cod_matrix_dir, "genome.mean_correlation.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

for (chromosome in chromosomes_with_cod) {
  take <- inter_cods$chromosome == chromosome
  chr_cod_matrix <- cod_matrix[take, take, drop = FALSE]
  write.table(chr_cod_matrix,
              file.path(cod_x_cod_matrix_dir, paste0(chromosome, ".mean_correlation.tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}

pair_index <- which(upper.tri(cod_matrix), arr.ind = TRUE)
cod_pairs <- data.frame(
  COD_1 = cod_ids[pair_index[, 1]],
  COD_2 = cod_ids[pair_index[, 2]],
  chromosome_1 = inter_cods$chromosome[pair_index[, 1]],
  chromosome_2 = inter_cods$chromosome[pair_index[, 2]],
  n_genes_1 = inter_cods$n_genes[pair_index[, 1]],
  n_genes_2 = inter_cods$n_genes[pair_index[, 2]],
  mean_r = cod_matrix[pair_index]
)
cod_pairs$comparison <- ifelse(cod_pairs$chromosome_1 == cod_pairs$chromosome_2,
                               "within_chromosome", "between_chromosomes")

write.table(cod_pairs, file.path(cod_x_cod_matrix_dir, "COD_pairs_mean_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cod_pairs[cod_pairs$comparison == "within_chromosome", ],
            file.path(cod_x_cod_matrix_dir, "COD_pairs_within_chromosomes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cod_pairs[cod_pairs$comparison == "between_chromosomes", ],
            file.path(cod_x_cod_matrix_dir, "COD_pairs_between_chromosomes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (set in c("genome", chromosomes_with_cod)) {
  if (set == "genome") {
    current_matrix <- cod_matrix
    current_cod <- inter_cods
  } else {
    take <- inter_cods$chromosome == set
    current_matrix <- cod_matrix[take, take, drop = FALSE]
    current_cod <- inter_cods[take, , drop = FALSE]
  }

  n_set <- nrow(current_matrix)
  plot_matrix <- pmin(pmax(current_matrix, -color_limit), color_limit)
  dim(plot_matrix) <- c(n_set, n_set)
  pix <- max(1, floor(1200 / n_set))
  set_rows <- rep(n_set:1, each = pix)
  set_cols <- rep(seq_len(n_set), each = pix)
  set_coords <- seq(0.5 + 0.5 / pix, n_set + 0.5 - 0.5 / pix,
                    length.out = n_set * pix)

  svg(file.path(cod_x_cod_figures_dir, paste0(set, ".mean_correlation.svg")),
      width = 10, height = 10, family = font_family)
  par(mar = c(3, 3, 5, 2))
  image(set_coords, set_coords,
        t(plot_matrix[set_rows, set_cols, drop = FALSE]),
        col = colors, zlim = c(-color_limit, color_limit), useRaster = TRUE,
        axes = FALSE, xlab = "", ylab = "",
        xlim = c(0.5, n_set + 0.5), ylim = c(0.5, n_set + 0.5),
        xaxs = "i", yaxs = "i")
  box()

  if (set == "genome") {
    chr_runs <- rle(current_cod$chromosome)
    chr_ends <- cumsum(chr_runs$lengths)
    chr_mids <- chr_ends - chr_runs$lengths / 2 + 0.5
    if (length(chr_ends) > 1) {
      abline(v = head(chr_ends, -1) + 0.5, col = "grey30", lwd = 0.4)
      abline(h = n_set - head(chr_ends, -1) + 0.5, col = "grey30", lwd = 0.4)
    }
    axis(3, at = chr_mids, labels = sub("chr", "", chr_runs$values),
         tick = FALSE, line = -0.8, cex.axis = 0.6)
    axis(2, at = n_set - chr_mids + 1,
         labels = sub("chr", "", chr_runs$values),
         tick = FALSE, line = -0.8, cex.axis = 0.6, las = 1)
  }
  title(main = paste0(set, ": mean Pearson r between COD"), line = 2)
  mtext(paste0("Blue: <= -", color_limit, "   White: 0   Red: >= ", color_limit),
        side = 1, line = 1)
  dev.off()
}


if (use_gene) {
  message("[14/15] CODs and genes related to ", gene_of_interest)

  gene_dir <- file.path(matrix_dir, paste0(gene_of_interest, "_analysis"))
  dir.create(gene_dir, showWarnings = FALSE)
  set.seed(random_seed)

  goi_rows <- which(gene_annot$Hugo_Symbol == gene_of_interest)
  if (length(goi_rows) > 1) {
    warning("Several rows found for ", gene_of_interest, ": using the first one")
  }
  goi_row <- goi_rows[1]
  goi_name <- gene_annot$matrix_name[goi_row]
  goi_chr <- gene_annot$chromosome[goi_row]
  goi_position <- match(goi_name, annot_list[[goi_chr]]$matrix_name)

  gene_r <- as.numeric(cor(exp_transposed[, goi_name], exp_transposed,
                           method = "pearson"))
  gene_results <- data.frame(
    matrix_name = gene_annot$matrix_name,
    gene = gene_annot$Hugo_Symbol,
    Entrez_Gene_Id = gene_annot$Entrez_Gene_Id,
    chromosome = gene_annot$chromosome,
    start = gene_annot$start,
    end = gene_annot$end,
    correlation_with_gene = gene_r,
    stringsAsFactors = FALSE
  )
  gene_results <- gene_results[gene_results$matrix_name != goi_name, ]
  gene_results <- gene_results[order(-gene_results$correlation_with_gene,
                                     gene_results$matrix_name), ]
  top_genes <- head(gene_results, 100)
  top_genes$rank <- seq_len(nrow(top_genes))
  top_genes$in_top_selected_COD <- NA
  top_genes$selected_COD_ID <- NA_character_

  goi_cod_rows <- which(inter_cods$chromosome == goi_chr &
                          inter_cods$START <= goi_position &
                          inter_cods$END >= goi_position)

  if (length(goi_cod_rows) != 1) {
    warning("CODs found for ", gene_of_interest, ": ", length(goi_cod_rows),
            ". The COD cannot be chosen unambiguously; the gene table is saved ",
            "and the COD comparison is skipped.")
  } else {
    goi_cod <- goi_cod_rows[1]
    goi_cod_id <- cod_ids[goi_cod]
    goi_annot <- annot_list[[goi_chr]]
    goi_genes <- goi_annot$matrix_name[
      inter_cods$START[goi_cod]:inter_cods$END[goi_cod]
    ]

    candidates <- setdiff(seq_len(n_inter), goi_cod)
    overlaps_goi <- logical(length(candidates))
    for (k in seq_along(candidates)) {
      candidate <- candidates[k]
      candidate_annot <- annot_list[[inter_cods$chromosome[candidate]]]
      candidate_genes <- candidate_annot$matrix_name[
        inter_cods$START[candidate]:inter_cods$END[candidate]
      ]
      overlaps_goi[k] <- any(candidate_genes %in% goi_genes)
    }
    candidates <- candidates[!overlaps_goi]

    if (length(candidates) == 0) {
      warning("No other non-overlapping CODs to compare with the COD of ", gene_of_interest)
    } else {
      observed <- cod_matrix[goi_cod, candidates]
      ordering <- order(-observed, cod_ids[candidates])
      candidates <- candidates[ordering]
      observed <- observed[ordering]
      same_candidates <- candidates[inter_cods$chromosome[candidates] == goi_chr]
      other_candidates <- candidates[inter_cods$chromosome[candidates] != goi_chr]
      top_indices <- c(head(same_candidates, 10), head(other_candidates, 10))
      n_top <- length(top_indices)
      top_observed <- cod_matrix[goi_cod, top_indices]

      goi_z <- scale(exp_transposed[, goi_genes, drop = FALSE])
      goi_profile <- rowMeans(goi_z)
      profiles <- matrix(NA_real_, nrow = nrow(exp_transposed),
                         ncol = length(candidates))
      for (k in seq_along(candidates)) {
        candidate <- candidates[k]
        candidate_annot <- annot_list[[inter_cods$chromosome[candidate]]]
        candidate_genes <- candidate_annot$matrix_name[
          inter_cods$START[candidate]:inter_cods$END[candidate]
        ]
        candidate_z <- scale(exp_transposed[, candidate_genes, drop = FALSE])
        profiles[, k] <- rowMeans(candidate_z)
      }
      profile_r <- as.numeric(crossprod(goi_profile, profiles) /
                                (nrow(exp_transposed) - 1))
      stopifnot(isTRUE(all.equal(as.numeric(observed), profile_r, tolerance = 1e-7)))

      permutation_max <- numeric(n_permutations)
      for (iteration in seq_len(n_permutations)) {
        shuffled <- sample(goi_profile, length(goi_profile), replace = FALSE)
        permuted_r <- as.numeric(crossprod(shuffled, profiles) /
                                   (nrow(exp_transposed) - 1))
        permutation_max[iteration] <- max(permuted_r)
      }
      p_adjusted <- numeric(n_top)
      for (k in seq_len(n_top)) {
        p_adjusted[k] <- (1 + sum(permutation_max >= top_observed[k])) /
          (n_permutations + 1)
      }

      top_cods <- inter_cods[top_indices, , drop = FALSE]
      top_cods$group <- ifelse(top_cods$chromosome == goi_chr,
                               "same_chromosome", "other_chromosomes")
      top_cods$rank <- ave(seq_len(n_top), top_cods$group, FUN = seq_along)
      top_cods$gene_of_interest_COD_ID <- goi_cod_id
      top_cods$mean_r_with_gene_COD <- top_observed
      top_cods$group_mean <- ifelse(top_cods$group == "same_chromosome",
                                    mean(cod_matrix[goi_cod, same_candidates]),
                                    mean(cod_matrix[goi_cod, other_candidates]))
      top_cods$difference_from_group_mean <- top_observed - top_cods$group_mean
      top_cods$p_maxT_one_sided <- p_adjusted
      top_cods$significant_p05 <- p_adjusted < 0.05
      top_cods$n_candidates <- length(candidates)
      top_cods$n_permutations <- n_permutations

      write.table(top_cods[top_cods$group == "same_chromosome", ],
                  file.path(gene_dir, paste0(gene_of_interest, "_top10_COD_same_chromosome.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(top_cods[top_cods$group == "other_chromosomes", ],
                  file.path(gene_dir, paste0(gene_of_interest, "_top10_COD_other_chromosomes.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      all_candidates <- inter_cods[candidates, , drop = FALSE]
      all_candidates$mean_r_with_gene_COD <- observed
      all_candidates$group <- ifelse(all_candidates$chromosome == goi_chr,
                                     "same_chromosome", "other_chromosomes")
      write.table(all_candidates,
                  file.path(gene_dir, paste0(gene_of_interest, "_all_candidate_COD.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      top_genes$in_top_selected_COD <- FALSE
      for (k in seq_len(n_top)) {
        candidate <- top_indices[k]
        candidate_annot <- annot_list[[inter_cods$chromosome[candidate]]]
        candidate_genes <- candidate_annot$matrix_name[
          inter_cods$START[candidate]:inter_cods$END[candidate]
        ]
        matches <- which(top_genes$matrix_name %in% candidate_genes)
        for (j in matches) {
          top_genes$in_top_selected_COD[j] <- TRUE
          if (is.na(top_genes$selected_COD_ID[j])) {
            top_genes$selected_COD_ID[j] <- cod_ids[candidate]
          } else {
            top_genes$selected_COD_ID[j] <- paste(top_genes$selected_COD_ID[j],
                                                  cod_ids[candidate], sep = ";")
          }
        }
      }
      print(top_cods[, c("rank", "COD_ID", "mean_r_with_gene_COD",
                         "p_maxT_one_sided", "significant_p05")])
    }
  }

  write.table(top_genes,
              file.path(gene_dir, paste0(gene_of_interest, "_top100_correlated_genes.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Genes saved: ", nrow(top_genes))
}


message("[15/15] Mirrored-boundary control plots")

control_figures_dir <- file.path(figures_dir, "mirrored_control")
control_matrix_dir <- file.path(matrix_dir, "mirrored_control")
dir.create(control_figures_dir, showWarnings = FALSE)
dir.create(control_matrix_dir, showWarnings = FALSE)

random_cods <- inter_cods
for (chr in chromosomes_with_cod) {
  ids <- which(inter_cods$chromosome == chr)
  n <- nrow(annot_list[[chr]])
  random_cods$START[ids] <- n + 1 - inter_cods$END[ids]
  random_cods$END[ids] <- n + 1 - inter_cods$START[ids]
}
real_intra <- inter_cods$intra_coexpression
random_intra <- numeric(n_inter)
for (k in seq_len(n_inter)) {
  chr <- random_cods$chromosome[k]
  ix <- random_cods$START[k]:random_cods$END[k]
  block <- correl_list[[chr]][ix, ix, drop = FALSE]
  random_intra[k] <- mean(block[upper.tri(block)])
}
write.table(data.frame(COD_ID = cod_ids, real = real_intra, mirrored = random_intra),
            file.path(control_matrix_dir, "intra_COD_real_mirrored.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
svg(file.path(control_figures_dir, "intra_COD_real_vs_mirrored.svg"), width = 7, height = 5)
hist(real_intra, breaks = 30, freq = FALSE, col = "#85c9e8", border = "white",
     xlim = c(-1, 1), xlab = "Mean intra-COD Pearson r", main = "Real and mirrored CODs")
if (length(random_intra) > 1 && sd(random_intra) > 0) lines(density(random_intra), col = "grey35", lwd = 2)
legend("topright", c("Real CODs", "Mirrored CODs"), fill = c("#85c9e8", NA),
       border = NA, lty = c(NA, 1), col = c("#85c9e8", "grey35"), bty = "n")
dev.off()

pair_values <- list(intra = numeric(), inter = numeric(), rest = numeric())
for (chr in intersect(genome_chromosomes, chromosomes)) {
  annot <- annot_list[[chr]]
  cods <- inter_cods[inter_cods$chromosome == chr, ]
  membership <- integer(nrow(annot))
  for (k in seq_len(nrow(cods))) membership[cods$START[k]:cods$END[k]] <- k
  indices <- which(upper.tri(correl_list[[chr]]), arr.ind = TRUE)
  a <- membership[indices[, 1]]
  b <- membership[indices[, 2]]
  values <- correl_list[[chr]][indices]
  pair_values$intra <- c(pair_values$intra, values[a > 0 & a == b])
  pair_values$inter <- c(pair_values$inter, values[a > 0 & b > 0 & a != b])
  pair_values$rest <- c(pair_values$rest, values[a == 0 | b == 0])
}
if (all(lengths(pair_values) > 0)) {
  svg(file.path(control_figures_dir, "gene_pair_correlations.svg"), width = 7, height = 5)
  boxplot(pair_values, outline = FALSE, col = c("#e87b73", "#92c5e8", "grey80"),
          ylab = "Gene-pair Pearson r", main = "Within-chromosome gene pairs")
  abline(h = 0, col = "grey60", lty = 2)
  dev.off()
}
write.table(data.frame(group = names(pair_values), n_pairs = lengths(pair_values),
                       median_r = sapply(pair_values, median)),
            file.path(control_matrix_dir, "gene_pair_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

within <- cod_pairs[cod_pairs$comparison == "within_chromosome", ]
if (nrow(within) > 0) {
  ia <- match(within$COD_1, cod_ids)
  ib <- match(within$COD_2, cod_ids)
  mirrored_r <- numeric(nrow(within))
  for (k in seq_len(nrow(within))) {
    chr <- within$chromosome_1[k]
    a <- random_cods$START[ia[k]]:random_cods$END[ia[k]]
    b <- random_cods$START[ib[k]]:random_cods$END[ib[k]]
    mirrored_r[k] <- mean(correl_list[[chr]][a, b, drop = FALSE])
  }
  within$mirrored_r <- mirrored_r
  center_a <- (inter_cods$start_bp[ia] + inter_cods$end_bp[ia]) / 2
  center_b <- (inter_cods$start_bp[ib] + inter_cods$end_bp[ib]) / 2
  within$distance_Mb <- abs(center_a - center_b) / 1e6
  write.table(within, file.path(control_matrix_dir, "within_COD_pairs_real_mirrored.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  svg(file.path(control_figures_dir, "inter_COD_distribution.svg"), width = 7, height = 5)
  hist(within$mean_r, breaks = 30, freq = FALSE, col = "#85c9e8", border = "white",
       xlim = c(-1, 1), xlab = "Mean inter-COD Pearson r", main = "Within-chromosome COD pairs")
  if (length(mirrored_r) > 1 && sd(mirrored_r) > 0) lines(density(mirrored_r), col = "grey35", lwd = 2)
  legend("topright", c("Real CODs", "Mirrored CODs"), fill = c("#85c9e8", NA),
         border = NA, lty = c(NA, 1), col = c("#85c9e8", "grey35"), bty = "n")
  dev.off()
  within <- within[order(within$distance_Mb), ]
  n_groups <- min(50, nrow(within))
  within$distance_bin <- ceiling(seq_len(nrow(within)) * n_groups / nrow(within))
  distance_summary <- aggregate(cbind(distance_Mb, mean_r, mirrored_r) ~ distance_bin,
                                data = within, FUN = mean)
  write.table(distance_summary, file.path(control_matrix_dir, "inter_COD_distance_bins.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  svg(file.path(control_figures_dir, "inter_COD_distance.svg"), width = 7, height = 5)
  plot(distance_summary$distance_Mb, distance_summary$mean_r, col = "firebrick", pch = 16,
       ylim = range(c(distance_summary$mean_r, distance_summary$mirrored_r)),
       xlab = "Distance between COD centers (Mb)", ylab = "Mean inter-COD Pearson r",
       main = "Coexpression and distance")
  points(distance_summary$distance_Mb, distance_summary$mirrored_r, col = "steelblue", pch = 16)
  legend("topright", c("Real CODs", "Mirrored CODs"), col = c("firebrick", "steelblue"), pch = 16, bty = "n")
  dev.off()
}

message("Done. Results are in: ", normalizePath(output_dir))
