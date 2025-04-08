#!/usr/bin/env Rscript

#' @title Genomic Feature Annotation and Visualization Tool
#' @description Create comprehensive genomic feature annotations and visualizations for ChIP-seq/ATAC-seq peaks
#' @param input_dir Directory containing narrowPeak files
#' @param output_dir Output directory for results
#' @param pattern File pattern to match in input directory
#' @param upstream Upstream distance from TSS for promoter annotation
#' @param downstream Downstream distance from TSS for promoter annotation
#' @param sample_regex Regular expression to extract sample name from filename
#' @param show_percentages Whether to show percentages in pie chart legends
#' @param save_default_piechart Whether to save default ChIPseeker pie charts
#' @author Surya Chhetri
#' @version 1.0
#' @date 2025-04-07

# Parse command line arguments
suppressPackageStartupMessages({
  library(optparse)
})

# Start timing
script_start_time <- Sys.time()

# Define command line options
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL,
              help="Directory containing narrowPeak files [required]"),
  make_option(c("-o", "--output_dir"), type="character", default="./annotation_results",
              help="Output directory for results [default: %default]"),
  make_option(c("-p", "--pattern"), type="character", default="*.narrowPeak",
              help="File pattern to match in input directory [default: %default]"),
  make_option(c("-u", "--upstream"), type="integer", default=3000,
              help="Upstream distance from TSS for promoter annotation [default: %default]"),
  make_option(c("-d", "--downstream"), type="integer", default=3000,
              help="Downstream distance from TSS for promoter annotation [default: %default]"),
  make_option(c("--sample_regex"), type="character", default="^(.*?)_peaks\\.narrowPeak$",
              help="Regular expression to extract sample name from filename [default: %default]"),
  make_option(c("--show_percentages"), type="logical", default=TRUE,
              help="Show percentages in pie chart legends [default: %default]"),
  make_option(c("--save_default_piechart"), type="logical", default=TRUE,
              help="Save default ChIPseeker pie charts [default: %default]")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list, 
                          description="Genomic annotation and visualization tool for ChIP-seq/ATAC-seq peaks",
                          usage="usage: %prog [options]")
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("Input directory must be specified", call.=FALSE)
}

# Load required libraries
packages <- c("GenomicRanges", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
              "org.Hs.eg.db", "ggplot2", "reshape2", "gridExtra", "grid", 
              "dplyr", "scales")

missing_packages <- c()
for (package in packages) {
  if (!suppressWarnings(suppressPackageStartupMessages(require(package, character.only = TRUE, quietly = TRUE)))) {
    missing_packages <- c(missing_packages, package)
  }
}

if (length(missing_packages) > 0) {
  stop("The following required packages are not installed: ", 
       paste(missing_packages, collapse = ", "), 
       "\nPlease install them with BiocManager::install() or install.packages()", call.=FALSE)
}

# Get input files from directory
input_dir <- normalizePath(opt$input_dir, mustWork=TRUE)
file_pattern <- opt$pattern
input_files <- list.files(path=input_dir, pattern=file_pattern, full.names=TRUE)

if (length(input_files) == 0) {
  stop("No files matching pattern '", file_pattern, "' found in directory: ", input_dir, call.=FALSE)
}

# Extract sample names from filenames using regex
sample_regex <- opt$sample_regex
sample_names <- sapply(basename(input_files), function(filename) {
  matches <- regmatches(filename, regexec(sample_regex, filename))[[1]]
  if (length(matches) > 1) {
    return(matches[2])
  } else {
    return(gsub("\\.[^.]*$", "", filename))
  }
})

# Print summary of files to be processed
cat("Processing", length(input_files), "files:\n")
for (i in 1:length(input_files)) {
  cat(sprintf("  %d. %s -> %s\n", i, basename(input_files[i]), sample_names[i]))
}

# Set other options from command line
output_dir <- opt$output_dir
upstream <- opt$upstream
downstream <- opt$downstream
show_percentages <- opt$show_percentages
save_default_piechart <- opt$save_default_piechart

# Define feature order and colors
feature_order <- c(
    "Promoter (<=1kb)",
    "Promoter (1-2kb)",
    "Promoter (2-3kb)",
    "5' UTR",
    "3' UTR",
    "1st Exon",
    "Other Exon",
    "1st Intron",
    "Other Intron",
    "Downstream (<=300)",
    "Distal Intergenic"
)

# Define color palette
feature_colors <- c(
    "Promoter (<=1kb)" = "#D0E0E3",
    "Promoter (1-2kb)" = "#4F81BD",
    "Promoter (2-3kb)" = "#C3D69B",
    "5' UTR" = "#76933C",
    "3' UTR" = "#FABF8F",
    "1st Exon" = "#C0504D",
    "Other Exon" = "#F9CB9C",
    "1st Intron" = "#F79646",
    "Other Intron" = "#D5A6BD",
    "Downstream (<=300)" = "#604A7B",
    "Distal Intergenic" = "#FFF2CC"
)

# Create output directory and subdirectories
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Create subdirectories for different output types
figures_dir <- file.path(output_dir, "figures")
if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
}

rds_dir <- file.path(output_dir, "rds")
if (!dir.exists(rds_dir)) {
    dir.create(rds_dir, recursive = TRUE)
}

logs_dir <- file.path(output_dir, "logs")
if (!dir.exists(logs_dir)) {
    dir.create(logs_dir, recursive = TRUE)
}

# Create subdirectories for visualization types
default_piecharts_dir <- file.path(figures_dir, "default_piecharts")
if (save_default_piechart && !dir.exists(default_piecharts_dir)) {
    dir.create(default_piecharts_dir, recursive = TRUE)
}

custom_piecharts_dir <- file.path(figures_dir, "custom_piecharts")
if (!dir.exists(custom_piecharts_dir)) {
    dir.create(custom_piecharts_dir, recursive = TRUE)
}

# Set up genome and annotation
cat("Loading TxDb annotation database...\n")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#' @title Process Peak File
#' @description Process a narrowPeak file and perform genomic annotation
#' @param file_path Path to the narrowPeak file
#' @param sample_name Name of the sample
#' @return Annotation object with genomic feature information
process_peak_file <- function(file_path, sample_name) {
    cat("Processing", sample_name, "from", basename(file_path), "\n")
    
    # Read narrowPeak format
    peak_data <- read.table(file_path, header = FALSE, 
                            col.names = c("chrom", "start", "end", "name", "score", 
                                          "strand", "signalValue", "pValue", "qValue", "peak"))
    
    # Convert to GRanges object
    gr <- GRanges(
        seqnames = peak_data$chrom,
        ranges = IRanges(start = peak_data$start, end = peak_data$end),
        name = peak_data$name,
        score = peak_data$score,
        signalValue = peak_data$signalValue,
        pValue = peak_data$pValue,
        qValue = peak_data$qValue,
        peak = peak_data$peak
    )
    
    metadata(gr)$sample <- sample_name
    
    # Annotate peaks
    cat("  Annotating peaks...\n")
    peakAnno <- annotatePeak(gr, 
                           tssRegion = c(-upstream, downstream),
                           TxDb = txdb,
                           annoDb = "org.Hs.eg.db")
    
    # Save annotation results
    anno_file <- file.path(rds_dir, paste0(sample_name, "_annotation.rds"))
    saveRDS(peakAnno, anno_file)
    
    # Save default ChIPseeker pie chart
    if (save_default_piechart) {
        pdf(file.path(default_piecharts_dir, paste0(sample_name, "_genomic_distribution_pie.pdf")), 
            width = 8, height = 8)
        plotAnnoPie(peakAnno)
        dev.off()
    }
    
    return(peakAnno)
}

# Process all files
anno_list <- list()
for (i in 1:length(input_files)) {
    anno_list[[sample_names[i]]] <- process_peak_file(input_files[i], sample_names[i])
}

#' @title Extract Annotation Statistics
#' @description Extract annotation statistics from multiple samples
#' @param anno_list List of annotation objects
#' @return Data frame with combined annotation statistics
extract_anno_stats <- function(anno_list) {
    stats_list <- lapply(names(anno_list), function(sample_name) {
        anno <- anno_list[[sample_name]]
        stats <- data.frame(
            Feature = anno@annoStat$Feature,
            Frequency = anno@annoStat$Frequency,
            Sample = sample_name
        )
        return(stats)
    })
    
    return(do.call(rbind, stats_list))
}

# Get combined annotation statistics
combined_stats <- extract_anno_stats(anno_list)

# Ensure feature order is consistent
combined_stats$Feature <- factor(combined_stats$Feature, levels = feature_order)

#' @title Create Color Palette
#' @description Create a color palette that scales with the number of samples
#' @param n Number of colors needed
#' @return Vector of colors
create_color_palette <- function(n) {
    if (n <= 8) {
        return(RColorBrewer::brewer.pal(max(3, n), "Set1")[1:n])
    } else {
        base_colors <- c(
            RColorBrewer::brewer.pal(8, "Set1"),
            RColorBrewer::brewer.pal(8, "Set2"),
            RColorBrewer::brewer.pal(8, "Set3"),
            RColorBrewer::brewer.pal(8, "Dark2"),
            RColorBrewer::brewer.pal(8, "Paired")
        )
        
        if (n <= length(base_colors)) {
            return(base_colors[1:n])
        } else {
            return(colorRampPalette(base_colors)(n))
        }
    }
}

# Create bar plot for genomic features
cat("Creating bar plot comparison of genomic features...\n")
n_samples <- length(unique(combined_stats$Sample))
sample_colors <- create_color_palette(n_samples)
names(sample_colors) <- unique(combined_stats$Sample)

p_bar <- ggplot(combined_stats, aes(x = Feature, y = Frequency, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Distribution of peaks across features",
         x = "Genomic Feature",
         y = "Frequency (%)") +
    scale_fill_manual(values = sample_colors)

# Save the bar plot
ggsave(file.path(figures_dir, "combined_genomic_distribution_bar.pdf"), 
       p_bar, width = 10, height = 6)

# Create combined TSS profile
cat("Creating combined TSS profile plot...\n")
pdf(file.path(figures_dir, "combined_TSS_distribution.pdf"), width = 10, height = 6)
plotDistToTSS(anno_list, title = "Distribution of peaks relative to TSS")
dev.off()

#' @title Create Pie Chart
#' @description Create a single pie chart for one sample
#' @param sample_data Data for a single sample
#' @param title Title for the chart
#' @return ggplot object
create_pie_chart <- function(sample_data, title) {
    pie_data <- sample_data %>%
        mutate(
            prop = Frequency / sum(Frequency)
        )
    
    feature_labels <- sapply(levels(pie_data$Feature), function(feat) {
        if (feat %in% pie_data$Feature) {
            freq <- pie_data$Frequency[pie_data$Feature == feat]
            return(paste0(feat, " (", sprintf("%.1f%%", freq), ")"))
        } else {
            return(feat)
        }
    })
    
    p <- ggplot(pie_data, aes(x = "", y = Frequency, fill = Feature)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 12),
            legend.title = element_blank(),
            legend.text = element_text(size = 8)
        ) +
        labs(title = title) +
        scale_fill_manual(values = feature_colors, 
                          labels = if(show_percentages) feature_labels else levels(pie_data$Feature),
                          breaks = levels(pie_data$Feature))
    
    return(p)
}

# Create pie charts grid
cat("Creating multi-sample pie chart grid...\n")
pie_charts <- list()

for (sample_name in names(anno_list)) {
    sample_data <- subset(combined_stats, Sample == sample_name)
    pie_charts[[sample_name]] <- create_pie_chart(sample_data, sample_name)
}

n_samples <- length(anno_list)
n_cols <- min(3, n_samples)
n_rows <- ceiling(n_samples / n_cols)

tryCatch({
    pdf(file.path(figures_dir, "multi_sample_pie_charts_grid_with_percentages.pdf"), 
        width = 3*n_cols + 4, height = 4*n_rows + 4)
    do.call(gridExtra::grid.arrange, c(pie_charts, ncol = n_cols, 
                                       top = "Genomic Distribution Across Samples"))
    dev.off()
}, error = function(e) {
    cat("Note: Could not create grid layout with grid.arrange:", e$message, "\n")
})

# Create individual pie chart files
tryCatch({
    for (sample_name in names(anno_list)) {
        sample_data <- subset(combined_stats, Sample == sample_name)
        p <- create_pie_chart(sample_data, sample_name)
        
        ggsave(file.path(custom_piecharts_dir, paste0(sample_name, "_custom_pie_chart.pdf")), 
               p, width = 10, height = 8)
    }
    
    cat("Created individual pie chart files with percentage labels\n")
}, error = function(e) {
    cat("Note: Could not create individual pie charts:", e$message, "\n")
})

#' @title Draw Arc
#' @description Helper function for drawing arcs in donut charts
#' @param x X center coordinate
#' @param y Y center coordinate
#' @param inner.r Inner radius
#' @param outer.r Outer radius
#' @param start.angle Start angle in radians
#' @param end.angle End angle in radians
#' @param col Fill color
#' @param border Border color
draw.arc <- function(x, y, inner.r, outer.r, start.angle, end.angle, col, border = NA) {
    start.angle.deg <- start.angle * 180 / pi
    end.angle.deg <- end.angle * 180 / pi
    
    segments <- 100
    angles <- seq(start.angle.deg, end.angle.deg, length.out = segments)
    angles.rad <- angles * pi / 180
    
    x.inner <- x + inner.r * cos(angles.rad)
    y.inner <- y + inner.r * sin(angles.rad)
    x.outer <- x + outer.r * cos(angles.rad)
    y.outer <- y + outer.r * sin(angles.rad)
    
    polygon(c(x.outer, rev(x.inner)), c(y.outer, rev(y.inner)), 
            col = col, border = border)
}

# Create stacked donut visualization
cat("Creating stacked donut visualization...\n")
tryCatch({
    # Prepare data for donut chart
    donut_data <- combined_stats
    samples <- unique(donut_data$Sample)
    
    # Dynamically calculate chart dimensions
    n_samples <- length(samples)
    
    if (n_samples <= 10) {
        base_r <- 0.8
        ring_width <- 0.18
        pdf_width <- 10
        pdf_height <- 10
        label_cex <- 0.8
        margin_factor <- 1.5
    } else if (n_samples <= 30) {
        base_r <- 0.5
        ring_width <- 0.1
        pdf_width <- 12
        pdf_height <- 12
        label_cex <- 0.6
        margin_factor <- 2.0
    } else {
        base_r <- 0.3
        ring_width <- max(0.02, 2 / n_samples)
        pdf_width <- 15
        pdf_height <- 15
        label_cex <- max(0.2, 5 / n_samples)
        margin_factor <- 2.5
    }
    
    # Create position columns for each sample
    for (i in seq_along(samples)) {
        sample_name <- samples[i]
        inner_r <- base_r + (i-1) * ring_width
        outer_r <- inner_r + ring_width
        
        donut_data$inner_r[donut_data$Sample == sample_name] <- inner_r
        donut_data$outer_r[donut_data$Sample == sample_name] <- outer_r
    }
    
    max_radius <- base_r + n_samples * ring_width
    plot_margin <- max_radius * margin_factor
    
    # Prepare angles
    donut_data <- donut_data %>%
        group_by(Sample) %>%
        arrange(Feature) %>%
        mutate(
            end_angle = cumsum(Frequency) * 2 * pi / 100,
            start_angle = lag(end_angle, default = 0)
        )
    
    # Create donut chart
    pdf(file.path(figures_dir, "concentric_donut_chart.pdf"), width = pdf_width, height = pdf_height)
    
    par(mar = c(max(1, n_samples/5), 1, 3, 1))
    plot(NULL, xlim = c(-plot_margin, plot_margin), ylim = c(-plot_margin, plot_margin), axes = FALSE, 
         xlab = "", ylab = "", main = "Genomic Feature Distribution\nAcross Samples")
    
    feature_vector <- unique(as.character(donut_data$Feature))
    segment_colors <- feature_colors[match(feature_vector, names(feature_colors))]
    
    # Draw segments
    for (i in 1:nrow(donut_data)) {
        row <- donut_data[i,]
        
        if (row$Frequency > 0) {
            feature_color <- feature_colors[as.character(row$Feature)]
            
            draw.arc(0, 0, row$inner_r, row$outer_r, row$start_angle, row$end_angle, 
                     col = feature_color, border = "white")
        }
    }
    
    # Add labels based on sample count
    if (n_samples <= 20) {
        for (i in seq_along(samples)) {
            sample_name <- samples[i]
            r <- base_r + (i-0.5) * ring_width
            text(0, 0 - r, sample_name, srt = 0, adj = c(0.5, 0.5), cex = label_cex)
        }
    } else {
        # Create a separate legend for many samples
        legend_x <- -plot_margin * 0.8
        legend_y_start <- plot_margin * 0.6
        legend_y_step <- (2 * plot_margin * 0.8) / n_samples
        
        # Left side labels
        for (i in 1:ceiling(n_samples/2)) {
            sample_name <- samples[i]
            r <- base_r + (i-0.5) * ring_width
            legend_y <- legend_y_start - (i-1) * legend_y_step
            
            lines(c(-r*0.9, legend_x + 0.1), c(0, legend_y), col = "gray70", lwd = 0.5)
            text(legend_x, legend_y, sample_name, adj = c(1, 0.5), cex = label_cex)
        }
        
        # Right side labels
        legend_x <- plot_margin * 0.8
        for (i in (ceiling(n_samples/2) + 1):n_samples) {
            sample_name <- samples[i]
            r <- base_r + (i-0.5) * ring_width
            legend_y <- legend_y_start - (i-ceiling(n_samples/2)-1) * legend_y_step
            
            lines(c(r*0.9, legend_x - 0.1), c(0, legend_y), col = "gray70", lwd = 0.5)
            text(legend_x, legend_y, sample_name, adj = c(0, 0.5), cex = label_cex)
        }
    }
    
    # Add feature legend
    if (n_samples <= 10) {
        legend("topright", legend = names(feature_colors), 
               fill = feature_colors, cex = 0.7, bty = "n")
    } else {
        legend("top", legend = names(feature_colors), 
               fill = feature_colors, cex = 0.6, bty = "n", 
               ncol = min(4, ceiling(length(feature_colors)/3)))
    }
    
    dev.off()
    
    # Create multi-page version for large sample counts
    if (n_samples > 30) {
        cat("Creating multi-page donut charts for large sample count...\n")
        samples_per_page <- 20
        num_pages <- ceiling(n_samples / samples_per_page)
        
        pdf(file.path(figures_dir, "concentric_donut_chart_multipage.pdf"), width = 10, height = 10)
        
        for (page in 1:num_pages) {
            start_idx <- (page-1) * samples_per_page + 1
            end_idx <- min(page * samples_per_page, n_samples)
            page_samples <- samples[start_idx:end_idx]
            
            page_data <- donut_data %>% filter(Sample %in% page_samples)
            
            par(mar = c(3, 3, 4, 3))
            plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), axes = FALSE, 
                xlab = "", ylab = "", 
                main = paste0("Genomic Feature Distribution\nSample Page ", page, " of ", num_pages))
            
            for (i in 1:nrow(page_data)) {
                row <- page_data[i,]
                
                if (row$Frequency > 0) {
                    feature_color <- feature_colors[as.character(row$Feature)]
                    
                    draw.arc(0, 0, row$inner_r, row$outer_r, row$start_angle, row$end_angle, 
                            col = feature_color, border = "white")
                }
            }
            
            for (i in seq_along(page_samples)) {
                sample_name <- page_samples[i]
                sample_idx <- which(samples == sample_name)
                r <- base_r + (sample_idx-0.5) * ring_width
                text(0, 0 - r, sample_name, srt = 0, adj = c(0.5, 0.5), cex = 0.8)
            }
            
            legend("topright", legend = names(feature_colors), 
                  fill = feature_colors, cex = 0.7, bty = "n")
        }
        
        dev.off()
    }
    
}, error = function(e) {
    cat("Note: Could not create concentric donut chart:", e$message, "\n")
})

# Create faceted pie chart visualization
cat("Creating faceted pie chart visualization...\n")
tryCatch({
    n_samples <- length(unique(combined_stats$Sample))
    
    if (n_samples <= 4) {
        n_rows <- 1
        n_cols <- n_samples
        strip_size <- 12
        pdf_width <- 10
        pdf_height <- 6
    } else if (n_samples <= 12) {
        n_rows <- 2
        n_cols <- ceiling(n_samples / n_rows)
        strip_size <- 10
        pdf_width <- min(16, 3 * n_cols)
        pdf_height <- 8
    } else if (n_samples <= 24) {
        n_rows <- 3
        n_cols <- ceiling(n_samples / n_rows)
        strip_size <- 8
        pdf_width <- min(20, 2.5 * n_cols)
        pdf_height <- 10
    } else {
        n_rows <- 4
        n_cols <- 6
        samples_per_page <- n_rows * n_cols
        strip_size <- 7
        pdf_width <- 16
        pdf_height <- 12
    }
    
    pie_theme <- theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            strip.text = element_text(size = strip_size),
            panel.spacing = unit(0.5, "lines"),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)
        )
    
    if (n_samples <= 24) {
        p_facet <- ggplot(combined_stats, aes(x = "", y = Frequency, fill = Feature)) +
            geom_bar(stat = "identity", width = 1, color = "white") +
            coord_polar("y", start = 0) +
            facet_wrap(~ Sample, nrow = n_rows) +
            pie_theme +
            labs(title = "Genomic Distribution Across Samples") +
            scale_fill_manual(values = feature_colors)
        
        ggsave(file.path(figures_dir, "faceted_pie_charts.pdf"), 
               p_facet, width = pdf_width, height = pdf_height)
    } else {
        all_samples <- unique(combined_stats$Sample)
        num_pages <- ceiling(n_samples / samples_per_page)
        
        pdf(file.path(figures_dir, "faceted_pie_charts.pdf"), 
            width = pdf_width, height = pdf_height)
        
        for (page in 1:num_pages) {
            start_idx <- (page-1) * samples_per_page + 1
            end_idx <- min(page * samples_per_page, n_samples)
            page_samples <- all_samples[start_idx:end_idx]
            
            page_data <- combined_stats %>% filter(Sample %in% page_samples)
            
            p_page <- ggplot(page_data, aes(x = "", y = Frequency, fill = Feature)) +
                geom_bar(stat = "identity", width = 1, color = "white") +
                coord_polar("y", start = 0) +
                facet_wrap(~ Sample, nrow = n_rows) +
                pie_theme +
                labs(title = paste0("Genomic Distribution - Page ", page, " of ", num_pages)) +
                scale_fill_manual(values = feature_colors)
            
            print(p_page)
        }
        
        dev.off()
        
        # Create compact overview
        if (n_samples > 50) {
            cat("Creating compact overview of all samples...\n")
            overview_rows <- ceiling(sqrt(n_samples))
            overview_cols <- ceiling(n_samples / overview_rows)
            
            p_overview <- ggplot(combined_stats, aes(x = "", y = Frequency, fill = Feature)) +
                geom_bar(stat = "identity", width = 1, color = "white") +
                coord_polar("y", start = 0) +
                facet_wrap(~ Sample, nrow = overview_rows) +
                theme_minimal() +
                theme(
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid = element_blank(),
                    strip.text = element_text(size = 6),
                    panel.spacing = unit(0.2, "lines"),
                    legend.position = "bottom",
                    legend.title = element_blank(),
                    legend.text = element_text(size = 7),
                    plot.title = element_text(hjust = 0.5, size = 10)
                ) +
                labs(title = "Genomic Distribution - Overview of All Samples") +
                scale_fill_manual(values = feature_colors)
            
            ggsave(file.path(figures_dir, "faceted_pie_charts_overview.pdf"), 
               p_overview, width = 16, height = 16, limitsize = FALSE)
        }
    }
    
}, error = function(e) {
    cat("Note: Could not create faceted pie chart:", e$message, "\n")
})

# Calculate script execution time
script_end_time <- Sys.time()
execution_time <- script_end_time - script_start_time

# Create a summary file with analysis information
summary_file <- file.path(logs_dir, "analysis_summary.txt")
cat("Creating summary file:", summary_file, "\n")

# Write summary information
sink(summary_file)
cat("# Genomic Annotation Analysis Summary\n\n")
cat("## Analysis Parameters\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Upstream region:", upstream, "bp\n")
cat("Downstream region:", downstream, "bp\n")
cat("Sample regex pattern:", sample_regex, "\n")
cat("Samples processed:", length(anno_list), "\n\n")

cat("## Processed Files\n")
for (i in 1:length(input_files)) {
  cat(sprintf("%d. %s -> %s\n", i, basename(input_files[i]), sample_names[i]))
}

cat("\n## Execution Time\n")
cat("Start time:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", format(execution_time, digits=2), "\n")
sink()

cat("All processing complete!\n")
cat("Results organized in the following directories:\n")
cat("  - RDS files saved to:      ", rds_dir, "\n")
cat("  - Figures saved to:        ", figures_dir, "\n")
cat("  - Log files saved to:      ", logs_dir, "\n")
cat("\nTotal execution time:", format(execution_time, digits=2), "\n")
cat("\nUse --help for more options and information.\n")
