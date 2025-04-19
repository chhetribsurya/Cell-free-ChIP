# 			Genomic Annotation Tool

## Overview

This doc explains the usage of Genomic Annotation Tool for analyzing cf-ChIP-seq narrowPeak files. The tool performs genomic feature annotation and creates visualizations to help understand the distribution of peaks across genomic features.

## Key Features and Outputs

The Genomic Annotation Tool produces a suite of visualizations and data outputs:

- **Genomic Feature Distribution**: Multiple visualization formats (bar charts, pie charts, donut charts) showing how peaks are distributed across genomic elements like promoters, exons, introns, and intergenic regions
- **TSS Proximity Analysis**: Plots showing peak distribution frequency relative to transcription start sites and promoter regions
- **Multi-sample Comparisons**: Specialized visualizations for comparing many samples simultaneously, with automatic scaling for large datasets
- **Detailed Annotation Data**: RDS files containing complete annotation information that can be used for further analysis
- **Analysis Logs**: Summary reports of processing parameters and execution details

These outputs provide both figures and detailed data for secondary analyses.

## Requirements

Before running the script, ensure the following prerequisites:

- R (version 3.6.0 or higher)
- Required R packages:
  - Bioconductor packages: GenomicRanges, ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db
  - CRAN packages: ggplot2, reshape2, gridExtra, grid, dplyr, scales, optparse

You can install the required packages using:

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GenomicRanges", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))

# Install CRAN packages
install.packages(c("ggplot2", "reshape2", "gridExtra", "grid", "dplyr", "scales", "optparse"))
```

Alternatively, you can use Conda/Mamba for installation:

```bash
conda create -n genomic_annotation r-base r-essentials bioconductor-genomicranges bioconductor-chipseeker bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-org.hs.eg.db r-ggplot2 r-reshape2 r-gridextra r-dplyr r-scales r-optparse
conda activate genomic_annotation
```

## Basic Usage

The script can be run directly from the command line:

```bash
Rscript genomic_annotation.R --input_dir /path/to/narrowPeaks/
```

This minimal command will:
1. Process all `.narrowPeak` files in the specified directory
2. Extract sample names using the default pattern
3. Create annotations and visualizations in `./annotation_results/` directory

## Complete Example

Here's a complete example with all parameters explicitly specified:

```bash
Rscript genomic_annotation.R \
  --input_dir /path/to/narrowPeaks/ \
  --output_dir ./annotation_results \
  --pattern "*.narrowPeak" \
  --upstream 3000 \
  --downstream 3000 \
  --sample_regex "^(.*?)_peaks\.narrowPeak$" \
  --show_percentages TRUE \
  --save_default_piechart TRUE
```

This command:

- Processes all narrowPeak files in the specified input directory
- Saves results to `./annotation_results/`
- Uses 3000bp upstream and downstream for TSS region
- Extracts sample names using the specified regex
- Shows percentages in custom built pie chart legends
- Saves default chip-seeker pie charts

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input_dir` | Directory containing narrowPeak files | (Required) |
| `-o, --output_dir` | Output directory for results | "./annotation_results" |
| `-p, --pattern` | File pattern to match in input directory | "*.narrowPeak" |
| `-u, --upstream` | Upstream distance from TSS (bp) | 3000 |
| `-d, --downstream` | Downstream distance from TSS (bp) | 3000 |
| `--sample_regex` | Regex to extract sample name from filename | "^(.*?)_peaks\\.narrowPeak$" |
| `--show_percentages` | Show percentages in pie charts | TRUE |
| `--save_default_piechart` | Save default ChIPseeker pie charts | TRUE |
| `--donut_threshold` | Number of samples at which to switch from donut to bar plot | 20 |
| `--force_donut` | Force donut plots for all sample counts | FALSE |
| `--force_barplot` | Force creation of bar plot regardless of other settings | FALSE |
| `--create_tss_plot` | Create TSS distribution plot | TRUE |
| `--help` | Show help message | |

## Examples

### Basic example with default parameters
```bash
Rscript genomic_annotation.R --input_dir /path/to/narrowPeaks/
```

### Custom output directory and TSS region
```bash
Rscript genomic_annotation.R --input_dir /path/to/narrowPeaks/ \
  --output_dir ./my_analysis_results \
  --upstream 5000 \
  --downstream 5000
```

### Custom sample name extraction for files. For example: "TP19-S91.MBD.PB11249_peaks.narrowPeak"
```bash
Rscript genomic_annotation.R --input_dir /path/to/narrowPeaks/ \
  --sample_regex "^(.*?)_peaks\.narrowPeak$"
```

### Process only specific files and disable percentage labels
```bash
Rscript genomic_annotation.R --input_dir /path/to/narrowPeaks/ \
  --pattern "*H3K27ac*narrowPeak" \
  --show_percentages FALSE
```

## Output Structure

The script organizes output files into the following directory structure:

```
output_dir/
├── figures/                # All visualizations
│   ├── default_piecharts/  # ChIPseeker default pie charts
│   │   ├── Sample1_genomic_distribution_pie.pdf
│   │   ├── Sample2_genomic_distribution_pie.pdf
│   │   └── ...
│   ├── custom_piecharts/   # Custom pie charts with percentages
│   │   ├── Sample1_custom_pie_chart.pdf
│   │   ├── Sample2_custom_pie_chart.pdf
│   │   └── ...
│   ├── combined_genomic_distribution_bar.pdf  # Bar chart comparing all samples
│   ├── combined_TSS_distribution.pdf          # TSS profiles for all samples
│   ├── concentric_donut_chart.pdf             # Concentric rings showing all samples
│   ├── concentric_donut_chart_multipage.pdf   # Multi-page version for many samples
├── rds/                    # Annotation data files
│   ├── Sample1_annotation.rds
│   ├── Sample2_annotation.rds
│   └── ...
└── logs/                   # Analysis logs
    └── analysis_summary.txt  # Summary of analysis parameters and timing
```

## Visualization Types

### 1. Bar Plot
- **File**: `figures/combined_genomic_distribution_bar.pdf`
- **Description**: Comparison of genomic feature distributions across all samples in a grouped bar chart.
- **Forced Bar Plot**: When `--force_barplot TRUE` is used, creates `figures/forced_genomic_distribution_barplot.pdf`

### 2. TSS Profile
- **File**: `figures/combined_TSS_distribution.pdf`
- **Description**: Peak distribution frequency relative to transcription start sites for all samples.
- **Control**: Can be enabled/disabled with `--create_tss_plot` parameter

### 3. Pie Charts
- **Files**: 
  - `figures/default_piecharts/*.pdf`: Standard ChIPseeker pie charts
  - `figures/custom_piecharts/*.pdf`: Custom pie charts with percentage labels
- **Description**: Individual pie charts showing genomic feature distribution for each sample.

### 4. Concentric Donut Chart
- **Files**: 
  - `figures/concentric_donut_chart.pdf`: Concentric visualization
  - `figures/concentric_donut_chart_multipage.pdf`: Multi-page version for many samples
- **Description**: Visualization showing all samples as concentric rings, allowing direct comparison of feature distribution across samples.

## Handling Large Sample Sets

The script is designed to handle a wide range of sample counts:

- For small sets (<10 samples): Creates standard visualizations
- For medium sets (10-30 samples): Adjusts layout and font sizes
- For large sets (30+ samples): Creates multi-page versions and overview visualizations
- For very large sets (100+ samples): Creates additional compact overview charts
- **Donut/Bar Threshold**: By default, switches from donut to bar plot at 20 samples (configurable with `--donut_threshold`)
- **Force Options**: 
  - Use `--force_donut TRUE` to always create donut charts
  - Use `--force_barplot TRUE` to always create bar plots

## Notes for Best Results

1. **Sample Name Extraction**: The default regex pattern extracts everything before "_peaks.narrowPeak". Adjust the `--sample_regex` parameter if the filenames have a different pattern.

2. **Memory Usage**: For very large datasets (100+ samples), ensure the system has sufficient memory (10GB+ recommended).

3. **Execution Time**: Processing time increases with the number of samples. For large datasets, the script displays progress information.

4. **Output Quality**: All visualizations are saved as PDF files for high-resolution output.

5. **Visualization Control**: 
   - Use `--force_barplot TRUE` to ensure bar plot creation
   - Use `--force_donut TRUE` to force donut charts even for large sample sets
   - Adjust `--donut_threshold` to control when the script switches from donut to bar plots
   - Use `--create_tss_plot FALSE` to disable TSS distribution plot creation

## Troubleshooting

If you encounter issues:

1. **Missing packages**: Ensure all required packages are installed
2. **Memory errors**: For large datasets, try increasing R memory allocation
3. **Sample name issues**: Verify that your `--sample_regex` correctly extracts sample names
4. **File access errors**: Check input file permissions and paths

## Contact

##### For questions, improvements or issues, please contact: chhetribsurya@gmail.com

This tool was developed to streamline the analysis of cf-ChIP-seq peak data to provide comprehensive genomic feature annotations and visualizations.
