shhh <- suppressPackageStartupMessages

shhh({
	library(ggplot2)
	library(data.table)
	library(dplyr)
	library(cowplot)
	library(tidyverse)
	library(stringr)
    library(optparse)
})

# Helper function to parse GTF attributes
extract_attributes <- function(gtf_attributes, att_of_interest) {
    att <- unlist(strsplit(gtf_attributes, " "))
    if (att_of_interest %in% att) {
        return(gsub("\"|;", "", att[which(att %in% att_of_interest) + 1]))
    } else {
        return(NA)
    }
}

make_gene_reference <- function(gtf_path, out_file) {
    # goal: chr	start	end	name	score	strand	Ensembl_ID	gene_type

    # extract protein coding genes from gtf
    gtf <- fread(gtf_path, header = FALSE, sep = "\t") %>% 
        setNames(c("chr","source","type","start","end","score","strand","phase", "attributes"))
    gtf$gene_type <- unlist(lapply(gtf$attributes, extract_attributes, "gene_type"))
    genes <- filter(gtf, type == "gene", gene_type == "protein_coding")
    message("# protein-coding genes: ", nrow(genes))

    # add gene symbol Ensembl ID
    genes$name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))
    genes$Ensembl_ID <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
    genes$score <- 0

    # remove genes with more than one gene symbol per Ensembl ID (no decimals), not on main chromosomes
    genes <- genes %>% select(chr, start, end, name, score, strand, Ensembl_ID, gene_type) %>% 
        mutate(Ensembl_ID = sub("\\.\\d+$", "", Ensembl_ID)) %>%
        filter(grepl("^chr", chr)) %>% # remove non "chr" chromosomes
        group_by(name) %>% # remove genes with more than one name per ensembl id 
        mutate(n_gene_names = n()) %>%
        filter(n_gene_names == 1) %>%
        select(-n_gene_names)
    fwrite(genes, out_file, sep = "\t", col.names = FALSE, quote = FALSE)

    # write bed with header
    header_line <- paste0("#", paste(names(genes), collapse = "\t"))
    writeLines(c(header_line, readLines(out_file)), out_file)

    message("# genes in final reference: ", nrow(genes))
}


find_one_tss <- function(this_Ensembl_ID, annot_tx, reference_regions = NULL) {
    tx_this <- annot_tx %>% filter(Ensembl_ID == this_Ensembl_ID)

    if (nrow(tx_this) == 0) {
        return(list(0, "NO transcripts!"))
    }

    if (length(unique(tx_this$tss)) == 1) {
        return(list(tx_this$tss[1], "only TSS"))
    }

    # Track whether we filtered by reference overlap
    used_ref_filter <- FALSE
    support_prefix <- ""

    # Check reference overlap if reference provided
    if (!is.null(reference_regions) && nrow(reference_regions) > 0) {
        gene_name <- unique(tx_this$gene_name)[1]
        chr <- tx_this$chr[1]

        # Filter reference to this gene by name
        ref_this_gene <- reference_regions[reference_regions$name == gene_name &
                                           reference_regions$chr == chr, ]

        if (nrow(ref_this_gene) > 0) {
            unique_tss <- unique(tx_this$tss)

            # Check which TSS candidates (±250bp) overlap reference regions
            overlaps <- sapply(unique_tss, function(tss_candidate) {
                tss_start <- tss_candidate - 250
                tss_end <- tss_candidate + 250
                any(ref_this_gene$start < tss_end & ref_this_gene$end > tss_start)
            })

            overlapping_tss <- unique_tss[overlaps]

            if (length(overlapping_tss) == 1) {
                return(list(overlapping_tss, "reference overlap"))
            }

            # If multiple overlap, filter to only those transcripts and continue
            if (length(overlapping_tss) > 1) {
                tx_this <- tx_this %>% filter(tss %in% overlapping_tss)
                used_ref_filter <- TRUE
                support_prefix <- "reference overlap + "
            }
            # If none overlap, continue with all transcripts (no filtering)
        }
    }

    # try level 1 and level 2 confidence
    for (i in 1:2){
        this_level <- filter(tx_this, annot_level == i)

        if (length(unique(this_level$tss)) == 1) {
            return(list(this_level$tss[1], paste0(support_prefix, "only level ", i, " TSS")))
        }

        if (length(unique(this_level$tss)) > 1) {
            this_level <- arrange(this_level, desc(tx_length))
            return(list(this_level$tss[1], paste0(support_prefix, "level ", i, " TSS from longest transcript")))
        }
    }

    # no level 1 or level 2, so turn to "transcript support level" from Ensembl
    for (i in 1:3){
        this_level <- filter(tx_this, tsl == i)

        if (length(unique(this_level$tss)) == 1) {
            return(list(this_level$tss[1], paste0(support_prefix, "only TSL ", i, " TSS")))
        }

        if (length(unique(this_level$tss)) > 1) {
            this_level <- arrange(this_level, desc(tx_length))
            return(list(this_level$tss[1], paste0(support_prefix, "TSL ", i, " TSS from longest transcript")))
        }
    }

    # resort to longest transcript
    long_to_short <- arrange(tx_this, desc(tx_length))
    suffix <- ifelse(used_ref_filter, "longest transcript from reference-overlapping TSS",
                     "no well-supported transcripts, returning TSS from longest transcript")
    return(list(long_to_short$tss[1], paste0(support_prefix, suffix)))
}

generate_qc_plots <- function(tss_res, genes, comparison_bed_path, qc_plot_path) {
    # Nature color palette
    nature_blue <- c("#c5e5fb", "#9bcae9", "#5496ce", "#006eae", "#00488d", "#002359")
    nature_green <- c("#d7e5c5", "#a1ca78", "#5eb342", "#429130", "#1c6e2b", "#0e3716")
    nature_teal <- c("#cae5ee", "#96ced3", "#49bcbc", "#0096a0", "#006479", "#003648")
    nature_yellow <- c("#ffeec1", "#f6dc87", "#e9c54e", "#ca9b23", "#9b740a", "#69540b")
    nature_red <- c("#f6ceca", "#e9a0a5", "#dc6464", "#c5373d", "#9b241c", "#730c0d")

    # Custom theme (no gridlines, black axes/text)
    theme_nature <- function() {
        theme_classic() +
        theme(
            text = element_text(color = "black"),
            axis.text = element_text(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.line = element_line(color = "black"),
            panel.grid = element_blank(),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9)
        )
    }

    # Plot 1: TSS selection path distribution
    tss_smry <- tss_res %>%
        group_by(support) %>%
        summarize(n_genes = n(), .groups = "drop") %>%
        mutate(support = factor(support, levels = support[order(n_genes)]))

    p1 <- ggplot(tss_smry, aes(x = support, y = n_genes)) +
        geom_bar(stat = "identity", fill = nature_blue[4]) +
        coord_flip() +
        labs(x = NULL, y = "Number of genes",
             title = "TSS selection method distribution") +
        theme_nature()

    plot_list <- list(p1)

    # Plot 2 and 3: Comparison with reference (if provided)
    if (!is.null(comparison_bed_path) && file.exists(comparison_bed_path)) {
        comparison <- fread(comparison_bed_path,
            col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"))

        # Plot 2: Gene overlap
        new_genes <- unique(genes$Ensembl_ID)
        ref_genes <- unique(comparison$Ensembl_ID)

        overlap_data <- data.frame(
            Category = c("In both", "New only", "Reference only"),
            Count = c(
                length(intersect(new_genes, ref_genes)),
                length(setdiff(new_genes, ref_genes)),
                length(setdiff(ref_genes, new_genes))
            )
        )
        overlap_data$Category <- factor(overlap_data$Category,
            levels = c("In both", "New only", "Reference only"))

        p2 <- ggplot(overlap_data, aes(x = Category, y = Count, fill = Category)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c(
                "In both" = nature_blue[4],
                "New only" = nature_green[4],
                "Reference only" = nature_red[4]
            ), guide = "none") +
            labs(x = NULL, y = "Number of genes",
                 title = "Gene overlap with reference") +
            theme_nature()

        plot_list <- c(plot_list, list(p2))

        # Plot 3: TSS distance distribution
        comparison <- comparison %>%
            mutate(ref_tss = (start + end) / 2)

        tss_comparison <- tss_res %>%
            inner_join(comparison %>% select(Ensembl_ID, ref_tss),
                       by = c("gene_Ensembl_ID" = "Ensembl_ID")) %>%
            mutate(distance = abs(tss - ref_tss))

        if (nrow(tss_comparison) > 0) {
            tss_comparison <- tss_comparison %>%
                mutate(distance_bin = cut(distance,
                    breaks = c(-1, 0, 100, 500, 1000, 5000, Inf),
                    labels = c("Identical (0)", "1-100", "101-500", "501-1000", "1001-5000", ">5000")))

            n_identical <- sum(tss_comparison$distance == 0)
            n_within_500 <- sum(tss_comparison$distance <= 500)
            n_total <- nrow(tss_comparison)
            median_dist <- median(tss_comparison$distance)

            p3 <- ggplot(tss_comparison, aes(x = distance_bin)) +
                geom_bar(fill = nature_teal[4]) +
                labs(x = "Distance to reference TSS (bp)", y = "Number of genes",
                     title = "TSS distance from reference",
                     subtitle = sprintf("Median: %.0f bp | Identical: %.1f%% | Within 500bp: %.1f%%",
                         median_dist, 100*n_identical/n_total, 100*n_within_500/n_total)) +
                theme_nature() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

            plot_list <- c(plot_list, list(p3))
        }
    }

    # Combine and save plots
    combined <- plot_grid(plotlist = plot_list, ncol = 1, align = "v")

    ggsave(qc_plot_path, combined, width = 8, height = 4 * length(plot_list),
           units = "in", device = "pdf")
    message("QC plots saved to: ", qc_plot_path)
}

make_tss_reference <- function(gtf_path, gene_ref_path, out_file_interm, out_file,
                               reference_bed_path = NULL,
                               comparison_bed_path = NULL,
                               qc_plot_path = NULL) {
    genes <- fread(gene_ref_path, col.names = c("chr", "start", "end", "name", "score", "strand", "Ensembl_ID", "gene_type"))

    # Load reference regions for TSS overlap check if provided
    reference_regions <- NULL
    if (!is.null(reference_bed_path) && file.exists(reference_bed_path)) {
        reference_regions <- fread(reference_bed_path,
            col.names = c("chr", "start", "end", "name", "score", "strand", "transcript_id", "source"))
        message("Loaded reference BED with ", nrow(reference_regions), " regions")
    }

    gtf <- fread(gtf_path, header = FALSE, sep = "\t") %>%
        setNames(c("chr","source","type","start","end","score","strand","phase", "attributes"))
    gtf$gene_type <- unlist(lapply(gtf$attributes, extract_attributes, "gene_type"))

    # filter gtf to protein-coding transcripts
    # compute tss and transcript length
    tx <- filter(gtf, type == "transcript", gene_type == "protein_coding") %>%
        mutate(tss = ifelse(strand == "+", start, end),
            tx_length = abs(start - end)) %>%
        filter(!is.na(tss))

    # add gene symbol, annotation "level", transcript support level, Ensembl ID (no decimals)
    tx$gene_name <- unlist(lapply(tx$attributes, extract_attributes, "gene_name"))
    tx$annot_level <- as.numeric(unlist(lapply(tx$attributes, extract_attributes, "level")))
    tx$tsl <- as.numeric(unlist(lapply(tx$attributes, extract_attributes, "transcript_support_level")))
    tx$Ensembl_ID <- unlist(lapply(tx$attributes, extract_attributes, "gene_id"))
    tx <- tx %>% mutate(Ensembl_ID = sub("\\.\\d+$", "", Ensembl_ID))

    print(unique(tx$tsl))
    print(unique(tx$annot_level))

    # define a single TSS for each gene in the gene reference file
    tss_res <- lapply(genes$Ensembl_ID, find_one_tss, tx, reference_regions) %>%
        rbindlist() %>%
        setNames(c("tss", "support"))
    message("# rows in TSS res: ", nrow(tss_res))
    message("# non-NA rows in TSS res: ", length(tss_res$tss[!is.na(tss_res$tss)]))
    tss_res$gene_Ensembl_ID <- genes$Ensembl_ID
    fwrite(tss_res, out_file_interm, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

    tss_smry <- group_by(tss_res, support) %>% summarize(n_entries = n())
    print(tss_smry)

    # define promoters as TSS +/- 250 bp
    tss_key <- setNames(tss_res$tss, tss_res$gene_Ensembl_ID)
    genes_out <- genes %>% mutate(tss = tss_key[Ensembl_ID],
        start = tss - 250, end = tss + 250) %>%
    select(chr, start, end, name, score, strand, Ensembl_ID, gene_type) %>%
    arrange(chr, start)

    fwrite(genes_out, out_file, sep = "\t", col.names = FALSE, quote = FALSE)

    header_line <- paste0("#", paste(names(genes_out), collapse = "\t"))
    writeLines(c(header_line, readLines(out_file)), out_file)

    # Generate QC plots if output path provided
    if (!is.null(qc_plot_path)) {
        generate_qc_plots(tss_res, genes, comparison_bed_path, qc_plot_path)
    }
}


##################
# Command line argument parsing
option_list <- list(
    make_option(c("-g", "--gtf"), type = "character", default = NULL,
                help = "Input GTF file path (GENCODE format)", metavar = "FILE"),
    make_option(c("-o", "--output-prefix"), type = "character", default = NULL,
                help = "Output file prefix (will generate .genes.bed, .TSS_choice_support.tsv, .TSS500bp.bed, .TSS_QC.pdf)", metavar = "PREFIX"),
    make_option(c("-r", "--reference"), type = "character", default = NULL,
                help = "Optional: Reference BED file for TSS overlap filtering", metavar = "FILE"),
    make_option(c("-c", "--comparison"), type = "character", default = NULL,
                help = "Optional: Comparison BED file for QC plots", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "Generate gene and TSS reference files from GENCODE GTF")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$gtf)) {
    print_help(opt_parser)
    stop("GTF file must be specified with -g/--gtf", call. = FALSE)
}

if (is.null(opt$`output-prefix`)) {
    print_help(opt_parser)
    stop("Output prefix must be specified with -o/--output-prefix", call. = FALSE)
}

# Construct output file paths from prefix
out_genes <- paste0(opt$`output-prefix`, ".genes.bed")
out_tss_interm <- paste0(opt$`output-prefix`, ".TSS_choice_support.tsv")
out_tss <- paste0(opt$`output-prefix`, ".TSS500bp.bed")
out_qc <- paste0(opt$`output-prefix`, ".TSS_QC.pdf")

# Main execution
make_gene_reference(opt$gtf, out_genes)
make_tss_reference(opt$gtf, out_genes, out_tss_interm, out_tss,
                   reference_bed_path = opt$reference,
                   comparison_bed_path = opt$comparison,
                   qc_plot_path = out_qc)