shhh <- suppressPackageStartupMessages

shhh({
    library(data.table)
    library(dplyr)
    library(optparse)
})

convert_refgene_to_tss <- function(refgene_path, out_file) {
    # Read RefGene file (UCSC format)
    # Columns: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
    #          exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames
    refgene <- fread(refgene_path, header = TRUE, sep = "\t")
    colnames(refgene) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd",
                           "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                           "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

    message("Read ", nrow(refgene), " transcripts from RefGene file")

    # Calculate TSS based on strand
    # RefGene uses 0-based coordinates for txStart, 1-based for txEnd
    # For + strand: TSS is at txStart (convert to 1-based by adding 1 for consistency)
    # For - strand: TSS is at txEnd (already 1-based)
    refgene <- refgene %>%
        mutate(
            tss = ifelse(strand == "+", txStart + 1, txEnd),
            start = tss - 250,
            end = tss + 250,
            gene_name = name2,
            transcript_id = name,
            source = "RefSeq"
        ) %>%
        filter(grepl("^chr[0-9XY]+$", chrom))

    message("After filtering to main chromosomes: ", nrow(refgene), " transcripts")
    message("Unique gene symbols: ", length(unique(refgene$gene_name)))

    # Format output as BED8
    out <- refgene %>%
        select(chr = chrom, start, end, name = gene_name,
               score, strand, transcript_id, source) %>%
        mutate(score = 0) %>%
        arrange(chr, start)

    # Write output without header first
    fwrite(out, out_file, sep = "\t", col.names = FALSE, quote = FALSE)

    # Add header line
    header_line <- paste0("#", paste(names(out), collapse = "\t"))
    writeLines(c(header_line, readLines(out_file)), out_file)

    message("RefGene TSS file created: ", out_file)
    message("Total entries: ", nrow(out))

    return(invisible(out))
}

##################
# Command line argument parsing
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input RefGene file path (UCSC format)", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file path for TSS BED file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "Convert RefGene transcript file to TSS ±250bp BED format")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Input file must be specified with -i/--input", call. = FALSE)
}

if (is.null(opt$output)) {
    print_help(opt_parser)
    stop("Output file must be specified with -o/--output", call. = FALSE)
}

# Main execution
convert_refgene_to_tss(opt$input, opt$output)