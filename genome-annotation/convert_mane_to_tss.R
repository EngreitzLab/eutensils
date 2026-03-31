shhh <- suppressPackageStartupMessages

shhh({
    library(data.table)
    library(dplyr)
    library(optparse)
})

convert_mane_to_tss <- function(mane_path, out_file) {
    # Read MANE file (BED12+ format from UCSC)
    # Key columns: chrom, chromStart, chromEnd, name, strand, geneName2 (gene symbol)
    mane <- fread(mane_path, header = TRUE, sep = "\t")

    message("Read ", nrow(mane), " transcripts from MANE file")

    # Calculate TSS based on strand
    # MANE uses 0-based coordinates for chromStart, 1-based for chromEnd
    # For + strand: TSS is at chromStart (convert to 1-based by adding 1)
    # For - strand: TSS is at chromEnd (already 1-based)
    mane <- mane %>%
        filter(grepl("^chr[0-9XY]+$", `#chrom`)) %>%
        mutate(
            tss = ifelse(strand == "+", chromStart + 1, chromEnd),
            start = tss - 250,
            end = tss + 250,
            gene_name = geneName2,
            transcript_id = name,
            source = "MANE",
            chr = `#chrom`
        )

    message("After filtering to main chromosomes: ", nrow(mane), " transcripts")
    message("Unique gene symbols: ", length(unique(mane$gene_name)))

    # Format output as BED8
    out <- mane %>%
        select(chr, start, end, name = gene_name,
               score, strand, transcript_id, source) %>%
        mutate(score = 0) %>%
        arrange(chr, start)

    # Write output without header first
    fwrite(out, out_file, sep = "\t", col.names = FALSE, quote = FALSE)

    # Add header line
    header_line <- paste0("#", paste(names(out), collapse = "\t"))
    writeLines(c(header_line, readLines(out_file)), out_file)

    message("MANE TSS file created: ", out_file)
    message("Total entries: ", nrow(out))

    return(invisible(out))
}

##################
# Command line argument parsing
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input MANE file path (BED12+ format from UCSC)", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file path for TSS BED file", metavar = "FILE")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "Convert MANE transcript file to TSS ±250bp BED format")
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
convert_mane_to_tss(opt$input, opt$output)