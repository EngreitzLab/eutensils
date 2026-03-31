# Overview
This directory contains the code used to derive one-promoter-per-gene TSS and gene reference files.

The following strategy is used to create the reference files:
1. Read in full GENCODE transcriptome annotations
2. Filter to “protein_coding” gene type
3. Retain genes with exactly 1 Ensembl ID (no decimals) per gene symbol, and on main chromosomes
4. Use start and end coordinates from the gene entry for gene reference file
5. Choose 1 TSS/gene for TSS reference file:
    - If gene has only one transcript, use that TSS
    - If gene has more than one transcript:
      - If additional reference annotations are provided, check for overlap of any TSS +/-250bp with the reference regions
      - If only one overlaps, use that TSS
      - If none overlap, continue to next step
      - If more than one overlap, continue to next step with just transcripts that have an overlapping TSS
    - Consider look at remaining transcripts with “level 1” and “level 2” annotation confidence (GENCODE annotation)
      - If there are level 1 transcripts, use TSS from longest option; otherwise use TSS from longest level 2 transcript
      - If there are no level 1 or 2 transcripts, look at any transcripts with Ensembl “transcript support level” (TSL) of 1–3
      - Use TSS from longest transcript with highest TSL available
      - If no level 1/2 or TSL 1–3 transcripts, use TSS from longest transcript 

## Reference files 
### GENCODE v43 + MANE annotations (default for hg38)
- GENCODE v43 GTF downloaded from [here](https://www.gencodegenes.org/human/release_43.html)
- MANE annotations downloaded from [here](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=3761611831_cxuXxnSLUAg8kHKdyeA9JteSAeyR&db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=0&hgta_regionType=genome&position=chr7%3A155%2C796%2C286-155%2C809%2C628&hgta_outputType=primaryTable&hgta_outFileName=) on March 3, 2026

- 500 bp TSS regions with header `#chr	start	end	name	score	strand	Ensembl_ID	gene_type` (use for ABC, ENCODE-rE2G, scE2G pipelines)
  - Path on Oak: `/oak/stanford/groups/engreitz/Data/hg38/genome-annotation/gencode.v43.protein_coding.TSS500bp.bed`
  - [Download](https://mitra.stanford.edu/engreitz/oak/Data/hg38/genome-annotation/gencode.v43.protein_coding.TSS500bp.bed)
- Gene body regions with header  `#chr	start	end	name	score	strand	Ensembl_ID	gene_type` (use for ABC, ENCODE-rE2G, scE2G pipelines)
  - Path on Oak: `/oak/stanford/groups/engreitz/Data/hg38/genome-annotation/gencode.v43.protein_coding.genes.bed`
  - [Download](https://mitra.stanford.edu/engreitz/oak/Data/hg38/genome-annotation/gencode.v43.protein_coding.genes.bed)
- 500 bp TSS regions bed6 format (use for CRISPR_comparison pipeline)
  - Path on Oak: `/oak/stanford/groups/engreitz/Data/hg38/genome-annotation/gencode.v43.protein_coding.TSS500bp.bed6`
  - [Download](https://mitra.stanford.edu/engreitz/oak/Data/hg38/genome-annotation/gencode.v43.protein_coding.TSS500bp.bed6)

## How to use
### Requirements
- `R`, `dplyr`, `data.table`, `optparse`, `ggplot2`, `cowplot`, `tidyverse` `stringr`

### Inputs
1. GENCODE GTF file
2. (optional) Additional transcript annotations to prioritize which TSS to choose for each gene (download from UCSC Genome Browser)

### Run
1. If using additional annotations, convert the raw UCSC download to format with header `#chr	start	end	name	score	strand	transcript_id	source`
  - `convert_mane_to_tss.R` does this for MANE annotation tables downloaded from UCSC
  - `convert_refgene_to_tss.R` does this for RefGene annotation tables downloaded from UCSC
2. Run `make_gene_reference_files.R` to generate TSS & gene reference files as well as some summary statistics & plots on TSS choice





