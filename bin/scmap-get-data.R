#!/usr/bin/env Rscript 

### This script provides a quick way to download and standardise 10X expression data 
### directories. It also allows to extract any metadata (inferred cell types) fro reference experiment    

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))

# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-m", "--matrix-file-url"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'URL to expression matrix file in .mtx format'
    ),
    make_option(
        c("-b", "--barcodes-file-url"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'URL to barcodes (matrix column names) file in .tsv format'
    ),
    make_option(
        c("-g", "--gene-id-file-url"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'URL to gene ID (matrix rows) file in .tsv format'
    ),
    make_option(
        c("-s", "--metadata-file-url"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'URL path of reference metadata file (known cell type annotations)'
    )
)
opt = wsc_parse_args(option_list, mandatory = c("matrix_file_url", "barcodes_file_url", "gene_id_file_url"))

# download matrix
download.file(url = opt$matrix_file_url, destfile="matrix.mtx")

# download barcodes 
download.file(url = opt$barcodes_file_url, destfile="barcodes.tsv")

# download gene IDs 
download.file(url = opt$gene_id_file_url, destfile="genes.tsv")

# download metadata 
if(!is.na(opt$metadata_file_url)){
    download.file(url = opt$metadata_file_url, destfile="reference_metadata.txt")
}


