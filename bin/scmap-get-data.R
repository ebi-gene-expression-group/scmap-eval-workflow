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
    ), 
    make_option(
        c("-o", "--output-dir-path"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the output directory'
    )

)
opt = wsc_parse_args(option_list, mandatory = c("matrix_file_url", "barcodes_file_url", 
                                                "gene_id_file_url", "output_dir_path"))

dir.create(opt$output_dir_path)
# download matrix
d = paste0(opt$output_dir_path, "matrix.mtx")
download.file(url = opt$matrix_file_url, destfile=d)
if(!file.exists(d)) stop("Expression matrix file failed to download")

# download barcodes 
d = paste0(opt$output_dir_path, "barcodes.tsv")
download.file(url = opt$barcodes_file_url, destfile=d)
if(!file.exists(d)) stop("Barcodes file failed to download")

# download gene IDs 
d = paste0(opt$output_dir_path, "genes.tsv")
download.file(url = opt$gene_id_file_url, destfile=d)
if(!file.exists(d)) stop("Gene names file failed to download")

if(!is.na(opt$metadata_file_url)){
    # download metadata 
    download.file(url = opt$metadata_file_url, destfile="reference_metadata.txt")
    if(!file.exists("reference_metadata.txt")) stop("Metadata file failed to download")
}
