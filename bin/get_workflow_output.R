#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

### create final workflow output in standard format
option_list = list(
    make_option(
        c("-i", "--predictions-file"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the predictions file in text format'
    ),
    make_option(
        c("-o", "--workflow-output"),
        action = 'store',
        default = NA,
        type = 'character',
        help = 'Path to the final output file in text format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("predictions_file", "workflow_output"))
data = read.csv(file=opt$predictions_file)
output = data[, c('cell', 'combined_labs')]
colnames(output) = c("cell_id", "predicted_label")
write.table(output, file = opt$workflow_output, sep="\t", row.names=FALSE)