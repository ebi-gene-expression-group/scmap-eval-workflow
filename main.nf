#!/usr/bin/env nextflow 

// data-extracting steps
QUERY_MAT_URL = Channel.from(params.query_matrix_url)
QUERY_BARCODES_URL = Channel.from(params.query_barcodes_url)
QUERY_GENES_URL = Channel.from(params.query_genes_url)

process get_query_data{
    publishDir "${baseDir}/data", mode: 'copy'

    conda 'envs/dropletutils.yaml'
    input:
        val matrix_url from QUERY_MAT_URL
        val barcodes_url from QUERY_BARCODES_URL
        val genes_url from QUERY_GENES_URL
       
    output:
        file("${params.query_raw_data}") into QUERY_DIR

    """
    scmap-get-data.R\
            --matrix-file-url ${matrix_url}\
            --barcodes-file-url ${barcodes_url}\
            --gene-id-file-url ${genes_url}\
            --output-dir-path ${params.query_raw_data}
    """
}

REF_MAT_URL = Channel.from(params.reference_matrix_url)
REF_BARCODES_URL = Channel.from(params.reference_barcodes_url)
REF_GENES_URL = Channel.from(params.reference_genes_url)
REF_METADATA_URL = Channel.from(params.reference_metadata_url)

process get_reference_data{
    publishDir "${baseDir}/data", mode: 'copy'

    conda 'envs/dropletutils.yaml'
    input:
        val matrix_url from REF_MAT_URL
        val barcodes_url from REF_BARCODES_URL
        val genes_url from REF_GENES_URL
        val metadata_url from REF_METADATA_URL

    output:
        file("${params.reference_raw_data}") into REF_DIR
        file("reference_metadata.txt") into REF_METADATA

    """
    scmap-get-data.R\
            --matrix-file-url ${matrix_url}\
            --barcodes-file-url ${barcodes_url}\
            --gene-id-file-url ${genes_url}\
            --metadata-file-url ${metadata_url}\
            --output-dir-path ${params.reference_raw_data}
    """
}

// produce sce object for query dataset
//query_exp_mat = params.query_raw_data
//QUERY_DIR = Channel.fromPath(params.query_raw_data)

process create_query_sce {
    conda 'envs/dropletutils.yaml'

    // resource handling 
    //errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }   
    //maxRetries 3
    //memory { 2.GB * task.attempt }
    
    input:
        file(query_dir) from QUERY_DIR

    output:
        file("query_sce.rds") into QUERY_SCE

    """
    dropletutils-read-10x-counts.R\
                        -s ${query_dir}\
                        -c ${params.col_names}\
                        -o query_sce.rds
    """ 
}

// produce sce object for reference dataset 
//REF_MAT = Channel.fromPath(params.reference_raw_data)
process create_reference_sce {
    conda 'envs/dropletutils.yaml'
    input:
        file(ref_metadata) from REF_METADATA
        file(ref_dir) from REF_DIR

    output:
        file("reference_sce.rds") into REF_SCE

    """
    dropletutils-read-10x-counts.R\
                -s ${ref_dir}\
                -c ${params.col_names}\
                -m ${ref_metadata}\
                -b ${params.cell_id_col}\
                -f ${params.cell_id_col},${params.cluster_col}\
                -o reference_sce.rds
    """ 
}

// pre-process query dataset 
process preprocess_query_sce {
    conda 'envs/scmap.yaml'
    input:
        file(query_sce) from QUERY_SCE

    output:
        file("query_sce_processed.rds") into QUERY_SCE_PROC 

    """
    scmap-preprocess-sce.R -i ${query_sce} -o query_sce_processed.rds

    """
}

// pre-process reference dataset
process preprocess_ref_sce {
    conda 'envs/scmap.yaml'
    input:
        file(ref_sce) from REF_SCE 

    output:
        file("ref_sce_processed.rds") into REF_SCE_PROC

    """
    scmap-preprocess-sce.R -i ${ref_sce} -o ref_sce_processed.rds
    """
}

// select relevant features for reference dataset 
process select_ref_features {
    publishDir "${baseDir}/data/output", mode: 'copy'
    conda 'envs/scmap.yaml'
    input:
        file(ref_sce) from REF_SCE_PROC
    output:
        file("ref_features.rds") into REF_FEATURES

    """
    scmap-select-features.R\
                -i ${ref_sce}\
                -p ${params.plot_file}\
                -o ref_features.rds
    """
}

projection_method = params.projection_method
REF_CLUSTER = Channel.create()
REF_CELL = Channel.create()
// make a map to re-direct input into correct channel 
channels = ["cluster":0, "cell":1]
REF_FEATURES.choice(REF_CLUSTER, REF_CELL){channels[projection_method]}

// obtain index for cluster-level projections 
process index_cluster {
    conda 'envs/scmap.yaml'
    input:
        file(ref_features_sce) from REF_CLUSTER

    output:
        file("ref_index_cluster.rds") into REF_CLUSTER_INDEX

    """
    scmap-index-cluster.R\
                    -i ${ref_features_sce}\
                    -c ${params.cluster_col}\
                    -o ref_index_cluster.rds
    """
}

// obtain index for cell-level projections 
process index_cell {
    conda 'envs/scmap.yaml'
    input:
        file(ref_features_sce) from REF_CELL

    output:
        file("ref_index_cell.rds") into REF_CELL_INDEX

    """
    scmap-index-cell.R\
                 -i ${ref_features_sce}\
                 -o ref_index_cell.rds
    """

}

// coerce queue channel into value channel for re-using 
QUERY_SCE_PROC = QUERY_SCE_PROC.first()

// obtain cluster-level projections 
process get_cluster_projections{
    conda 'envs/scmap.yaml'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file(ref_cluster_index) from REF_CLUSTER_INDEX
        file(query_sce) from QUERY_SCE_PROC

    output:
        file("cluster_projection_result_sce.rds") into PROJECTED_CLUSTERS_SCE
        file("cluster_projection_result_file.txt") into PROJECTED_CLUSTERS_TXT
        
            """
             scmap-scmap-cluster.R\
                                -i ${ref_cluster_index}\
                                -p ${query_sce}\
                                -r ${params.threshold}\
                                -t cluster_projection_result_file.txt\
                                -o cluster_projection_result_sce.rds
            """
}

// obtain cell-level projections 
process get_cell_projections {
    conda 'envs/scmap.yaml'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file(ref_cell_index) from REF_CELL_INDEX
        file(query_sce) from QUERY_SCE_PROC

    output:
        file("cell_projection_result_file.txt") into PROJECTED_CELLS_TXT 
        file("closest_cells.txt") into CLOSEST_CELLS
        file("closest_cell_similarity.txt") into CLOSEST_CELLS_SIMILARITY
        file("cell_projection_result_object.rds") into PROJECTED_CELLS_SCE


    """
    scmap-scmap-cell.R\
                 -i ${ref_cell_index}\
                 -p ${query_sce}\
                 -c ${params.cluster_col}\
                 -t cell_projection_result_file.txt\
                 -l closest_cells.txt\
                 -s closest_cell_similarity.txt\
                 -o cell_projection_result_object.rds
    """
}

