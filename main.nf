#!/usr/bin/env nextflow 


// produce sce object for query dataset
QUERY_DIR = Channel.from(params.query_10x_dir)
process create_query_sce {
    conda "${baseDir}/envs/dropletutils.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }
    
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
REF_DIR = Channel.from(params.reference_10x_dir)
REF_METADATA = Channel.from(params.reference_metadata)
process create_reference_sce {
    conda "${baseDir}/envs/dropletutils.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }
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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }
    maxRetries 10
    memory { 16.GB * task.attempt }

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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

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
    conda "${baseDir}/envs/scmap.yaml"

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

    publishDir "${params.output_dir_cluster}", mode: 'copy'

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
    conda "${baseDir}/envs/scmap.yaml"
    publishDir "${params.output_dir_cell}", mode: 'copy'

    errorStrategy { task.attempt < 10  ? 'retry' : 'finish' }   
    maxRetries 10
    memory { 16.GB * task.attempt }

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

// combine the output of two channels 
//PROJECTED_CELLS_TXT
//    .concat(PROJECTED_CLUSTERS_TXT)
//    .set(ALL_PROJECTIONS)


process get_final_output_cluster{
    conda "${baseDir}/envs/dropletutils.yaml"
    publishDir "${params.results_dir}", mode: 'copy'  

    input:
        file(predictions) from PROJECTED_CLUSTERS_TXT

    output:
        file("scmap-cluster_output.txt") into FINAL_TABLE_CLUSTERS

    """
    get_workflow_output.R\
                --predictions-file ${predictions}\
                --workflow-output scmap-cluster_output.txt
    """

}

process get_final_output_cell{
    conda "${baseDir}/envs/dropletutils.yaml"
    publishDir "${params.results_dir}", mode: 'copy'  

    input:
        file(predictions) from PROJECTED_CELLS_TXT

    output:
        file("scmap-cell_output.txt") into FINAL_TABLE_CELLS

    """
    get_workflow_output.R\
                --predictions-file ${predictions}\
                --workflow-output scmap-cell_output.txt
    """

}
