/*--------------PROCESSES----------------*/
process write_params{ 
    // write input parameters to file
    label 'tiny'
    output:
        path "*.txt"
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_sc_p_sample_${suffix}.txt
    echo "!{params}" | tr , '\n' >> parameters_sc_p_sample_${suffix}.txt
    '''
}

process set_default_assay {
    // Set default assay to the assay requested by the user in the params file
    label 'medium_memory'
    tag "${meta.name}"
    input:
    tuple val(meta), path("sc_obj.rds")
    output:
    tuple val(meta), path("sc_obj_default_assay_set.rds")
    shell:
    '''
    set_default_assay.R --assay=!{params.assay}
    '''
}

process subset_seurat_object { // To subset the Seurat object based on key value pairs
    label 'medium_memory'
    tag "${meta.name}"
    input:
        tuple val(meta), path("sc_obj.rds")
        tuple val(subset_key), val(subset_values)
    output:
        tuple val(meta), path("*_subsetted_*.rds" )
    """
    subset_seurat_object.R --file_sc_obj="sc_obj.rds" --subset_key=${subset_key} --subset_values="${subset_values}"
    """
}

process separate_seurat_object { // To subset the Seurat object based on key value pairs
    label 'medium_memory'
    tag "${separate_rds_value}"
    input:
        tuple val(meta), path("sc_obj.rds"), val(separate_rds_value)
    output:
        tuple val(meta), val(separate_rds_value), path("*_subsetted_*.rds" )
    shell:
    '''
    subset_seurat_object.R --file_sc_obj="sc_obj.rds" --subset_key=!{params.separate_rds_key} --subset_values=!{separate_rds_value}
    '''
}

process differential_analysis_on_single_cell { // Perform differential analysis on the input Seurat object to compare differential cell signalling
    // debug true
    label 'medium_memory'
    tag "${meta.name}"
    input:
    tuple val(meta), path("sc_obj.rds")
    env group_by
    env cell_types_column
    output:
    path "*.pdf"
    path "*.tsv"
    shell:
    '''
    differential_single_cell_analysis_n_plotting.R
    '''
}
