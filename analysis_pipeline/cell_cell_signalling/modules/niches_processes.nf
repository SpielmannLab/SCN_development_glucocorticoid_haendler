// To make sure the total cell numbers are the same across multiple conditions for differential testing
process downsample_seurat_object { 
    label 'medium_memory'
    tag "${meta.name}_${meta.part}"
    input:
    tuple val(meta), path("sc_obj.rds")
    output:
    tuple val(meta), path("*_downsampled*.rds")
    """
    downsample_seurat_object.R --file_sc_obj="sc_obj.rds" --downsample=${params.downsample}
    """
}

process impute_by_alra {
    // Pass a Seurat object containing normalized count matrix ("data") slot
    label 'medium_memory'
    tag "${meta.name}_${meta.part}"
    input:
    tuple val(meta), path("sc_obj.rds")
    output:
    tuple val(meta), path("sc_obj_imputed.rds")
    shell:
    '''
    impute_by_alra.R --assay=!{params.assay}
    '''
}

process do_niches {
    // Convert single-cell object to single-cell-signalling object using basic NICHES
    // Will use the DefaultAssay
    label 'medium_memory'
    tag "${meta.name}_${meta.part}"
    input:
    tuple val(meta), path("sc_obj.rds")
    env species
    env LRdatabase
    env cell_types_column
    env min_cells_p_gene
    env min_cells_p_cluster
    env mode_of_analysis 
    env signalling_npcs 
    output:
    tuple val(meta), path("CellToCell_obj.rds"), emit: cell2cell_ch
    path "*.p??"
    shell:
    '''
    do_niches.R
    '''
}

process do_diff_niches {
    // debug true
    label 'medium_memory'
    tag "${meta.name}"
    input:
    tuple val(meta), path("cell2cell_obj?.rds")
    env group_by
    env signalling_npcs 
    env min_feat_cell2cell
    output:
    path "*.pdf"
    path "*.tsv"
    path "*.rds"
    shell:
    '''
    do_diff_niches.R
    '''
}
