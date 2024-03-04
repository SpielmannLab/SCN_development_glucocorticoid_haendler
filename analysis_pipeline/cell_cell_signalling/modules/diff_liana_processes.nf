// To convert from Seurat to a scanpy h5ad file
process diff_seurat_to_scanpy { 
    label 'medium_memory'
    tag "${meta.name}"
    input:
    tuple val(meta), path("sc_obj.rds")
    output:
    tuple val(meta), path("*.h5ad")
    """
    rds_h5ad_no_scale_data.R --scobj="sc_obj.rds"
    """
}

process do_diff_liana {
    // Convert single-cell object to single-cell-signalling object using basic NICHES
    debug true
    label 'medium_memory'
    tag "${meta.name}"
    input:
    tuple val(meta), path("scanpy_obj.h5ad")
    output:
    tuple val(meta), path("figures")
    tuple val(meta), path("scores")
    shell:
    '''
    do_diff_liana.py --input_file "scanpy_obj.h5ad" --celltype_column !{params.cell_types_column} --sample_column !{params.sample_column} --condition_column !{params.separate_rds_key} --condition_ref_level !{params.separate_rds_values[0]} --condition_alt_level=!{params.separate_rds_values[1]} --resource_name !{params.diff_liana_resource_name} --tileplot_width !{params.diff_liana_figure_width}  --tileplot_height !{params.diff_liana_figure_height}
    '''
}
