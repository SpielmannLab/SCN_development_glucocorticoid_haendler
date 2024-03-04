// To convert from Seurat to a scanpy h5ad file
process seurat_to_scanpy { 
    label 'medium_memory'
    tag "${meta.name}_${meta.part}"
    input:
    tuple val(meta), path("sc_obj.rds")
    output:
    tuple val(meta), path("*.h5ad")
    """
    tree
    rds_h5ad.R --scobj="sc_obj.rds"
    tree
    """
}

process do_liana {
    // Convert single-cell object to single-cell-signalling object using basic NICHES
    debug true
    label 'medium_memory'
    tag "${meta.name}_${meta.part}"
    input:
    tuple val(meta), path("scanpy_obj.h5ad")
    output:
    tuple val(meta), path("figures")
    tuple val(meta), path("scores")
    path "liana.log" optional true
    path "liana_help.out" optional true
    shell:
    '''
    do_liana.py --input_file "scanpy_obj.h5ad" --celltype_column !{params.cell_types_column} --source_celltypes "!{params.source_celltypes_to_plot}" --target_celltypes "!{params.target_celltypes_to_plot}" --methodsToRun "!{params.methodsToRun}" --resource_name !{params.liana_resource_name} --de_method !{params.de_method} --dotplot_width !{params.liana_figure_width} --dotplot_height !{params.liana_figure_height} --tileplot_width !{params.liana_figure_width}  --tileplot_height !{params.liana_figure_height}
    '''
}
