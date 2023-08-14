// Usage : 
// module load nextflow/v22.04.1
// nextflow run /data/humangen_mouse/2023_glucocorticoid/code/TOME.nf -params-file /data/humangen_mouse/2023_glucocorticoid/code/TOME.yaml
// do this either after launching an srun or sbatch, and going into $SCRATCH
// do this in the conda environment AstizSCN

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

workflow {
    preprocess(channel.fromList(params.annotated_file_sc_obj), params.annotations_column, params.sex)
    
    // Form pairs for integration and lineage connections and create a single channel with file pairs.
    gd17_pnd2 = preprocess.out.gd17
        .concat( preprocess.out.pnd2 )
        .collect()
    pnd2_pnd10 = preprocess.out.pnd2
        .concat( preprocess.out.pnd10 )
        .collect()
    pnd10_pnd30 = preprocess.out.pnd10
        .concat( preprocess.out.pnd30 )
        .collect()
    channelFilePairs = gd17_pnd2
        .mix( pnd2_pnd10,pnd10_pnd30 )
        .view()

    integrate_pairs(channelFilePairs, params.file_functions)

    print(params.ks)

    connect_pairs(integrate_pairs.out.ch_integrated, params.file_functions, params.ks)

    //collect all the sankey links to make combined sankey plots
    collected_sankey_links=connect_pairs.out.ch_connection_sankey_links
        .collect() 
    
    combined_sankeyplot(collected_sankey_links, params.ks)
}

process preprocess {
    conda '/work/sreenivasan/.omics/anaconda3/envs/AstizSCN'
    publishDir params.outfolder+"TOMEby"+params.annotations_column+"/"+'processed', mode: 'copy', overwrite: true 
    input:
        path "file_sc_obj.rds"
        val annotations_column 
        val sex
    output:
        path "gd17_processed.rds", optional: true, emit: gd17
        path "pnd2_processed.rds", optional: true, emit: pnd2
        path "pnd10_processed.rds", optional: true, emit: pnd10
        path "pnd30_processed.rds", optional: true, emit: pnd30
        path "*.pdf", optional: true
        path "*.tsv", optional: true
        path "output.log"
    """
    ls -lh *
    which Rscript
    preprocess.R ${annotations_column} ${sex} |& tee output.log
    tree
    """
}

process integrate_pairs {
    conda '/work/sreenivasan/.omics/anaconda3/envs/AstizSCN'
    publishDir params.outfolder+"TOMEby"+params.annotations_column+"/"+'integrated_pairs', mode: 'copy', overwrite: true // The output figures will be saved here.
    input:
        tuple path("file_sc_obj_1.rds"), path("file_sc_obj_2.rds") 
        path "functions.R"
    output:
        path "*.png"
        tuple path("*_integrated.rds"), path("*_umap_coords.rds"), path("*_pca_coords.rds"), emit: ch_integrated
    """
    ls -lh *
    integrate_pairs.R 
    tree
    """
}

// connect pairs using the integrated pairs - Create knn lineage matrix
process connect_pairs {
        conda '/work/sreenivasan/.omics/anaconda3/envs/AstizSCN'
    publishDir params.outfolder+"TOMEby"+params.annotations_column+"/"+'connected', mode: 'copy', overwrite: true // The output figures will be saved here.
    input:
        tuple path('file_sc_obj.rds'), path('file_umap_coords.rds'), path('file_pca_coords.rds')
        path 'functions.R'
        each k // testing if results are robust by different k parameter
    output:
        path "*lineage.rds"
        path "*sankey_links.rds", emit:ch_connection_sankey_links
        path "*.txt"
    """
    echo "Running for k-value of $k"
    connect_pairs.R $k
    tree
    """
}

// Step2_3_Combined_SankeyPlots
process combined_sankeyplot {
    conda '/work/sreenivasan/.omics/anaconda3/envs/AstizSCN'
    publishDir params.outfolder+"TOMEby"+params.annotations_column+"/"+'combined_sankeyplots', mode: 'copy', overwrite: true // The output figures will be saved here.
    input:
        file collected_sankey_links
        each k // do this individually for each k
    output:
        path "*.html"
        path "*.pdf"
    """
    tree # Check if all the files have been added by Nextflow
    combined_sankeyplot.R $k ${params.min_link_strength}
    """
}

/*
// Step2_4_Connection with permutation of annotations
code="/data/humangen_mouse/AstizSCN/code/Step2_4_connection_w_permute.R"
ks=[5]
process step2_4_connection_w_permute {
  publishDir path_outfolder+'2_4_connection_w_permute', mode: 'copy', overwrite: true // The output figures will be saved here.
  input:
  tuple path('file_sc_obj.rds'), path('file_umap_coords.rds'), path('file_pca_coords.rds') from channelIntegrated_copy2
  path 'functions.R' from file_functions
  path 'code.R' from code
  each k from ks // testing if results are robust by different k parameter
  output:
  path "*lineage.rds" into channelConnection_lineages_w_permute
  path "*sankey_links.rds" into channelConnection_sankey_links_w_permute
  path "*.txt" into channelConnection_edge_weight_prop_w_permute
  """
  echo "Running for k-value of $k"
  Rscript code.R $k
  tree
  """
}

// Step2_5_Combined_SankeyPlots with permutation
collected_sankey_links_w_permute=channelConnection_sankey_links_w_permute.collect() //collect all the sankey links to make combined sankey plots
code="/data/humangen_mouse/AstizSCN/code/Step2_3_combined_sankeyplots.R"
process step2_5_combined_sankeyplots_w_permute {
  publishDir path_outfolder+'2_5_combined_sankeyplots_w_permute', mode: 'copy', overwrite: true // The output figures will be saved here.
    input:
    path 'code.R' from code
    file collected_sankey_links_w_permute
    each k from ks // do this individually for each k
    output:
    path "*.html" into channelCollectedsankeylinks_w_permute
    """
    tree # Check if all the files have been added by Nextflow
    Rscript code.R $k
    """
}
*/