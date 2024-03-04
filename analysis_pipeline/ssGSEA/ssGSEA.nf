// Launch an srun (slurm job), and go into $SCRATCH
// module load nextflow/v22.04.1
// Submit example: nextflow run ssGSEA.nf -params-file ssGSEA_params.yaml --id ${SLURM_JOB_ID}

/*--------------MAIN WORKFLOW----------------*/
workflow {
    write_params()

    sc_obj_channel = Channel.fromList(params.in_seurat_rds)
        .view()

    downsample_seurat_object(sc_obj_channel)

    gsea(downsample_seurat_object.out, params.geneSets)

    heatmap_gsea(gsea.out, params.groupby)

    if (params.do_statistics){
        statistics_gsea(gsea.out)
        if (params.do_volcano) {
            volcano_gsea(statistics_gsea.out)
        }
    }

    if (params.plot_heatmap_per_pathway) {
        heatmap_gsea_per_pathway(gsea.out)
    }

    if (params.plot_vlnplot_per_pathway) {
        vlnplot_gsea_per_pathway(gsea.out)
    }
}

/*--------------PROCESSES----------------*/
process write_params{ // write input parameters to file
    output:
        path "*.txt"
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Parameters of the job" >> parameters_sc_p_sample_${suffix}.txt
    echo "!{params}" | tr , '\n' >> parameters_sc_p_sample_${suffix}.txt
    '''
}

// To downsample the Seurat object to run faster
process downsample_seurat_object { 
    input:
        path file_sc_obj
    output:
        path "*_downsampled*.rds"
    """
    downsample_seurat_object.R --file_sc_obj=${file_sc_obj} --downsample=${params.downsample}
    """
}

// ssGSEA score calculation - only run if it hasn't already been run this runs for two days
process gsea {
    input:
        path file_sc_obj
        each geneSet
    output:
        file "*.rds"
    script:
    """
    analysis_gsea.R --file_sc_obj=${file_sc_obj} --geneSet="${geneSet}" --assay=${params.assay} --species="${params.species}" --outfolder="${params.outfolder}gsea"
    """
}

// GSEA Plots
// heatmap. 1 plot for all pathways together
process heatmap_gsea {
    input:
        file file_sc_gsea_metadata
        each groupby
    output:
        file "*.p??"
    script:
        """
        heatmap_gsea.R --file_sc_gsea_metadata=${file_sc_gsea_metadata} --metadata_to_plot="ESCAPE_" --groupby=${groupby}
        """
}

// GSEA statistics
// Do statistics and output TSV file
process statistics_gsea {
    input:
        file file_sc_gsea_metadata
    output:
        file "*.tsv"
    script:
        """
        statistics_gsea.R --file_sc_gsea_metadata=${file_sc_gsea_metadata} --metadata_to_analyse="ESCAPE_" --groupby=${params.statistics_groupby} --groups_to_compare=${params.statistics_groups_to_compare} --splitby=${params.statistics_splitby} --multiple_comparison_adjustment_method=${params.statistics_multiple_comparisons_adjustment_method}
        """
}

// GSEA Plots
process volcano_gsea {
    input:
        file file_gsea_statistics_with_median
    output:
        file "*.pdf"
    script:
        """
        volcano_gsea.R --file_gsea_statistics_with_median=${file_gsea_statistics_with_median} --width=${params.volcano_width} --height=${params.volcano_height} --pt_size=${params.volcano_pt_size}
        """
}

// GSEA Plots
// heatmap. 1 plot per pathway
process heatmap_gsea_per_pathway {
    input:
        file file_sc_gsea_metadata
    output:
        file "*.p??"
    script:
        """
        heatmap_gsea_per_pathway.R --file_sc_gsea_metadata=${file_sc_gsea_metadata} --metadata_to_plot="ESCAPE_" --heatmap_column=${params.heatmap_column} --heatmap_row=${params.heatmap_row}
        """
}

// GSEA Plots
// violin plot. 1 plot per pathway
process vlnplot_gsea_per_pathway {
    input:
        file file_sc_gsea_metadata
    output:
        file "*.p??"
        file "*.tsv"
    script:
        """
        vlnplot_gsea_per_pathway.R --file_sc_gsea_metadata=${file_sc_gsea_metadata} --metadata_to_plot="ESCAPE_" --vln_groupby=${params.vln_groupby} --vln_splitby=${params.vln_splitby} --pt_size=${params.vln_pt_size} --vln_width=${params.vln_width} --vln_height=${params.vln_height} --vln_colors=${params.vln_colors}
        """
}
