// Launch an srun (slurm job), and go into $SCRATCH
// module load nextflow/v22.04.1
// Submit example: nextflow run ssGSEA.nf -params-file ssGSEA_params.yaml --id ${SLURM_JOB_ID}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------MAIN WORKFLOW----------------*/
workflow {
    write_params()

    sc_obj_channel = Channel.fromList(params.in_seurat_rds)
        .view()

    downsample_seurat_object(sc_obj_channel)

    gsea(downsample_seurat_object.out, params.geneSets)

    heatmap_gsea(gsea.out, params.groupby)
}

/*--------------PROCESSES----------------*/
process write_params{ // write input parameters to file
    cpus 1
    memory '1 GB'
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        path "*.txt"
    """
    echo \$(date) > parameters_sc_p_sample_${params.id}.txt
    echo "Parameters used to analyze the seurat object for each sample:" >> parameters_sc_p_sample_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_sc_p_sample_${params.id}.txt
    """
}

// To downsample the Seurat object to run faster
process downsample_seurat_object { 
    conda "${WORK}/.omics/anaconda3/envs/scVelocity"
    cpus 1
    memory '50 GB'
    publishDir path: "${params.outfolder}downsampled_seurat_object", mode: 'copy', overwrite: true
    input:
        path file_sc_obj
    output:
        path "*_downsampled.rds"
    """
    downsample_seurat_object.R --file_sc_obj=${file_sc_obj} --downsample=${params.downsample}
    """
}

// ssGSEA score calculation - only run if it hasn't already been run this runs for two days
process gsea {
    conda "${WORK}/.omics/anaconda3/envs/ssGSEA"
    cpus 1
    memory '300 GB'
    publishDir path: "${params.outfolder}gsea", mode: 'copy', overwrite: true
    input:
        path file_sc_obj
        each geneSet
    output:
        file "*.rds"
    script:
    """
    analysis_gsea.R --file_sc_obj=${file_sc_obj} --geneSet=${geneSet} --assay=${params.assay} --species="${params.species}" --outfolder="${params.outfolder}gsea"
    """
}

// GSEA Plots
process heatmap_gsea {
  conda "${WORK}/.omics/anaconda3/envs/ssGSEA"
  publishDir path: "${params.outfolder}heatmap", mode: 'copy', overwrite: true
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
