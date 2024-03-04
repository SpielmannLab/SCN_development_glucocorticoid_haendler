workflow CELLCELLSIGNALLINGR {
    sc_obj_channel = Channel.fromList(params.in_seurat_rds)
        .view()
    
    if (params.run_SingleCellSignalR) {
        // Extract normalized count from the Seurat object
        seurat_to_matrix(sc_obj_channel, "RNA", "data")

        // perform the SingleCellSignalR analysis
        singlecellsignalr(seurat_to_matrix.out, species, min_logFC)

        // make heatmap of LR interaction strenth for every LR-pair and cell-cell pair
        singlecellsignalr_plots(singlecellsignalr.out.file_signal, min_LRscore, sampling)
    }

}

process seurat_to_matrix { 
    conda "${WORK}/opt/anaconda3/envs/scVelocity"
    cpus 1
    memory '2 GB'
    publishDir path: "${params.outfolder}${task.process.replaceAll(':', '-')}", mode: 'copy', overwrite: true
    input:
        path file_sc_obj
    output:
        tuple path("*_matrix.rds"), path("*_clusters.rds"), path("*_annotations.rds")
    """
    seurat_to_matrix.R --file_sc_obj=${file_sc_obj} --assay=${params.assay} --slot=${params.slot} --cluster_column=${cluster_column} --annotation_column=${annotation_column}
    """
}

// Run the singlecellsignal network analysis
process singlecellsignalr { 
    conda "${WORK}/opt/anaconda3/envs/SingleCellSignalR"
    cpus 1
    memory '2 GB'
    publishDir path: "${params.outfolder}${task.process.replaceAll(':', '-')}", mode: 'copy', overwrite: true
    input:
        tuple path(file_matrix), path(file_clusters), path(file_annotations)
        val species
        val min_logFC
    output:
        path "*_all.pdf"
        path "*_LRinteractions.rds", emit: file_signal

    """
    singlecellsignalr.R --file_matrix=${file_matrix} --file_clusters=${file_clusters} --file_annotations=${file_annotations} --species="${species}" --min_logFC=${min_logFC}
    mv images/* ./
    """
}

process singlecellsignalr_plots {
    conda "${WORK}/opt/anaconda3/envs/MegaBundle"
    cpus 1
    memory '1 GB'
    publishDir path: "${params.outfolder}${task.process.replaceAll(':', '-')}", mode: 'copy', overwrite: true
    input:
        path file_signal
        val min_LRscore
        val sampling
    output:
        path("*_heatmap.pdf")
    """
    singlecellsignalr_plots.R --file_signal=${file_signal} --min_LRscore=${min_LRscore} --sampling=${sampling}
    """
}

