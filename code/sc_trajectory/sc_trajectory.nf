// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scTrajectory
// Submit example: nextflow run sc_trajectory.nf -params-file sc_trajectory_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------MAIN WORKFLOW----------------*/
workflow {
    write_params()

    seurat2monocle(params.in_seurat_rds)

    trajectory(seurat2monocle.out)

    // do pseudotime only if parameters to set the root node has been provided
    if ((params.root_node!=null) || ((params.root_metadata.key!=null) && (params.root_metadata.val!=null))) {
        println(params.root_node)
        println(params.root_metadata.key)
        println(params.root_metadata.val)

        pseudotime(trajectory.out.cds_obj)

        pseudotime_in_seurat(params.in_seurat_rds, 
            pseudotime.out.cds_obj_with_pseudotime)
    }
}

/*--------------PROCESSES----------------*/

process write_params{ // write input parameters to file
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        path "*.txt"
    """
    echo \$(date) > parameters_sc_trajectory_${params.id}.txt
    echo "Parameters used to perform trajectory analysis:" >> parameters_sc_trajectory_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_sc_trajectory_${params.id}.txt
    """
}

// Convert seurat objects to monocle object. Choose whether or not to keep embeddings
process seurat2monocle {
    input:
        path file_sc_obj
    output:
        path "*_cds.rds"
    """
    seurat2monocle.R --file_sc_obj=${file_sc_obj}  --assay=${params.assay} --keep_embeddings=${params.keep_embeddings} --npcs=${params.npcs}
    """
}

// Perform trajectory analyis on the monocle object
process trajectory {
  publishDir params.outfolder+"trajectory", mode: 'copy', overwrite: true
  input:
        path file_cds_obj
  output:
        path "*.rds", emit: cds_obj
        path "*.p??"
        path "*.tsv" optional true
  """
  trajectory.R --file_cds_obj=${file_cds_obj} --filename_prefix=${params.filename_prefix} --clustering_k=${params.clustering_k} --seperate_trajectory_by_partition=${params.seperate_trajectory_by_partition} --close_loop=${params.close_loop} --group_bys=${params.group_bys} --pt_size=${params.pt_size} --width=${params.width} --height=${params.height} --do_de_genenes_in_trajectory=${params.do_de_genenes_in_trajectory}
  """
}

// Make trajectory related paths
process pseudotime {
  publishDir params.outfolder+"pseudotime", mode: 'copy', overwrite: true
  input:
        path file_cds_obj
  output:
        path "*_pseudotime.rds", emit: cds_obj_with_pseudotime
        path "*.p??"
        path "*.tsv" optional true
        path "*.txt" optional true
  """
  pseudotime.R --file_cds_obj=${file_cds_obj} --filename_prefix=${params.filename_prefix} --root_node=${params.root_node} --root_metadata_key=${params.root_metadata.key} --root_metadata_val=${params.root_metadata.val}  --group_bys=${params.group_bys} --pt_size=${params.pt_size} --width=${params.width} --height=${params.height} --genes=${params.genes} --gex_genes_per_file=${params.gex_genes_per_file} --gex_pt_size=${params.gex_pt_size} --gex_width=${params.gex_width} --gex_height=${params.gex_height}
  """
}

// Make trajectory related paths in Seurat as per my taste
// Note: coloring of the cells ins done 
process pseudotime_in_seurat {
  publishDir params.outfolder+"pseudotime_in_seurat", mode: 'copy', overwrite: true
  input:
        path file_sc_obj
        path file_cds_obj
  output:
        path "*.rds"
        path "*.p??"
        path "*.tsv" optional true
        path "*.txt" optional true
  """
  pseudotime_in_seurat.R --file_sc_obj=${file_sc_obj} --filename_prefix=${params.filename_prefix} --assay=${params.assay} --file_cds_obj=${file_cds_obj} --genes=${params.genes} --gex_pt_size=${params.gex_pt_size} --gex_width=${params.gex_width} --gex_height=${params.gex_height}
  """
}
