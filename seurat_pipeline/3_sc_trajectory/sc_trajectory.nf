// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scTrajectory
// Submit example: nextflow run sc_trajectory.nf -params-file sc_trajectory_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------DEFINE GLOBAL PARAMETERS----------------*/
// scripts
code_seurat2monocle =  "$SCRATCH/src/seurat2monocle.R"
code_trajectory = "$SCRATCH/src/trajectory.R"
code_pseudotime = "$SCRATCH/src/pseudotime.R"

/*--------------MAIN WORKFLOW----------------*/
workflow {
    write_params()

    seurat2monocle(code_seurat2monocle,
        params.in_seurat_rds)

    trajectory(code_trajectory,
        seurat2monocle.out)

    // do pseudotime only if parameters to set the root node has been provided
    if ((params.root_node!=null) || ((params.root_metadata.key!=null) && (params.root_metadata.val!=null))) {
        println(params.root_node)
        println(params.root_metadata.key)
        println(params.root_metadata.val)

        pseudotime(code_pseudotime,
            trajectory.out.cds_obj)
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
        path "code.R"
        path file_sc_obj
    output:
        path "*_cds.rds"
    """
    Rscript code.R --file_sc_obj=${file_sc_obj}  --assay=${params.assay} --keep_embeddings=${params.keep_embeddings} --npcs=${params.npcs}
    """
}

// Perform trajectory analyis on the monocle object
process trajectory {
  publishDir params.outfolder+"trajectory", mode: 'copy', overwrite: true
  input:
        path "code.R"
        path file_cds_obj
  output:
        path "*.rds", emit: cds_obj
        path "*.p??"
        path "*.tsv" optional true
  """
  Rscript code.R --file_cds_obj=${file_cds_obj} --filename_prefix=${params.filename_prefix} --clustering_k=${params.clustering_k} --seperate_trajectory_by_partition=${params.seperate_trajectory_by_partition} --close_loop=${params.close_loop} --group_bys=${params.group_bys} --pt_size=${params.pt_size} --width=${params.width} --height=${params.height} --do_de_genenes_in_trajectory=${params.do_de_genenes_in_trajectory}
  """
}

// Make trajectory related paths
process pseudotime {
  publishDir params.outfolder+"pseudotime", mode: 'copy', overwrite: true
  input:
        path "code.R"
        path file_cds_obj
  output:
        path "*.p??"
        path "*.tsv" optional true
        path "*.txt" optional true
  """
  Rscript code.R --file_cds_obj=${file_cds_obj} --filename_prefix=${params.filename_prefix} --root_node=${params.root_node} --root_metadata_key=${params.root_metadata.key} --root_metadata_val=${params.root_metadata.val}  --group_bys=${params.group_bys} --pt_size=${params.pt_size} --width=${params.width} --height=${params.height} --genes=${params.genes} --gex_genes_per_file=${params.gex_genes_per_file} --gex_pt_size=${params.gex_pt_size} --gex_width=${params.gex_width} --gex_height=${params.gex_height}
  """
}
