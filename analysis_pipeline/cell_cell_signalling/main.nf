// Launching an srun (slurm job), and go into $SCRATCH
// module load NextFlow/v22.04.1
// Submit example: nextflow run main.nf -params-file params.yaml -profile codon

include { write_params; set_default_assay; separate_seurat_object} from "./modules/generic_process.nf"
include { differential_analysis_on_single_cell } from "./modules/generic_process.nf"
include { NICHES } from "./workflows/niches/main.nf"
include { LIANA } from "./workflows/LIANA/main.nf"
include { diff_LIANA } from "./workflows/diff_LIANA/main.nf"

workflow {
    write_params()

    // create channel from input and make a meta, file tuple
    sc_obj_ch_merged = Channel.fromList(params.in_seurat_rds)
    | map { it ->
        meta = it.subMap('name')
        [meta, it.file]
    }
    | set_default_assay 
    | view

    // do DE analysis and some metadata plotting on the input Seurat object
    // differential_analysis_on_single_cell(sc_obj_ch_merged, params.separate_rds_key, params.cell_types_column)

    // Run Liana and NICHES on the two conditions separately. Differential NICHES needs separating the two conditions first
    if ( params.run_NICHES || params.run_liana )  {
        // separate the merged morris object into seurat objects on which NICHES will be run separately.
        sc_separate_ch = sc_obj_ch_merged
        | combine(Channel.fromList(params.separate_rds_values))
        | view()
        sc_obj_ch = separate_seurat_object(sc_separate_ch)
        | map {meta, part, filename ->
            meta_new = [name:meta.name, part:part]
            [meta_new, filename]
        }
        // if asked to run niches, call the NICHES workflow
        if ( params.run_NICHES )  {
            NICHES(sc_obj_ch) 
        }
        // if asked to run LIANA, call the NICHES workflow
        if ( params.run_liana ) {
            LIANA(sc_obj_ch)
        }
    }
    
    // Run differential LIANA on the merged dataset
    if ( params.run_diff_liana ) {
        diff_LIANA(sc_obj_ch_merged)
    }

}
