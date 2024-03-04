include { seurat_to_scanpy; do_liana } from "../../modules/liana_processes.nf"

workflow LIANA {
    take:
    sc_obj_ch

    main:
    sc_obj_ch | seurat_to_scanpy | do_liana 

}
