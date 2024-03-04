include { diff_seurat_to_scanpy; do_diff_liana } from "../../modules/diff_liana_processes.nf"

workflow diff_LIANA {
    take:
    sc_obj_ch_merged

    main:
    sc_obj_ch_merged | diff_seurat_to_scanpy | do_diff_liana

}
