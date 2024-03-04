include { downsample_seurat_object; impute_by_alra; do_niches; do_diff_niches } from "../../modules/niches_processes.nf"
include { subset_seurat_object } from "../../modules/generic_process.nf"

workflow NICHES {
    take:
    sc_obj_ch

    main:
    if ( params.do_downsample ) {
        sc_obj_ch_ds =  downsample_seurat_object(sc_obj_ch)
    } else {
        sc_obj_ch_ds =  sc_obj_ch
    }

    // subset only the required clusters
    if ( params.subset ) {
        sc_obj_ch_ds_sb = subset_seurat_object(sc_obj_ch_ds, params.key_value_pairs)
    } else {
        sc_obj_ch_ds_sb = sc_obj_ch_ds
    }

    if ( params.use_ALRA ) {
       sc_obj_ch_ds_sb_alra  = impute_by_alra(sc_obj_ch_ds_sb) 
    } else {
       sc_obj_ch_ds_sb_alra  = sc_obj_ch_ds_sb 
    }

    // perform the NICHES analysis
    do_niches( sc_obj_ch_ds_sb_alra, params.species, params.LRdatabase, params.cell_types_column, params.min_cells_p_gene, params.min_cells_p_cluster, params.mode_of_analysis, params.signalling_npcs)

    cell2cell_ch = do_niches.out.cell2cell_ch
    | map {meta, filename -> 
        meta_new = meta.subMap('name')
        [meta_new, filename]
    }
    | groupTuple
    | view

    // Combine the two separated RDS files and analyse them
    do_diff_niches(cell2cell_ch, params.separate_rds_key, params.signalling_npcs, params.min_feat_cell2cell)

}
