###############################################
### Function: doing clustering using Seurat ###
###############################################
doSeuratIntegration <- function(obj, nfeatures = 2500, resolution = 1, k.filter = 200, correctCC = FALSE, n_dim = 30, min.dist = 0.75, n.comp = 3){

    if(length(table(obj$orig.ident))!=1){

        obj.list <- SplitObject(object = obj, split.by = "orig.ident")

        for (i in 1:length(x = obj.list)) {
            obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
            obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                                  selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }

        reference.list <- obj.list[names(table(obj$orig.ident))]
        obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter)

        obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)

        # switch to integrated assay. The variable features of this assay are
        # automatically set during IntegrateData
        DefaultAssay(object = obj.integrated) <- "integrated"

        obj <- obj.integrated

    } else {

        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)

    }

    if(correctCC == TRUE){
        obj <- ScaleData(object = obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj), verbose = FALSE)
    } else {
        obj <- ScaleData(object = obj, verbose = FALSE)
    }

    obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:n_dim, min.dist = min.dist, n.components = n.comp)
    obj <- RunTSNE(object = obj, reduction = "pca", dims = 1:30)

    return(obj)
}

#####################################################
### Function: finding ancestor node for each node ###
#####################################################

createLineage_Knn <- function(emb, pd, replication_times=500, removing_cells_ratio=0.2, k_neigh = 5, permute=FALSE){

    print(dim(emb))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}



    res = list()
    rep_i = 1

    while(rep_i < (replication_times+1)){

      sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))

      emb_sub = emb[sampling_index,]
      pd_sub = pd[sampling_index,]

      irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
      irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
      pd_sub1 <- pd_sub[pd_sub$day == "pre",]
      pd_sub2 <- pd_sub[pd_sub$day == "nex",]

      if(permute){
        pd_sub1$state <- pd_sub1$Anno[sample(1:nrow(pd_sub1))]
        pd_sub2$state <- pd_sub2$Anno[sample(1:nrow(pd_sub2))]
      } else {
        pd_sub1$state <- pd_sub1$Anno
        pd_sub2$state <- pd_sub2$Anno
      }

      pre_state_min = min(table(as.vector(pd_sub1$state)))

      if (pre_state_min < k_neigh & pre_state_min >= 3){
          k_neigh = pre_state_min
          print(k_neigh)
      }

      if (pre_state_min < 3){
          next
      }

      neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index

      tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
      for(i in 1:k_neigh){
          tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
      }
      state1 <- names(table(as.vector(pd_sub1$state)))
      state2 <- names(table(as.vector(pd_sub2$state)))

      tmp2 <- matrix(NA,length(state2),length(state1))
      for(i in 1:length(state2)){
          x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
          for(j in 1:length(state1)){
              tmp2[i,j] <- sum(x==state1[j])
          }
      }
      tmp2 <- tmp2/apply(tmp2,1,sum)
      tmp2 <- data.frame(tmp2)
      row.names(tmp2) = state2
      names(tmp2) = state1

      res[[rep_i]] = tmp2

      rep_i = rep_i + 1

    }

    # Calculate the median from the replications

    #### creating the median value of the matrix
    state_1 = row.names(res[[1]])
    state_2 = names(res[[1]])
    tmp_1 = matrix(NA,nrow(res[[1]]),ncol(res[[1]]))
    for(i in 1:nrow(res[[1]])){
        for(j in 1:ncol(res[[1]])){
            xx = NULL
            for(k in 1:replication_times){
                xx = c(xx, res[[k]][i,j])
            }
            tmp_1[i,j] = median(xx[!is.na(xx)])
        }
    }
    tmp_1 = data.frame(tmp_1)
    row.names(tmp_1) = state_1
    names(tmp_1) = state_2

    res <- tmp_1

    return(res)
}


# Get links for sankey for the knn-lineage. Deoes not run on omics cluster
make_sankey_links <- function(knn_matrix){
  library(webshot)
  library(networkD3)

  # Transform matrix to connection data frame with tidyr from the tidyverse:
  links <- knn_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="target") %>%
    tidyr::gather(key="source", value="value", -1) %>%
    filter(value != 0)

  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>%
    unique()
  )

  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1

  return(links)
}

# Merge multiple seurat objects and perform scaling.
perform_merge <- function(sc_list, project="MyProject", nhvg, covars=NULL){
	sc <- merge(x=sc_list[[1]],
		y=sc_list[-1],
		project=project,
		merge.data=TRUE
		)
	sc <- FindVariableFeatures(sc,
		assay = "RNA",
		selection.method = "vst",
		nfeatures = nhvg,
		verbose=FALSE
		)
	sc <- ScaleData(sc,
		features = rownames(sc),
		assay = "RNA",
		vars.to.regress = covars,
		verbose=FALSE
		)
	return(sc)
}
