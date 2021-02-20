library(monocle3)
library(Seurat)
library(argparse)
library(stringr)
####################################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}



# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, state=NULL,index="celltype"){
  cell_ids <- which(colData(cds)[, index] == state)

  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}


reduceDim_and_cluster <- function(cds,alignment_group=NULL){
	message("INFO : Preprocess PCA ...")
        cds <- preprocess_cds(cds,
                      num_dim = 50,
                      method="PCA",
                      norm_method="log")


        message("INFO : Align ...")
        cds <- align_cds(cds,
                 preprocess_method="PCA",
                 alignment_k=20,
                 residual_model_formula_str="~Size_Factor+num_genes_expressed",
                 alignment_group=alignment_group)
	
	message("INFO : reduce dimension ...")
	message("INFO : run tSNE ...")
        cds <- reduce_dimension(cds,reduction_method="tSNE",preprocess_method="Aligned",cores=8)

	message("INFO : run UMAP ...")
        cds <- reduce_dimension(cds,reduction_method="UMAP",preprocess_method="Aligned",cores=8)
	
	message("INFO : Cluster")
	cds<-cluster_cells(cds,
                   reduction_method="UMAP",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)
	
	cds<-cluster_cells(cds,
                   reduction_method="tSNE",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

	message("INFO : Learn graph ...")
	cds<-learn_graph(cds,
                 use_partition=TRUE,
                 close_loop=TRUE)

	return(cds)
}


reduceDim3D <- function(cds){
	cds_3d <- reduce_dimension(cds, reduction_method="UMAP",max_components = 3,preprocess_method="Aligned",cores=8)
        cds_3d <- cluster_cells(cds_3d,reduction_method="UMAP")
        cds_3d <- learn_graph(cds_3d)
	return(cds_3d)
}

#####################################
parser <- ArgumentParser(description='Program to handle monocle3 pstime')
parser$add_argument("--cds",
                    type="character",
                    default=NULL,
                    help="monocle3 object")

parser$add_argument("--outdir",
                    type="character",
                    default="./outdir")

parser$add_argument("--again",
		    action="store_true",
		    default=FALSE)

parser$add_argument("--alignment",
                    type="character",
                    default=NULL,
		    help="alignment group in align_cds function")


parser$add_argument("--index",
		    type="character",
		    default=NULL,
		    help="column in colData to run Pseudotime")

parser$add_argument("--state",
                    type="character",
                    default=NULL,
                    help="select the state in column as root state")

args <- parser$parse_args()

outDir=args$outdir
flag=args$again
group=args$alignment
index=args$index
state=args$state

makedir(outDir)
message("INFO : Loading dataset")
cds=readRDS(args$cds)

if(flag){
	message("INFO : run Cluster and ReduceDim again ...")
	cds=reduceDim_and_cluster(cds,group)
}

assertthat::assert_that(!is.null(state) | !is.null(index))
nodes=get_earliest_principal_node(cds,index=index,state=state)

message("INFO : Order Cells and Generate Pseudotime Variable")
cds = order_cells(cds, root_pr_nodes=nodes)

message("INFO : Run UMAP 3D")
cds_3d=reduceDim3D(cds)
nodes=get_earliest_principal_node(cds_3d,index=index,state=state)
cds_3d = order_cells(cds_3d, root_pr_nodes=nodes)

message("INFO : Done and Save ...")
saveRDS(cds,file.path(outDir,"Pseudotime-cds.rds"))
saveRDS(cds_3d,file.path(outDir,"Pseudotime-cds_3d.rds"))
