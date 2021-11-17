#!/usr/bin/env Rscript
library(Seurat)
library(feather)
args <- commandArgs(trailingOnly = TRUE)

dat <- readRDS(file = args[1])
metadata <- as.data.frame(x = read_feather(path = args[2]))
rownames(x = metadata) <- metadata$sample_id
metadata <- metadata[colnames(x = dat), ]

seurat_object <- CreateSeuratObject(counts = dat, meta.data = metadata)

if(length(args) < 3) {
	outfile <- paste0(tools::file_path_sans_ext(args[1]), ".h5ad")
} else {
    outfile <- args[3]
}

library(sceasy)
sceasy::convertFormat(seurat_object, from="seurat", to="anndata", outFile=outfile)
