import Severo: NamedCountMatrix, filter_counts, normalize_cells, find_variable_features, scale_features,
    embedding, shared_nearest_neighbours, cluster, umap, find_all_markers

import RCall: @rimport, @R_str, RObject, getclass
@rimport Seurat

function matrix_to_R(X::NamedCountMatrix)
  x = X.array
	R"""
(function(i, p, x, dim, rn, cn) {
  A <- new("dgCMatrix")
  A@i <- i
  A@p <- p
  A@x <- x
  A@Dim <- dim
  rownames(A) <- rn
  colnames(A) <- cn
  Matrix::t(A)
})($(x.rowval .- 1), $(x.colptr .- 1), $(convert(Vector{Float64},x.nzval)), $(collect(size(x))), $(names(X,1)), $(names(X, 2)))
	 """
end

umap_coords = R"""
function(data) data[["umap"]]@cell.embeddings
"""

function filter_counts(A::RObject; min_cells=3, min_features=200)
    #Seurat.CreateSeuratObject(counts=A, var"min.cells"=min_cells, var"min.features"=min_features)
    R"Seurat::CreateSeuratObject(counts=$A, min.cells=$min_cells, min.features=$min_features)"
end

function normalize_cells(data::RObject; method=:lognormalize, scale_factor=1e4)
    normalization_method = if method == :lognormalize
        "LogNormalize"
    elseif method == :relativecounts
        "R"
    else
        error("unknown normalization method: $method")
    end

    #Seurat.NormalizeData(object=data, var"normalization.method"=normalization_method, var"scale.factor"=scale_factor)
    R"Seurat::NormalizeData(object=$data, normalization.method=$normalization_method, scale.factor=$scale_factor)"
end

function find_variable_features(data::RObject, nfeatures=2000; method=:vst)
    selection_method = if method == :vst
        "vst"
    elseif method == :dispersion
        "disp"
    elseif method == :meanvarplot
        "mvp"
    else
        error("unknown selection method: $method")
    end

    #data = Seurat.FindVariableFeatures(object=data, var"selection.method"=selection_method, nfeatures=nfeatures)
    d = R"Seurat::FindVariableFeatures(object=$data, selection.method=$selection_method, nfeatures=$nfeatures)"
    data.p = d.p
    R"Seurat::VariableFeatures(object=$data)"
end

function scale_features(data::RObject; scale_max::Real=Inf, features::Union{Nothing, AbstractArray, RObject}=nothing)
    #Seurat.ScaleData(object=data, features=features, var"scale.max"=scale_max)
    R"Seurat::ScaleData(object=$data, features=$features, scale.max=$scale_max)"
end

function embedding(data::RObject, ncomponents::Int64=50; method=:pca, kw...)
    @assert method == :pca
    #Seurat.RunPCA(object=data, npcs=ncomponents, kw...)
    R"Seurat::RunPCA(object=$data, npcs=$ncomponents)"
end

function shared_nearest_neighbours(data::RObject, k=20; dims=1:50, ntables=50)
    #Seurat.FindNeighbors(object=data, var"k.param"=k, dims=dims, var"n.trees"=ntables)
    R"Seurat::FindNeighbors(object=$data, k.param=$k, dims=$dims, n.trees=$ntables)"
end

function cluster(data::RObject; algorithm=:louvain, resolution=0.8, nrandomstarts=10, niterations=10, verbose=false, group_singletons::Bool=true)
    #data = Seurat.FindClusters(object=data, resolution=resolution, var"n.start"=nrandomstarts,
    #    var"n.iter"=niterations, verbose=verbose, var"group.singletons"=group_singletons)
    #Seurat.Idents(object=data)
    data = R"Seurat::FindClusters(object=$data, resolution=$resolution, n.start=$nrandomstarts,
        n.iter=$niterations, verbose=$verbose, group.singletons=$group_singletons)"
    R"Seurat::Idents(object=$data)"
end

function umap(data::RObject, ncomponents::Int64=2; dims=:, metric=:cosine,
        nneighbours::Integer=30, min_dist::Real=0.3, nepochs::Union{Nothing,Integer}=nothing)
    if metric isa Symbol
        metric = String(metric)
    end

    #Seurat.RunUMAP(object=data, dims=dims, var"n.neighbors"=nneighbours, var"n.epochs"=nepochs,
    #    metric=metric, var"n.components"=ncomponents, var"min.dist"=min_dist)
    data = R"Seurat::RunUMAP(object=$data, dims=$dims, n.neighbors=$nneighbours, n.epochs=$nepochs,
        metric=$metric, n.components=$ncomponents, min.dist=$min_dist)"
    umap_coords(data)
end

function find_all_markers(data::RObject, idents::RObject; logfc_threshold::Real=0.25, min_pct::Real=0.1,
         min_diff_pct::Real=-Inf, only_pos::Bool=false, log::Bool=false, method=:wilcoxon)

    test_use = if method == :wilcoxon
        "wilcox"
    elseif method == :t
        "t"
    else
        error("unknown differential expression method: $method")
    end

    #Seurat.FindAllMarkers(data, var"only.pos"=only_pos, var"min.pct"=min_pct, var"logfc.threshold"=logfc_threshold,
    #    var"min.diff.pct"=min_diff_pct, var"test.use"=test_use)
    R"Seurat::FindAllMarkers($data, only.pos=$only_pos, min.pct=$min_pct, logfc.threshold=$logfc_threshold,
        min.diff.pct=$min_diff_pct, test.use=$test_use)"
end
