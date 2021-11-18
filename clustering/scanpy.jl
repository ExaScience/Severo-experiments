import Severo: NamedCountMatrix, filter_counts, normalize_cells, find_variable_features, scale_features,
    embedding, shared_nearest_neighbours, cluster, umap, find_all_markers
using DataFrames

using PyCall
pushfirst!(PyVector(pyimport("sys")."path"), "/shared/common/software/Python/3.8.2-GCCcore-9.3.0/lib/python3.8/site-packages")
sc = pyimport("scanpy")

py"""
def convert_matrix(i, p, x, dim, rn, cn):
  from scipy.sparse import csc_matrix
  from anndata import AnnData
  X = csc_matrix((x, i, p), shape=dim)
  return AnnData(X, {'row_names':rn}, {'col_names':cn}, dtype=X.dtype.name)

def variable_features(X):
   hvf = X.var.highly_variable
   rank = X.var.highly_variable_rank[hvf].sort_values()
   return rank.index

def set_idents(X, lbls):
  import pandas as pd
  X.obs['louvain'] = pd.Categorical(lbls)
  return X

def de_to_df(X):
  import pandas as pd
  de = X.uns['rank_genes_groups']
  frames = [pd.melt(pd.DataFrame.from_records(de[n]), var_name='group', value_name=n).set_index('group', append=True) for n in ["names", "logfoldchanges", "pvals", "pvals_adj"]]
  df = pd.DataFrame.join(frames[0], frames[1:]).reset_index('group')
  df.group = df.group.astype('int')
  return df
"""

function matrix_to_py(X::NamedCountMatrix)
  x = X.array
  #py"convert_matrix"(x.rowval .- 1, x.colptr .- 1, convert(Vector{Float32},x.nzval), size(X), names(X,1), names(X, 2))
  py"convert_matrix"(x.rowval .- 1, x.colptr .- 1, x.nzval, size(X), names(X,1), names(X, 2))
end

function filter_counts(X::PyObject; min_cells=3, min_features=200)
    sc.pp.filter_cells(X, min_genes=min_features)
    sc.pp.filter_genes(X, min_cells=min_cells)
    X
end

function normalize_cells(X::PyObject; method=:lognormalize, scale_factor=1e4)
    normalization_method = if method == :lognormalize
        "LogNormalize"
    else
        error("unknown normalization method: $method")
    end

    sc.pp.normalize_total(X, target_sum=scale_factor)
    sc.pp.log1p(X)
    X
end

function find_variable_features(X::PyObject, nfeatures=2000; method=:vst)
    selection_method = if method == :vst
        "seurat_v3"
    else
        error("unknown selection method: $method")
    end

    sc.pp.highly_variable_genes(X, flavor=selection_method, n_top_genes=nfeatures)
    py"variable_features"(X)
end

function scale_features(X::PyObject; scale_max::Real=Inf, features::Union{Nothing, AbstractArray,PyObject}=nothing)
    Y = if features !== nothing
        py"$X[:, $features]"
    else
        X
    end

    sc.pp.scale(Y, max_value=10)
    Y
end

function embedding(X::PyObject, ncomponents::Int64=50; method=:pca, kw...)
    @assert method == :pca
    sc.tl.pca(X, svd_solver="arpack", n_comps=ncomponents)
    X
end

function shared_nearest_neighbours(X::PyObject, k=20; dims=1:50, ntables::Int=50, seed=nothing)
    if seed !== nothing
      sc.pp.neighbors(X, n_neighbors=k, n_pcs=maximum(dims), random_state=seed)
    else
      sc.pp.neighbors(X, n_neighbors=k, n_pcs=maximum(dims))
    end
    X
end

function cluster(X::PyObject; algorithm=:louvain, resolution=0.8, nrandomstarts=10, niterations=10, verbose=false, group_singletons::Bool=true, seed=nothing)
    @assert algorithm == :louvain
    if seed !== nothing
        sc.tl.louvain(X, resolution=resolution, random_state=seed)
    else
        sc.tl.louvain(X, resolution=resolution)
    end
    X.obs.louvain.astype("int")
end

function umap(X::PyObject, ncomponents::Int64=2; dims=:, metric=:cosine,
        nneighbours::Integer=30, min_dist::Real=0.3, nepochs::Union{Nothing,Integer}=nothing)
    if metric isa Symbol
        metric = String(metric)
    end

    sc.pp.neighbors(X, n_neighbors=nneighbours, metric=cosine, n_pcs=maximum(dims))
    sc.tl.umap(X, n_components=ncomponents, maxiter=nepochs, min_dist=min_dist, init_pos="spectral")
    py"$(X).obsm['X_umap']"
end

function pd_to_df(df_pd)
    df= DataFrame()
    for col in df_pd.columns
        df[!, col] = if col == "names"
          convert(Vector{String}, getproperty(df_pd, col).values)
        else
          getproperty(df_pd, col).values
        end
    end
    df
end

function find_all_markers(X::PyObject, idents; logfc_threshold::Real=0.25, min_pct::Real=0.1,
         min_diff_pct::Real=-Inf, only_pos::Bool=false, log::Bool=false, method=:wilcoxon)

    test_use = if method == :wilcoxon
        "wilcox"
    else
        error("unknown differential expression method: $method")
    end

    if idents !== nothing
        py"set_idents"(X, idents)
    end

    sc.tl.rank_genes_groups(X, "louvain", method="wilcoxon", tie_correct=true)
    df = py"de_to_df"(X)
    pd_to_df(df)
end
