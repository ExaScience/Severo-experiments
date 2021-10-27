using HDF5
using CSV
using DataFrames

const steps_jl = ["filter_counts", "normalize_cells", "find_variable_features", "scale_features", "embedding", "shared_nearest_neighbours", "cluster", "umap", "find_all_markers"]
const steps_R = ["CreateSeuratObject", "NormalizeData", "FindVariableFeatures", "ScaleData", "RunPCA", "FindNeighbors", "FindClusters", "RunUMAP", "FindAllMarkers"]

function read_dict(ds)
    k = read(ds, "keys")
    v = read(ds, "vals")
    replace!(k, "normalize"=>"normalize_cells")
    Dict(zip(k,v))
end

function read_time(ds, steps)
    if haskey(ds, "t")
        d = read_dict(ds["t"])
        [get(d, k, missing) for k in steps]
    else
        [missing for k in steps]
    end
end

function read_mem(ds, steps)
    if haskey(ds, "mem")
        d = read_dict(ds["mem"])
        [get(d, k, missing) for k in steps]
    else
        [missing for k in steps]
    end
end

function find_cellcount(f)
    for ds in f
        if haskey(ds, "lbls")
            return length(ds["lbls"])
        end
    end

    error("cell count not found")
end

df = DataFrame(dataset=String[], it=Int64[], size=Int64[], implementation=String[],
     step=String[], time=Union{Missing, Float64}[], mem=Union{Missing, Float64}[], clusters=Int64[])
for fname in filter(endswith(".h5"), readdir())
    try
        name, _ = splitext(fname)
        h5open(fname, "r") do io
            it = if startswith(name, "subsample")
                parts = split(name, "_")
                it = parse(Int64, parts[3])
                size = parse(Int64, parts[5])
                it
            else
                1
            end

            size = find_cellcount(io)

            for (prefix, steps) in [("R", steps_R), ("jl", steps_jl), ("jl32", steps_jl), ("py", steps_R), ("jl32_opt", steps_jl), ("jl_opt", steps_jl)]
                if haskey(io, prefix)
                    ds = io[prefix]
                    lbls = read(ds, "lbls")
                    time = read_time(ds, steps)
                    mem = read_mem(ds, steps)
                    append!(df, DataFrame(dataset=name, it=it, size=size, implementation=prefix, step=steps_R, time=time, mem=mem, clusters=length(unique(lbls))))
                end
            end
            df
        end
    catch e
        println(stderr, "Failed to process $fname")
        showerror(stderr, e)
    end
end
CSV.write(stdout, df)
