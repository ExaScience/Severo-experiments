using ArgParse
using Severo
using HDF5
using DataFrames
import CategoricalArrays: CategoricalVector, CategoricalValue
import Clustering: randindex, counts

include("seurat.jl")
include("scanpy.jl")
include("named_product.jl")
using .NamedProduct

function pipeline(X; resolution=0.5, k=20, ntables=50, dims=50, seed=nothing)
    X = filter_counts(X, min_cells=3, min_features=200)
    hvf = find_variable_features(X, 2000, method=:vst)

    Y = normalize_cells(X, scale_factor=1e4, method=:lognormalize)#, dtype=Float32)

    S = scale_features(Y, scale_max=10, features=hvf)
    em = embedding(S, dims, method=:pca, algorithm=:tssvd)
    #@assert eltype(em.coordinates) == Float32
    snn = shared_nearest_neighbours(em, k, dims=1:dims, ntables=ntables, seed=seed)
    lbls = cluster(snn, resolution=resolution, seed=seed)
    lbls
end

function read_dataframe(g::Union{HDF5.File, HDF5.Group})
    attrs = HDF5.attributes(g)
    @assert read(attrs, "encoding-version") == "0.1.0"
    @assert read(attrs, "encoding-type") == "dataframe"

    function read_col(n::AbstractString)
        ds = g[n]
        iscategorical = haskey(HDF5.attributes(ds), "categories")
        if iscategorical
            cats = read(HDF5.attributes(ds), "categories")
            cats = getindex(g, cats)
            isordered = read(HDF5.attributes(cats), "ordered") != 0

            vals = read(ds)
            lvls = read(cats)
            v = CategoricalVector(undef, length(vals); levels=lvls, ordered=isordered)
            v .= CategoricalValue.(vals .+ 1, Ref(v.pool))
            v
        else
            read(ds)
        end
    end

    columns = [n => read_col(n) for n in read(attrs, "column-order")]
    df = DataFrame(columns)

    if haskey(attrs, "_index")
        index = read(attrs, "_index")
        index = read_col(index)
        insertcols!(df, 1, :index=>index)
    end

    df
end

function read_labels(fname::AbstractString)
    h5open(fname, "r") do io
        df = read_dataframe(io["obs"])
        df[:,:class_id], df[:,:subclass_id]
    end
end

estimate_mapping(from, to) = map(I->I[2], argmax(counts(from, to), dims=2))

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--output"
            help = "output file"
            default = "results.h5"
        "--ntables"
            help = "number of tables to use"
            arg_type = Union{Int, Nothing}
            default = [nothing]
            nargs = '+'
        "--k"
            arg_type = Int
            default = [20]
            nargs = '+'
        "--dims"
            arg_type = Int
            default = [50]
            nargs = '+'
        "--resolution"
            arg_type = Float64
            default = [0.5]
            nargs = '+'
        "--seed"
            arg_type = Int
            default = [0]
            nargs = '+'
        "--backend"
            arg_type = String
            default = "severo"
            range_tester = in(["severo", "scanpy", "seurat", "severo32"])
        "datasets"
            help = "datasets to process"
            required = true
            arg_type = String
            nargs = '+'
    end

    return parse_args(s)
end

function as_factors(lbls)
    z = Dict(reverse.(enumerate(unique(lbls))))
    getindex.([z], lbls)
end

args = parse_arguments()
for fname in args["datasets"]
#    fname = "/data/thaber/M1/Human_M1_10xV3-D1.h5ad"
    name, _ = splitext(basename(fname))

    X = read_h5ad(fname)
    class_true, lbls_true = read_labels(fname)

    backend = args["backend"]
    for (dims, k, resolution, ntables, seed) in named_product(dims=args["dims"], k=args["k"], resolution=args["resolution"], ntables=args["ntables"], seed=args["seed"])
        println("$(dims) $(k) $(ntables) $(resolution) $(seed)")

        Y = if backend == "seurat"
            matrix_to_R(X)
        elseif backend == "scanpy"
            matrix_to_py(X)
        else
            X
        end

        if ntables === nothing
            ntables = dims * 2
        end

        lbls, time, stats = @timed pipeline(Y, resolution=resolution, k=k, ntables=ntables, dims=dims, seed=seed)
        lbls = convert(Vector{Int}, lbls)

        h5open(args["output"], "cw") do io
            g = if haskey(io, "$backend/$name")
                io["$backend/$name"]
            else
                create_group(io, "$backend/$name")
            end

            dsname = if seed == nothing
                "$(dims)_$(k)_$(ntables)_$(resolution)"
            else
                "$(dims)_$(k)_$(ntables)_$(resolution)_$(seed)"
            end

            if !haskey(g, dsname)
                ds = create_group(g, dsname)
                write(ds, "lbls", lbls)

                write_attribute(ds, "dims", dims)
                write_attribute(ds, "resolution", resolution)
                write_attribute(ds, "k", k)
                write_attribute(ds, "tables", ntables)
                if seed !== nothing
                    write_attribute(ds, "seed", seed)
                end

                write_attribute(ds, "time", time)
                write_attribute(ds, "stats", collect(stats))

                ri = randindex(lbls_true, lbls)
                write_attribute(ds, "randindex", collect(ri))

                p = purity(lbls, lbls_true)
                write_attribute(ds, "purity", p)

                p = purity(lbls, class_true)
                write_attribute(ds, "classpurity", p)

                clusters = length(unique(read(ds)))
                write_attribute(ds, "clusters", clusters)
            else
                println("$dsname already present")
            end
        end
    end
end
