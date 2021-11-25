using HDF5
using CSV
using DataFrames

function each_attribute(parent::Union{HDF5.Dataset, HDF5.Datatype, HDF5.File, HDF5.Group})
    attrs = attributes(parent)
    Dict(Symbol(k) => read_attribute(parent, k) for k in keys(attrs))
end

df = DataFrame(dataset=String[], size=Int[], implementation=String[], dims=Int[], k=Int[], resolution=Float64[], tables=Int[],
    time=Float64[], purity=Union{Missing,Float64}[], classpurity=Union{Missing,Float64}[], clusters=Int[],
    ari=Union{Missing,Float64}[], ri=Union{Missing,Float64}[], avgsilh=Union{Missing,Float64}[], seed=Union{Missing,Int}[])
for results in ARGS
    h5open(results, "r") do io
        for (implementation, g) in zip(keys(io), io)
            for (dataset, d) in zip(keys(g), g)
                for ds in d
                    size = if ds isa HDF5.Dataset
                        size = length(ds)
                    else
                        size = length(ds["lbls"])
                    end

                    attrs = each_attribute(ds)
                    ari, ri = if haskey(attrs, :randindex)
                        attrs[:randindex][1:2]
                    else
                        missing, missing
                    end

                    push!(df, (
                        dataset=dataset, size=size, implementation=implementation,
                        dims=attrs[:dims], k=attrs[:k], resolution=attrs[:resolution], tables=attrs[:tables],
                        time=attrs[:time], purity=get(attrs, :purity, missing), classpurity=get(attrs,:classpurity,missing),
                        clusters=attrs[:clusters], ari=ari, ri=ri,
                        avgsilh=get(attrs, :avgsilh, missing), seed=get(attrs,:seed, missing)
                    ))
                end
            end
        end
    end
end
CSV.write(stdout, df)
