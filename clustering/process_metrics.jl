using HDF5
using CSV
using DataFrames

function each_attribute(parent::Union{HDF5.Dataset, HDF5.Datatype, HDF5.File, HDF5.Group})
    attrs = attributes(parent)
    NamedTuple((Symbol(k), read_attribute(parent, k)) for k in keys(attrs))
end

df = DataFrame()
h5open("results.h5", "cw") do io
    for (implementation, g) in zip(keys(io), io)
        for (dataset, d) in zip(keys(g), g)
            for ds in d
                size = length(ds)

                attrs = each_attribute(ds)
                push!(df, (
                    dataset=dataset, size=size, implementation=implementation,
                    dims=attrs[:dims], k=attrs[:k], resolution=attrs[:resolution], tables=attrs[:tables],
                    time=attrs[:time], purity=attrs[:purity], classpurity=attrs[:classpurity],
                    clusters=attrs[:clusters], ari=attrs[:randindex][1], ri=attrs[:randindex][2]
                ))
            end
        end
    end
end
CSV.write(stdout, df)
