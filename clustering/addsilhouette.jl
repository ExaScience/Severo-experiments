using HDF5
using Severo
using Statistics

CACHE = Dict{String, Any}()

for result in ARGS
    println("Processing $result")
    h5open(result, "cw") do io
        for (implementation, g) in zip(keys(io), io)
            for (dataset, d) in zip(keys(g), g)
                X = get!(CACHE, dataset) do
                    println("calculating for $dataset")
                    if dataset == "sce_full_Zhengmix4uneq"
                        X = read_h5("/data/thaber/DuoClustering2018/sce_full/sce_full_Zhengmix4uneq.h5")
                        X = filter_counts(X, min_cells=3, min_features=0)
                    elseif dataset == "3k"
                        X = read_10X("/data/thaber/PBMC/3k")
                        X = filter_counts(X, min_cells=3, min_features=200)
                    else
                        X = read_h5ad(string("/data/thaber/M1/",dataset,".h5ad"))
                        X = filter_counts(X, min_cells=3, min_features=200)
                    end

                    hvf = find_variable_features(X, 2000, method=:vst)
                    Y = normalize_cells(X, scale_factor=1e4, method=:lognormalize)

                    S = scale_features(Y, scale_max=10, features=hvf)
                    em = embedding(S, 50, method=:pca, algorithm=:tssvd)
                    em.coordinates
                end

                for ds in d
                    if ds isa HDF5.Dataset
                        name = basename(HDF5.name(ds))
                        move_link(d, name, d, "tmp")
                        x = create_group(d, name)
                        move_link(d, "tmp", x, "lbls")

                        for k in keys(attributes(ds))
                            v = read_attribute(ds, k)
                            write_attribute(x, k, v)
                            delete_attribute(ds, k)
                        end

                        ds = x
                    end

                    haskey(ds, "width") && continue

                    lbls = read(ds, "lbls")
                    if implementation == "scanpy"
                        lbls .+= 1
                    end

                    dims = read_attribute(ds, "dims")
                    width = approx_silhouette(X[:, 1:dims], lbls)
                    write(ds, "width", width)
                    write_attribute(ds, "avgsilh", mean(width))
                end
            end
        end
    end
end
