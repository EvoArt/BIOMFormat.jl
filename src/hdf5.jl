# BIOM 2.1: an HDF5 file. Specification:
# https://biom-format.org/documentation/format_versions/biom-2.1.html
#
# The values are stored twice, as CSR under observation/matrix and as CSC under
# sample/matrix. HDF5 is row-major, so an (N, k) dataset reads into Julia as a
# k × N array: a 2-D metadata field (a taxonomy) has one column per id here.

_attr(a, k, default = "") = haskey(a, k) ? _scalar(a[k]) : default      # attrs(x)[k] is the value
_scalar(x::AbstractArray) = isempty(x) ? "" : _scalar(first(x))
_scalar(x::AbstractString) = String(x)
_scalar(x) = x

function _read_hdf5(path::AbstractString)
    HDF5.h5open(path, "r") do f
        a = HDF5.attrs(f)
        shape = haskey(a, "shape") ? Int.(a["shape"]) : nothing
        obs = String.(read(f["observation/ids"]))
        samp = String.(read(f["sample/ids"]))
        n, m = length(obs), length(samp)
        shape === nothing || Tuple(shape) == (n, m) ||
            throw(ArgumentError("BIOM shape $(Tuple(shape)) does not match $n observation and $m sample ids"))
        # CSC over samples, as triplets: sparse() sorts rows within each column
        # and sums duplicates, whatever order the writer used.
        indptr = Int.(read(f["sample/matrix/indptr"]))
        J = Int[j for j in 1:m for _ in (indptr[j] + 1):indptr[j + 1]]
        I = Int.(read(f["sample/matrix/indices"])) .+ 1
        V = Float64.(read(f["sample/matrix/data"]))
        S = sparse(I, J, V, n, m)
        BIOMTable(_narrow(S), obs, samp;
                  id = string(_attr(a, "id")), type = string(_attr(a, "type")),
                  generated_by = string(_attr(a, "generated-by")),
                  creation_date = string(_attr(a, "creation-date")),
                  observation_metadata = _h5_metadata(f, "observation", n),
                  sample_metadata = _h5_metadata(f, "sample", m),
                  observation_group_metadata = _h5_group_metadata(f, "observation"),
                  sample_group_metadata = _h5_group_metadata(f, "sample"))
    end
end

# One Dict per id, from the datasets under <axis>/metadata; `nothing` if none.
function _h5_metadata(f, axis, n)
    path = "$axis/metadata"
    haskey(f, path) || return nothing
    g = f[path]
    names = collect(keys(g))
    isempty(names) && return nothing
    md = [Dict{String,Any}() for _ in 1:n]
    for name in names
        v = read(g[name])
        if ndims(v) == 1
            length(v) == n || throw(ArgumentError("metadata \"$name\" has $(length(v)) values for $n ids"))
            for i in 1:n
                md[i][name] = v[i] isa AbstractString ? String(v[i]) : v[i]
            end
        elseif ndims(v) == 2
            size(v, 2) == n || throw(ArgumentError("metadata \"$name\" has $(size(v, 2)) rows for $n ids"))
            for i in 1:n
                col = String.(v[:, i])
                while !isempty(col) && isempty(col[end])      # rows are padded to one length
                    pop!(col)
                end
                md[i][name] = col
            end
        else
            throw(ArgumentError("metadata \"$name\" has $(ndims(v)) dimensions"))
        end
    end
    md
end

function _h5_group_metadata(f, axis)
    out = Dict{String,Tuple{String,String}}()
    path = "$axis/group-metadata"
    haskey(f, path) || return out
    g = f[path]
    for name in keys(g)
        ds = g[name]
        dt = haskey(HDF5.attrs(ds), "data_type") ? string(_scalar(HDF5.attrs(ds)["data_type"])) : ""
        out[name] = (dt, string(_scalar(read(ds))))
    end
    out
end

function _write_hdf5(path::AbstractString, t::BIOMTable)
    S = SparseMatrixCSC{Float64,Int}(t.data)
    St = sparse(transpose(S))              # CSR of S: columns of St are rows of S
    HDF5.h5open(path, "w") do f
        a = HDF5.attrs(f)
        a["id"] = t.id
        a["type"] = t.type
        a["format-url"] = "http://biom-format.org"
        a["format-version"] = Int64[2, 1]
        a["generated-by"] = t.generated_by
        a["creation-date"] = t.creation_date
        a["shape"] = Int64[size(t, 1), size(t, 2)]
        a["nnz"] = Int64(nnz(S))
        for (axis, ids, M, md, gmd) in (("observation", t.observations, St, t.observation_metadata,
                                         t.observation_group_metadata),
                                        ("sample", t.samples, S, t.sample_metadata, t.sample_group_metadata))
            g = HDF5.create_group(f, axis)
            g["ids"] = ids
            mg = HDF5.create_group(g, "matrix")
            mg["data"] = Vector{Float64}(nonzeros(M))
            mg["indices"] = Int32.(rowvals(M) .- 1)
            mg["indptr"] = Int32.(M.colptr .- 1)
            _h5_write_metadata(HDF5.create_group(g, "metadata"), md)
            gg = HDF5.create_group(g, "group-metadata")
            for (name, (dt, value)) in gmd
                gg[name] = [value]              # a one-element dataset, as the reference writes
                HDF5.attrs(gg[name])["data_type"] = dt
            end
        end
    end
    nothing
end

# One dataset per field: a lineage-like field (vectors) as a padded 2-D string
# dataset, numbers as Float64 (NaN where an id has none), anything else as
# strings ("" where an id has none).
function _h5_write_metadata(g, md)
    md === nothing && return
    for key in sort!(unique!(String[k for d in md for k in keys(d)]))
        vals = [get(d, key, missing) for d in md]
        present = filter(!ismissing, vals)
        if !isempty(present) && all(x -> x isa AbstractVector, present)
            L = maximum(length, present)
            g[key] = [ismissing(x) || j > length(x) ? "" : string(x[j]) for j in 1:L, x in vals]
        elseif !isempty(present) && all(x -> x isa Real && !(x isa Bool), present)
            g[key] = [ismissing(x) ? NaN : Float64(x) for x in vals]
        else
            g[key] = [ismissing(x) || x === nothing ? "" : string(x) for x in vals]
        end
    end
end
