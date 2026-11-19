# The table: what a BIOM file holds, in the file's own orientation
# (observations as rows, samples as columns).

"""
    BIOMTable

A BIOM table: counts (or other values) of observations (OTUs, ASVs, genes, ...)
in samples, with optional metadata on each.

- `id`, `type`, `generated_by`, `creation_date`: the file's header. `type` is
  BIOM's controlled vocabulary, e.g. `"OTU table"`.
- `observations`, `samples`: the ids, in row and column order.
- `data`: the values, a sparse `observations × samples` matrix (`Int` when
  every value is a whole number, else `Float64`).
- `observation_metadata`, `sample_metadata`: one `Dict{String,Any}` per id,
  or `nothing` when the file has none. A lineage is stored as a vector of
  ranks under `"taxonomy"`.
- `observation_group_metadata`, `sample_group_metadata`: values describing a
  whole axis, e.g. a phylogeny, as `name => (data_type, value)`.

Construct one with `BIOMTable(data, observations, samples; kwargs...)`, or
read one with [`read_biom`](@ref).
"""
struct BIOMTable{T<:Real}
    id::String
    type::String
    generated_by::String
    creation_date::String
    observations::Vector{String}
    samples::Vector{String}
    data::SparseMatrixCSC{T,Int}
    observation_metadata::Union{Nothing,Vector{Dict{String,Any}}}
    sample_metadata::Union{Nothing,Vector{Dict{String,Any}}}
    observation_group_metadata::Dict{String,Tuple{String,String}}
    sample_group_metadata::Dict{String,Tuple{String,String}}

    function BIOMTable{T}(id, type, generated_by, creation_date, observations, samples, data,
                          observation_metadata, sample_metadata, observation_group_metadata,
                          sample_group_metadata) where {T<:Real}
        n, m = size(data)
        length(observations) == n ||
            throw(DimensionMismatch("$(length(observations)) observation ids for $n rows"))
        length(samples) == m ||
            throw(DimensionMismatch("$(length(samples)) sample ids for $m columns"))
        allunique(observations) || throw(ArgumentError("observation ids are not unique"))
        allunique(samples) || throw(ArgumentError("sample ids are not unique"))
        observation_metadata === nothing || length(observation_metadata) == n ||
            throw(DimensionMismatch("$(length(observation_metadata)) observation metadata entries for $n observations"))
        sample_metadata === nothing || length(sample_metadata) == m ||
            throw(DimensionMismatch("$(length(sample_metadata)) sample metadata entries for $m samples"))
        new{T}(id, type, generated_by, creation_date, observations, samples, data,
               observation_metadata, sample_metadata, observation_group_metadata, sample_group_metadata)
    end
end

"""
    BIOMTable(data, observations, samples; id = "", type = "OTU table",
              generated_by = "BIOMFormat.jl", creation_date = now,
              observation_metadata = nothing, sample_metadata = nothing,
              observation_group_metadata = Dict(), sample_group_metadata = Dict())

A table from an `observations × samples` matrix (dense or sparse) and the ids.
Metadata is one `Dict` per id, or a vector of `Dict`s keyed by anything
`string` turns into a key.
"""
function BIOMTable(data::AbstractMatrix{<:Real}, observations::AbstractVector, samples::AbstractVector;
                   id::AbstractString = "", type::AbstractString = "OTU table",
                   generated_by::AbstractString = "BIOMFormat.jl",
                   creation_date::AbstractString = _now(),
                   observation_metadata = nothing, sample_metadata = nothing,
                   observation_group_metadata = Dict{String,Tuple{String,String}}(),
                   sample_group_metadata = Dict{String,Tuple{String,String}}())
    S = sparse(data)
    T = eltype(S) <: Integer ? Int : Float64
    BIOMTable{T}(String(id), String(type), String(generated_by), String(creation_date),
                 String.(string.(observations)), String.(string.(samples)),
                 SparseMatrixCSC{T,Int}(S), _metadata(observation_metadata), _metadata(sample_metadata),
                 Dict{String,Tuple{String,String}}(observation_group_metadata),
                 Dict{String,Tuple{String,String}}(sample_group_metadata))
end

_metadata(::Nothing) = nothing
_metadata(md::AbstractVector) = [Dict{String,Any}(string(k) => v for (k, v) in pairs(d)) for d in md]

# ISO 8601, to the second, without a Dates dependency
_now() = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())

# Values that are all whole numbers become Int, as BIOM's "int" element type.
function _narrow(S::SparseMatrixCSC)
    all(isinteger, nonzeros(S)) || return SparseMatrixCSC{Float64,Int}(S)
    SparseMatrixCSC{Int,Int}(size(S)..., copy(S.colptr), copy(S.rowval), Int.(nonzeros(S)))
end

Base.size(t::BIOMTable) = size(t.data)
Base.size(t::BIOMTable, d) = size(t.data, d)
Base.:(==)(a::BIOMTable, b::BIOMTable) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(BIOMTable))

"""
    counts(t::BIOMTable) -> SparseMatrixCSC

The values, observations × samples.
"""
counts(t::BIOMTable) = t.data

"""
    observations(t::BIOMTable) -> Vector{String}
    samples(t::BIOMTable) -> Vector{String}

The observation (row) and sample (column) ids.
"""
observations(t::BIOMTable) = t.observations
samples(t::BIOMTable) = t.samples

"""
    observation_metadata(t, key) -> Vector
    sample_metadata(t, key) -> Vector

The value of one metadata field for every observation or sample, with
`missing` where an id has none.
"""
observation_metadata(t::BIOMTable, key) = _field(t.observation_metadata, length(t.observations), key)
sample_metadata(t::BIOMTable, key) = _field(t.sample_metadata, length(t.samples), key)

_field(::Nothing, n, key) = fill(missing, n)
_field(md, n, key) = [get(d, string(key), missing) for d in md]

"""
    metadata(t, axis) -> Vector{String}

The metadata fields present on `axis` (`:observation` or `:sample`).
"""
function metadata(t::BIOMTable, axis::Symbol)
    md = axis === :observation ? t.observation_metadata :
         axis === :sample ? t.sample_metadata :
         throw(ArgumentError("axis is :observation or :sample, got :$axis"))
    md === nothing && return String[]
    sort!(unique!(String[k for d in md for k in keys(d)]))
end

"""
    taxonomy(t::BIOMTable) -> Vector{Vector{String}}

Each observation's lineage from its `"taxonomy"` metadata, as a vector of ranks
(empty where it has none). A lineage stored as one `;`-separated string is split.
"""
function taxonomy(t::BIOMTable)
    map(observation_metadata(t, "taxonomy")) do x
        x === missing || x === nothing ? String[] :
        x isa AbstractString ? String.(strip.(split(x, ';'))) : String.(string.(x))
    end
end

"""
    group_metadata(t, axis) -> Dict{String,Tuple{String,String}}

Metadata describing a whole axis (`:observation` or `:sample`), such as a
phylogeny: `name => (data_type, value)`, e.g. `"phylogeny" => ("newick", "(...);")`.
"""
group_metadata(t::BIOMTable, axis::Symbol) =
    axis === :observation ? t.observation_group_metadata :
    axis === :sample ? t.sample_group_metadata :
    throw(ArgumentError("axis is :observation or :sample, got :$axis"))

"""
    t[obs, samp]

The table restricted to some observations and samples: indices, ranges,
logical masks, or vectors of ids. Metadata follows the ids it belongs to.
"""
function Base.getindex(t::BIOMTable{T}, obs, samp) where {T}
    oi = _index(obs, t.observations)
    si = _index(samp, t.samples)
    BIOMTable{T}(t.id, t.type, t.generated_by, t.creation_date, t.observations[oi], t.samples[si],
                 t.data[oi, si],
                 t.observation_metadata === nothing ? nothing : t.observation_metadata[oi],
                 t.sample_metadata === nothing ? nothing : t.sample_metadata[si],
                 copy(t.observation_group_metadata), copy(t.sample_group_metadata))
end

_index(::Colon, ids) = eachindex(ids)
_index(i::Union{Integer,AbstractRange,AbstractVector{<:Integer},AbstractVector{Bool}}, ids) =
    collect(eachindex(ids))[i isa Integer ? [i] : i]
function _index(names::AbstractVector{<:AbstractString}, ids)
    pos = Dict(id => k for (k, id) in enumerate(ids))
    [haskey(pos, n) ? pos[n] : throw(KeyError(n)) for n in names]
end
_index(name::AbstractString, ids) = _index([name], ids)
