# BIOM 1.0: one JSON object. Specification:
# https://biom-format.org/documentation/format_versions/biom-1.0.html

# Parsed JSON as plain Julia values: Dict{String,Any}, Vector{String} where
# every element is a string, Vector{Any} otherwise.
_plain(x) = x
_plain(x::AbstractDict) = Dict{String,Any}(string(k) => _plain(v) for (k, v) in x)
function _plain(x::AbstractVector)
    v = Any[_plain(y) for y in x]
    !isempty(v) && all(y -> y isa AbstractString, v) ? String.(v) : v
end
_plain(x::AbstractString) = String(x)

_str(x) = x === nothing ? "" : String(x)

function _read_json(io::IO)
    j = _plain(JSON.parse(read(io, String)))
    get(j, "format", nothing) === nothing && throw(ArgumentError("not a BIOM file: no \"format\" field"))
    rows, cols = j["rows"], j["columns"]
    n, m = Int.(j["shape"])
    length(rows) == n && length(cols) == m ||
        throw(ArgumentError("BIOM shape $((n, m)) does not match $(length(rows)) rows and $(length(cols)) columns"))
    mtype = get(j, "matrix_type", "sparse")
    eltype_name = get(j, "matrix_element_type", "float")
    eltype_name == "unicode" && throw(ArgumentError("BIOM tables of unicode values are not supported"))
    T = eltype_name == "int" ? Int : Float64
    data = j["data"]
    S = if mtype == "sparse"
        I = Int[]; J = Int[]; V = T[]
        for e in data
            push!(I, Int(e[1]) + 1); push!(J, Int(e[2]) + 1); push!(V, T(e[3]))
        end
        sparse(I, J, V, n, m)
    elseif mtype == "dense"
        sparse(T[data[r][c] for r in 1:n, c in 1:m])
    else
        throw(ArgumentError("unknown BIOM matrix_type \"$mtype\""))
    end
    BIOMTable(_narrow(S), [_str(r["id"]) for r in rows], [_str(c["id"]) for c in cols];
              id = _str(get(j, "id", "")), type = _str(get(j, "type", "")),
              generated_by = _str(get(j, "generated_by", "")), creation_date = _str(get(j, "date", "")),
              observation_metadata = _json_metadata(rows), sample_metadata = _json_metadata(cols))
end

# One Dict per id; `nothing` for the axis when no id has any.
function _json_metadata(entries)
    md = [get(e, "metadata", nothing) for e in entries]
    all(x -> x === nothing || isempty(x), md) && return nothing
    [x === nothing ? Dict{String,Any}() : Dict{String,Any}(x) for x in md]
end

function _write_json(io::IO, t::BIOMTable{T}) where {T}
    entries(ids, md) = [Dict("id" => id, "metadata" => (md === nothing || isempty(md[k]) ? nothing : md[k]))
                        for (k, id) in enumerate(ids)]
    I, J, V = findnz(t.data)
    o = Dict(
        "id" => t.id,
        "format" => "1.0.0",
        "format_url" => "http://biom-format.org",
        "type" => t.type,
        "generated_by" => t.generated_by,
        "date" => t.creation_date,
        "matrix_type" => "sparse",
        "matrix_element_type" => T <: Integer ? "int" : "float",
        "shape" => [size(t, 1), size(t, 2)],
        "rows" => entries(t.observations, t.observation_metadata),
        "columns" => entries(t.samples, t.sample_metadata),
        "data" => [[i - 1, j - 1, v] for (i, j, v) in zip(I, J, V)])
    isempty(t.observation_group_metadata) && isempty(t.sample_group_metadata) ||
        @warn "BIOM 1.0 (JSON) has no group metadata; it was not written. Use format = :hdf5 to keep it."
    write(io, JSON.json(o))
    nothing
end
