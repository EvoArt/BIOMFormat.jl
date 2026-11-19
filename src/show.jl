# Displays, in DenoiseDNA's style: `show(io, t)` is one line; the REPL's
# `MIME"text/plain"` display is a headline, a few key-value lines, the
# samples with a strip of their most abundant observations, and the most
# abundant observations. Helpers are in style.jl.

Base.show(io::IO, t::BIOMTable) =
    print(io, "BIOMTable(", _fmt(size(t, 2)), " samples × ", _fmt(size(t, 1)), " observations)")

# a sample's values over the top observations: shade by share of the sample
function _strip(io, col::AbstractVector{<:Real}, tot)
    b = IOBuffer()
    for x in col
        f = tot > 0 ? x / tot : 0.0
        print(b, x == 0 ? '·' : f >= 0.1 ? '█' : f >= 0.01 ? '▓' : f >= 0.001 ? '▒' : '░')
    end
    s = String(take!(b))
    _hascolor(io) || return s
    return replace(s, "·" => _sty(io, "·"; color=_PAL.tail), r"[░▒▓█]+" => m -> _sty(io, m; color=_PAL.spark))
end

# ids that are DNA sequences (DenoiseDNA's ASVs) are shown as sequences
_isdna(s) = length(s) >= 20 && all(c -> c in "ACGTN", s)

# an observation's id, shortened to `w`
_id(io, s, w) = _isdna(s) ? _seq(io, s, w) : length(s) > w ? first(s, w - 1) * "…" : s

# the most specific rank of a lineage that names something ("g__Escherichia",
# not an empty "s__")
function _rank(lineage)
    for r in Iterators.reverse(lineage)
        name = replace(r, r"^[a-zA-Z]__" => "")
        isempty(strip(name)) || return r
    end
    ""
end

function _fields(io, t, axis)
    f = metadata(t, axis)
    isempty(f) && return _note(io, "none")
    k = 6
    s = join(f[1:min(k, end)], _note(io, ", "))
    length(f) > k ? s * _note(io, ", … $(length(f) - k) more") : s
end

function Base.show(io::IO, ::MIME"text/plain", t::BIOMTable{T}) where {T}
    n, m = size(t)
    tot = sum(t.data)
    what = T <: Integer ? "counts" : "total"
    _title(io, "BIOMTable", "$(_fmt(m)) samples × $(_fmt(n)) observations" * _note(io, " · $(_fmt(tot)) $what"))
    _kv(io, "type", isempty(t.type) ? _note(io, "–") : t.type)
    src = join(filter(!isempty, [t.generated_by, t.creation_date]), _note(io, " · "))
    isempty(src) || _kv(io, "made by", src)
    isempty(t.id) || _kv(io, "id", t.id)
    (n == 0 || m == 0) && return
    _kv(io, "non-zero", "$(_fmt(nnz(t.data))) values" * _note(io, " · $(_pct(nnz(t.data) / (n * m))) of the table"))
    colsum = vec(sum(t.data; dims=1))
    rowsum = vec(sum(t.data; dims=2))
    prev = vec(sum(t.data .!= 0; dims=2))
    sc = sort(colsum)
    _kv(io, what, m == 1 ? _fmt(tot) :
                  "$(_fmt(sc[1]))–$(_fmt(sc[end])) per sample" * _note(io, " · median $(_fmt(sc[cld(m, 2)]))"))
    m > 1 && _kv(io, "prevalence", "$(_fmt(count(==(1), prev))) observations in 1 sample" *
                                   _note(io, " · $(_fmt(count(==(m), prev))) in all $m"))
    _kv(io, "metadata", "observations: " * _fields(io, t, :observation) *
                        _note(io, " · ") * "samples: " * _fields(io, t, :sample))
    gm = [(axis, k, dt) for (axis, d) in (("observations", t.observation_group_metadata),
                                         ("samples", t.sample_group_metadata)) for (k, (dt, _)) in d]
    isempty(gm) || _kv(io, "group", join(("$k" * _note(io, " ($dt, $axis)") for (axis, k, dt) in gm), "  "))

    W = _termwidth(io)
    # samples, with a strip over the most abundant observations
    ord = sortperm(rowsum; rev=true)
    k = _nshow(io, m, 12)
    nmw = min(20, maximum(length, t.samples[1:k]))
    rows = Vector{String}[]
    for j in 1:k
        nm = length(t.samples[j]) > nmw ? first(t.samples[j], nmw - 1) * "…" : t.samples[j]
        push!(rows, [nm, _fmt(colsum[j]), _fmt(count(!=(0), view(t.data, :, j)))])
    end
    hdr = ["sample", what, "observations"]
    used = 2 + sum(c -> maximum(_vis(r[c]) for r in vcat([hdr], rows)) + 2, eachindex(hdr))
    nstrip = min(n, W - used - 1)
    for (j, row) in enumerate(rows)
        push!(row, _strip(io, t.data[ord[1:nstrip], j], colsum[j]))
    end
    _section(io, "Samples", min(W, max(used + nstrip, 60)))
    _table(io, vcat(hdr, "observations by total →"), rows; align=[:l, :r, :r, :l])
    _more(io, m, k, "samples")
    println(io, _note(io, "  share of the sample: ") *
                join((_sty(io, g; color=_PAL.spark) * _note(io, " " * s) for (g, s) in
                      (("█", "≥10%"), ("▓", "≥1%"), ("▒", "≥0.1%"), ("░", "less"))), "  ") *
                "  " * _sty(io, "·"; color=_PAL.tail) * _note(io, " absent"))

    # the most abundant observations, with their lineage where there is one
    tax = taxonomy(t)
    hastax = any(!isempty, tax)
    ka = _nshow(io, n, 8)
    _section(io, "Most abundant observations", min(W, 100))
    rows = [[_fmt(rowsum[ord[r]]), _pct(rowsum[ord[r]] / tot), string(prev[ord[r]])] for r in 1:ka]
    hdr = ["$what", "%", "samples"]
    if hastax
        push!(hdr, "taxon")
        for (r, row) in enumerate(rows)
            push!(row, _rank(tax[ord[r]]))
        end
    end
    used = 2 + sum(c -> maximum(_vis(r[c]) for r in vcat([hdr], rows)) + 2, eachindex(hdr))
    idw = max(12, min(W - used - 2, maximum(length, t.observations[ord[1:ka]])))
    rows = [vcat([_id(io, t.observations[ord[r]], idw)], rows[r]) for r in 1:ka]
    _table(io, vcat("id", hdr), rows; align=vcat(:l, fill(:r, 3), hastax ? [:l] : Symbol[]))
    _more(io, n, ka, "observations")
end
