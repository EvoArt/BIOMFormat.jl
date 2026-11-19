# Conversions between DenoiseDNA's sequence tables and BIOM tables. A
# SequenceTable is samples × sequences; BIOM puts observations in rows, so
# each conversion transposes.
module BIOMFormatDenoiseDNAExt

using BIOMFormat
using DenoiseDNA
using SparseArrays

"""
    BIOMTable(st::SequenceTable; ids = :sequence, kwargs...)

A BIOM table of DenoiseDNA's sequence table. With `ids = :sequence` (the
default, as dada2 tables are used in phyloseq and QIIME 2) the sequences are
the observation ids; with `ids = :asv` the observations are `ASV1`, `ASV2`, ...
and each sequence is kept as `"sequence"` metadata. Other keywords go to the
`BIOMTable` constructor (`id`, `type`, `sample_metadata`, ...).
"""
function BIOMFormat.BIOMTable(st::DenoiseDNA.SequenceTable; ids::Symbol = :sequence, kwargs...)
    obs, md = if ids === :sequence
        st.sequences, nothing
    elseif ids === :asv
        ["ASV$k" for k in eachindex(st.sequences)], [Dict{String,Any}("sequence" => s) for s in st.sequences]
    else
        throw(ArgumentError("ids is :sequence or :asv, got :$ids"))
    end
    kw = merge((; generated_by = "DenoiseDNA.jl", observation_metadata = md), kwargs)
    BIOMTable(sparse(permutedims(st.table)), obs, st.samples; kw...)
end

"""
    SequenceTable(t::BIOMTable)

DenoiseDNA's sequence table of a BIOM table of counts. The sequences are the
observations' `"sequence"` metadata where every observation has one, and
otherwise the observation ids.
"""
function DenoiseDNA.SequenceTable(t::BIOMTable)
    eltype(counts(t)) <: Integer ||
        throw(ArgumentError("a SequenceTable holds counts; this BIOM table has non-integer values"))
    md = observation_metadata(t, "sequence")
    seqs = all(x -> x isa AbstractString, md) ? String.(md) : copy(observations(t))
    DenoiseDNA.SequenceTable(copy(samples(t)), seqs, Matrix{Int}(permutedims(counts(t))))
end

end
