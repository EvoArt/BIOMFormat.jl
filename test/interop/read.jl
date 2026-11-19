# Reads a table the reference implementation wrote (interop.py) and compares it
# with BIOMFormat's own reading of the original JSON fixture.
#   julia --project=. test/interop/read.jl <outdir>
using BIOMFormat
fx(name) = joinpath(@__DIR__, "..", "fixtures", name)
ref = read_biom(fx("rich_sparse_otu_table.biom"))
p = read_biom(joinpath(ARGS[1], "rich_from_python_hdf5.biom"))
@assert observations(p) == observations(ref)
@assert samples(p) == samples(ref)
@assert counts(p) == counts(ref)
@assert taxonomy(p) == taxonomy(ref)
for k in metadata(ref, :sample)
    @assert sample_metadata(p, k) == sample_metadata(ref, k) k
end
println("ok: rich_from_python_hdf5.biom")
