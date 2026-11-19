# Writes tables for the reference implementation to read (interop.py).
#   julia --project=. test/interop/write.jl <outdir>
using BIOMFormat
out = ARGS[1]
fx(name) = joinpath(@__DIR__, "..", "fixtures", name)
t = read_biom(fx("rich_sparse_otu_table.biom"))
write_biom(joinpath(out, "rich_from_julia_hdf5.biom"), t)
write_biom(joinpath(out, "rich_from_julia_json.biom"), t; format = :json)
g = read_biom(fx("rich_sparse_otu_table_hdf5_group_metadata.biom"))
write_biom(joinpath(out, "group_from_julia_hdf5.biom"), g)
f = BIOMTable([0.5 0.0 2.0; 1.25 2.0 0.0], ["obsA", "obsB"], ["s1", "s2", "s3"];
              sample_metadata = [Dict("depth" => 10.0, "site" => "x"), Dict("depth" => 2.5, "site" => "y"),
                                 Dict("depth" => 4.0, "site" => "z")])
write_biom(joinpath(out, "float_from_julia_hdf5.biom"), f)
println("written")
