"""
    BIOMFormat

Read and write BIOM tables, the Biological Observation Matrix format of QIIME,
QIIME 2, phyloseq and most other amplicon tools: both BIOM 1.0 (JSON) and BIOM
2.1 (HDF5).

```julia
using BIOMFormat
t = read_biom("feature-table.biom")     # format detected from the file
counts(t)                               # observations × samples, sparse
taxonomy(t)                             # per-observation lineage, if present
write_biom("out.biom", t)               # BIOM 2.1; `format = :json` for 1.0
```

With DenoiseDNA loaded, `BIOMTable(seqtab)` and `SequenceTable(t)` convert to
and from its sequence tables.
"""
module BIOMFormat

using SparseArrays
using Printf: @sprintf
import HDF5
import JSON
import Tables
import DataAPI: metadata    # extend DataAPI's (and DataFrames') metadata, so they don't clash

export BIOMTable, read_biom, write_biom
export observations, metadata, observation_metadata, sample_metadata, group_metadata
# `counts`, `samples` and `taxonomy` are not exported: CommunityEcology and other
# packages export different functions with these names. Use
# `using BIOMFormat: counts, samples, taxonomy`, or qualify them (`BIOMFormat.counts(t)`).

include("style.jl")
include("table.jl")
include("json.jl")
include("hdf5.jl")
include("io.jl")
include("tables.jl")
include("show.jl")

end
