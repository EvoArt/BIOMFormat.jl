# BIOMFormat.jl

[![CI](https://github.com/EvoArt/BIOMFormat.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/BIOMFormat.jl/actions/workflows/CI.yml)

Read and write [BIOM](https://biom-format.org) tables, the Biological
Observation Matrix format used by QIIME 2, QIIME, phyloseq and most other
amplicon and metagenomics tools. Both versions in use are supported: BIOM 2.1
(HDF5, what QIIME 2 exports) and BIOM 1.0 (JSON).

```julia
using Pkg
Pkg.add("BIOMFormat")
```

```julia
using BIOMFormat

t = read_biom("feature-table.biom")    # either version; the format is read from the file
counts(t)                              # observations × samples, a sparse matrix
observations(t), samples(t)            # the ids
taxonomy(t)                            # each observation's lineage, if the file has one
sample_metadata(t, "BODY_SITE")        # one metadata field across the samples

gut = t[:, sample_metadata(t, "BODY_SITE") .== "gut"]   # subset by index, mask or id
write_biom("gut.biom", gut)                             # BIOM 2.1; format = :json for 1.0
```

A table is also a [Tables.jl](https://github.com/JuliaData/Tables.jl) table,
one row per non-zero value, so `DataFrame(t)` gives its long form.

## With DenoiseDNA

[DenoiseDNA.jl](https://github.com/EvoArt/DenoiseDNA.jl) sequence tables
convert both ways:

```julia
using DenoiseDNA, BIOMFormat

t = BIOMTable(nochim)                  # sequences as observation ids
t = BIOMTable(nochim; ids = :asv)      # ASV1, ASV2, ...; sequences kept as metadata
write_biom("asvs.biom", t)             # ready for `qiime tools import`

seqtab = SequenceTable(read_biom("asvs.biom"))
```

In the REPL a table shows a summary in the same style as DenoiseDNA's
displays: totals, metadata fields, each sample with a strip of its most
abundant observations, and the most abundant observations with their taxa.

The example tables in `test/fixtures` are from the
[biom-format](https://github.com/biocore/biom-format) project, under its BSD
licence (`test/fixtures/LICENSE-biom-format.txt`).
