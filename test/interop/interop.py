"""Checks BIOMFormat.jl's files with the reference implementation, biom-format.

    python interop.py <fixtures dir> <dir written by write.jl>

Reads each Julia-written table with biom.load_table and compares it with the
reference's own reading of the original fixture; then writes BIOM 2.1 tables
for read.jl to check from the Julia side.
"""
import sys, os
import numpy as np
import biom
from biom import Table
from biom.util import biom_open

fixtures, out = sys.argv[1], sys.argv[2]
ref = biom.load_table(os.path.join(fixtures, "rich_sparse_otu_table.biom"))

def same(a, b, what):
    assert list(a.ids(axis="observation")) == list(b.ids(axis="observation")), what + ": observation ids"
    assert list(a.ids(axis="sample")) == list(b.ids(axis="sample")), what + ": sample ids"
    assert np.array_equal(a.matrix_data.toarray(), b.matrix_data.toarray()), what + ": values"
    ta = [list(m["taxonomy"]) for m in a.metadata(axis="observation")]
    tb = [list(m["taxonomy"]) for m in b.metadata(axis="observation")]
    assert ta == tb, what + ": taxonomy"
    for k in ("BODY_SITE", "BarcodeSequence"):
        assert [m[k] for m in a.metadata(axis="sample")] == [m[k] for m in b.metadata(axis="sample")], what + ": " + k
    print("ok:", what)

for name in ("rich_from_julia_hdf5.biom", "rich_from_julia_json.biom"):
    same(biom.load_table(os.path.join(out, name)), ref, name)

with biom_open(os.path.join(out, "group_from_julia_hdf5.biom")) as h:
    g = Table.from_hdf5(h)          # not load_table, whose JSON fallback hides the HDF5 error
gm = g.group_metadata(axis="observation")
tree = gm["tree"] if isinstance(gm["tree"], str) else gm["tree"][1]      # biom returns the value alone
assert tree.startswith("(GG_OTU_5:0.1,") and tree.strip().endswith(";"), gm
print("ok: group_from_julia_hdf5.biom (newick tree)")

f = biom.load_table(os.path.join(out, "float_from_julia_hdf5.biom"))
assert np.allclose(f.matrix_data.toarray(), [[0.5, 0.0, 2.0], [1.25, 2.0, 0.0]])
assert [m["depth"] for m in f.metadata(axis="sample")] == [10.0, 2.5, 4.0]
assert [m["site"] for m in f.metadata(axis="sample")] == ["x", "y", "z"]
print("ok: float_from_julia_hdf5.biom")

# for read.jl: the reference writing BIOM 2.1
with biom_open(os.path.join(out, "rich_from_python_hdf5.biom"), "w") as h:
    ref.to_hdf5(h, "biom-format reference")
print("written: rich_from_python_hdf5.biom")
