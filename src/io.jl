# Reading and writing, with the format worked out from the file.

const _HDF5_MAGIC = UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]

_is_hdf5(path) = open(io -> read(io, 8) == _HDF5_MAGIC, path)

"""
    read_biom(path) -> BIOMTable
    read_biom(io::IO) -> BIOMTable

Read a BIOM table. Both BIOM 2.1 (HDF5, what QIIME 2 exports) and BIOM 1.0
(JSON) are read from a path, the format told from the file's first bytes; from
an `IO`, only JSON.
"""
read_biom(path::AbstractString) = _is_hdf5(path) ? _read_hdf5(path) : open(_read_json, path)
read_biom(io::IO) = _read_json(io)

"""
    write_biom(path, t::BIOMTable; format = :hdf5)
    write_biom(io::IO, t::BIOMTable)

Write `t` as BIOM 2.1 (`format = :hdf5`, the default) or BIOM 1.0
(`format = :json`). An `IO` gets JSON. BIOM 1.0 has no group metadata, so a
table with some warns that it is dropped.
"""
function write_biom(path::AbstractString, t::BIOMTable; format::Symbol = :hdf5)
    format === :hdf5 ? _write_hdf5(path, t) :
    format === :json ? open(io -> _write_json(io, t), path, "w") :
    throw(ArgumentError("format is :hdf5 or :json, got :$format"))
    path
end
write_biom(io::IO, t::BIOMTable) = (_write_json(io, t); io)
