# Tables.jl: the long form, one row per non-zero value, so `DataFrame(t)` and
# any other Tables.jl sink work. Zeros are left out, as in the file.

Tables.istable(::Type{<:BIOMTable}) = true
Tables.columnaccess(::Type{<:BIOMTable}) = true

function Tables.columns(t::BIOMTable)
    I, J, V = findnz(t.data)
    (observation = t.observations[I], sample = t.samples[J], value = V)
end
