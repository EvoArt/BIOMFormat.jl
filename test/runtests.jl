using BIOMFormat
using BIOMFormat: counts, samples, taxonomy
using SparseArrays
using Tables
using Test
import DenoiseDNA

# The example tables of the reference implementation (biocore/biom-format,
# BSD licence: fixtures/LICENSE-biom-format.txt), each in both formats.
fixture(name) = joinpath(@__DIR__, "fixtures", name)

@testset "BIOMFormat" begin

@testset "BIOM 1.0 (JSON)" begin
    t = read_biom(fixture("min_sparse_otu_table.biom"))
    @test size(t) == (5, 6)
    @test observations(t) == ["GG_OTU_$k" for k in 1:5]
    @test samples(t) == ["Sample$k" for k in 1:6]
    @test eltype(counts(t)) == Int
    @test counts(t)[1, 3] == 1 && counts(t)[2, 1] == 5 && counts(t)[3, 4] == 4 && counts(t)[1, 1] == 0
    @test sum(counts(t)) == 27
    @test nnz(counts(t)) == 15
    @test t.type == "OTU table"
    @test t.generated_by == "QIIME revision XYZ"
    @test t.observation_metadata === nothing

    r = read_biom(fixture("rich_sparse_otu_table.biom"))
    @test taxonomy(r)[1] == ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
                             "o__Enterobacteriales", "f__Enterobacteriaceae", "g__Escherichia", "s__"]
    @test "taxonomy" in metadata(r, :observation)
    @test "BODY_SITE" in metadata(r, :sample)
    @test sample_metadata(r, "BODY_SITE")[1] == "gut"
    @test read_biom(IOBuffer(read(fixture("rich_sparse_otu_table.biom")))) == r
end

@testset "BIOM 2.1 (HDF5) agrees with BIOM 1.0" begin
    for name in ("min_sparse_otu_table", "rich_sparse_otu_table")
        j = read_biom(fixture("$name.biom"))
        h = read_biom(fixture("$(name)_hdf5.biom"))
        @test observations(h) == observations(j)
        @test samples(h) == samples(j)
        @test counts(h) == counts(j)
        @test eltype(counts(h)) == Int                 # stored as float64, narrowed
        @test taxonomy(h) == taxonomy(j)
        for key in metadata(j, :sample)
            @test sample_metadata(h, key) == sample_metadata(j, key)
        end
    end
    g = read_biom(fixture("rich_sparse_otu_table_hdf5_group_metadata.biom"))
    gm = group_metadata(g, :observation)
    @test haskey(gm, "tree")
    @test first(gm["tree"]) == "newick"
    @test startswith(last(gm["tree"]), "(GG_OTU_5:0.1,") && endswith(strip(last(gm["tree"])), ";")
end

@testset "round trips" begin
    dir = mktempdir()
    for name in ("min_sparse_otu_table.biom", "rich_sparse_otu_table.biom",
                 "rich_sparse_otu_table_hdf5_group_metadata.biom")
        t = read_biom(fixture(name))
        @test read_biom(write_biom(joinpath(dir, "a.biom"), t)) == t
        if isempty(group_metadata(t, :observation)) && isempty(group_metadata(t, :sample))
            @test read_biom(write_biom(joinpath(dir, "b.biom"), t; format = :json)) == t
            @test read_biom(seekstart(write_biom(IOBuffer(), t))) == t
        else
            @test_logs (:warn, r"no group metadata") write_biom(joinpath(dir, "b.biom"), t; format = :json)
        end
    end
    # non-integer values stay Float64, numeric metadata survives HDF5
    f = BIOMTable([0.5 0.0; 1.25 2.0], ["a", "b"], ["s1", "s2"];
                  sample_metadata = [Dict("depth" => 10.0), Dict("depth" => 2.5)])
    @test eltype(counts(f)) == Float64
    g = read_biom(write_biom(joinpath(dir, "f.biom"), f))
    @test counts(g) == counts(f)
    @test sample_metadata(g, "depth") == [10.0, 2.5]
    @test_throws ArgumentError write_biom(joinpath(dir, "x.biom"), f; format = :csv)
end

@testset "construction, subsetting, Tables" begin
    t = read_biom(fixture("rich_sparse_otu_table.biom"))
    s = t[["GG_OTU_2", "GG_OTU_4"], 2:4]
    @test size(s) == (2, 3)
    @test observations(s) == ["GG_OTU_2", "GG_OTU_4"]
    @test counts(s) == counts(t)[[2, 4], 2:4]
    @test taxonomy(s) == taxonomy(t)[[2, 4]]
    @test samples(t[:, sample_metadata(t, "BODY_SITE") .== "gut"]) ==
          samples(t)[sample_metadata(t, "BODY_SITE") .== "gut"]
    @test_throws KeyError t["nope", :]

    cols = Tables.columntable(t)
    @test length(cols.value) == nnz(counts(t))
    @test sum(cols.value) == sum(counts(t))
    @test Tables.istable(t)

    @test_throws ArgumentError BIOMTable([1 2; 3 4], ["a", "a"], ["s", "t"])
    @test_throws DimensionMismatch BIOMTable([1 2; 3 4], ["a"], ["s", "t"])
    @test_throws DimensionMismatch BIOMTable([1 2; 3 4], ["a", "b"], ["s", "t"];
                                             sample_metadata = [Dict("x" => 1)])
    @test_throws ArgumentError metadata(t, :row)
    @test taxonomy(BIOMTable([1;;], ["a"], ["s"]; observation_metadata = [Dict("taxonomy" => "k__A; p__B")])) ==
          [["k__A", "p__B"]]
end

@testset "DenoiseDNA" begin
    seqs = ["ACGT"^10, "TTGCA"^8, "GGGACC"^7]
    st = DenoiseDNA.SequenceTable(["s1", "s2"], seqs, [10 0 3; 4 7 0])
    t = BIOMTable(st)
    @test size(t) == (3, 2)
    @test observations(t) == seqs
    @test counts(t) == [10 4; 0 7; 3 0]
    back = DenoiseDNA.SequenceTable(t)
    @test back.samples == st.samples && back.sequences == st.sequences && back.table == st.table

    a = BIOMTable(st; ids = :asv, id = "run1")
    @test observations(a) == ["ASV1", "ASV2", "ASV3"]
    @test observation_metadata(a, "sequence") == seqs
    @test a.id == "run1"
    h = read_biom(write_biom(joinpath(mktempdir(), "st.biom"), a))
    @test DenoiseDNA.SequenceTable(h).sequences == seqs
    @test_throws ArgumentError BIOMTable(st; ids = :hash)
    @test_throws ArgumentError DenoiseDNA.SequenceTable(BIOMTable([0.5;;], ["a"], ["s"]))
end

@testset "display" begin
    t = read_biom(fixture("rich_sparse_otu_table.biom"))
    @test sprint(show, t) == "BIOMTable(6 samples × 5 observations)"
    plain = sprint(show, MIME"text/plain"(), t)
    @test occursin("BIOMTable · 6 samples × 5 observations", plain)
    @test occursin("Most abundant observations", plain)
    @test occursin("g__Escherichia", plain)
    @test !occursin("\e[", plain)
    colour = sprint(show, MIME"text/plain"(), t; context = :color => true)
    @test occursin("\e[", colour)
    empty = BIOMTable(spzeros(Int, 0, 2), String[], ["a", "b"])
    @test occursin("2 samples × 0 observations", sprint(show, MIME"text/plain"(), empty))
    st = DenoiseDNA.SequenceTable(["s1"], ["ACGT"^10], fill(5, 1, 1))
    @test occursin("ACGT", sprint(show, MIME"text/plain"(), BIOMTable(st)))
end

end
