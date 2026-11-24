using Documenter, BIOMFormat

makedocs(
    sitename = "BIOMFormat.jl",
    modules = [BIOMFormat],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true",
                             size_threshold = nothing),
    pages = ["Home" => "index.md", "API" => "api.md"],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(repo = "github.com/EvoArt/BIOMFormat.jl.git", devbranch = "master")
