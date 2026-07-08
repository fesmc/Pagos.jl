# Standard stuff
cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using CairoMakie, Documenter, Literate
using DocumenterTools: Themes
using DocumenterCitations
ENV["JULIA_DEBUG"] = "Documenter"

# Packages specific to these docs
using Pagos

bib = CitationBibliography(
    joinpath(@__DIR__, "src/assets", "pagos.bib");
    style=:authoryear
)

Literate.markdown("src/physics/material.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/basal_friction.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/calving.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/topography.jl", "src/physics"; credit = false)
Literate.markdown("src/numerics/integrators.jl", "src/numerics"; credit = false)
Literate.markdown("src/numerics/staggered_grids.jl", "src/numerics"; credit = false)

# %% Build docs
PAGES = [
    "index.md",
    "Quickstart" => [],
    "Physics" => [
        "physics/material.md",
        "physics/basal_friction.md",
        "physics/calving.md",
        "physics/topography.md",
    ],
    "Numerics" => [
        "numerics/integrators.md",
        "numerics/staggered_grids.md",
    ],
    "Guidelines" => [
        "guidelines/naming.md",
        "guidelines/performance.md",
    ],
    "API reference" => ["API_public.md", "API_private.md"],
    "references.md",
    # "examples/tutorial.md",
    # "Examples" => example_pages,
    # "References" => ref_pages,
]

include("style.jl")

makedocs(
    modules = [Pagos],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        collapselevel = 2,
        ),
    sitename = "Pagos.jl",
    authors = "Jan Swierczek-Jereczek, Alexander Robinson",
    pages = PAGES,
    doctest = CI,
    draft = false,
    plugins = [bib],
    checkdocs = :none,
    warnonly = true,
)

deploydocs(;
    repo="https://github.com/fesmc/Pagos.jl",
)