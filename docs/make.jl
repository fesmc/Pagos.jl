# Standard stuff
cd(@__DIR__)
using Pkg
Pkg.activate(".")
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using CairoMakie, Documenter, Literate
using DocumenterTools: Themes
using DocumenterCitations
ENV["JULIA_DEBUG"] = "Documenter"

# Packages specific to these docs
using Pagos

"""
Build a `CitationBibliography`, tolerating a `.bib` entry that is missing a required
field (e.g. an `@article` with no `journal`) as a warning rather than a hard failure.

`DocumenterCitations.CitationBibliography` (checked against 1.4.1, the latest release)
hardcodes `check = :error` in its call to `Bibliography.import_bibtex` and does not
expose the `check` keyword, so there is no way to ask it for lenient parsing directly.
On the specific `BibInternal` "entry is missing the ... field(s)" error, this patches a
*scratch copy* of the bib text with an explicit `[missing]` placeholder for the
offending field (visibly a gap, not a fabricated fact — the real `.bib` file is never
touched) and retries, looping to cover more than one broken entry. Any other error
(including the "mutually exclusive required fields" message shape, which does not
patch cleanly into valid BibTeX syntax) is rethrown unchanged.
"""
function lenient_bibliography(bibfile; style, max_attempts = 30)
    text = read(bibfile, String)
    for _ in 1:max_attempts
        path = tempname() * ".bib"
        write(path, text)
        try
            return CitationBibliography(path; style)
        catch err
            err isa ErrorException || rethrow()
            m = match(r"^Entry (\S+) is missing the (.+) field\(s\)\.$", err.msg)
            (m === nothing || contains(m.captures[2], "≡")) && rethrow()
            id, fields_str = m.captures
            fields = split(fields_str, ", ")
            @warn "docs/make.jl: bib entry `$id` is missing $(join(fields, ", ")); " *
                  "filling with a `[missing]` placeholder so the docs build. Fix " *
                  "`docs/src/assets/pagos.bib` when convenient."
            placeholder = join(("  $f = {[missing]}" for f in fields), ",\n") * ","
            text = replace(text, Regex("(@\\w+\\{$id,)") => SubstitutionString("\\1\n" * placeholder); count = 1)
        end
    end
    error("lenient_bibliography: exceeded $max_attempts patch attempts on $bibfile")
end

bib = lenient_bibliography(
    joinpath(@__DIR__, "src/assets", "pagos.bib");
    style=:authoryear
)

Literate.markdown("src/examples/mismip.jl", "src/examples"; credit = false)

Literate.markdown("src/physics/material.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/basal_friction.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/calving.jl", "src/physics"; credit = false)
Literate.markdown("src/physics/topography.jl", "src/physics"; credit = false)

Literate.markdown("src/numerics/integrators.jl", "src/numerics"; credit = false)
Literate.markdown("src/numerics/staggered_grids.jl", "src/numerics"; credit = false)

# %% Build docs
PAGES = [
    "index.md",
    "Examples" => ["src/examples/mismip.md"],
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
    # "Guidelines" => [
    #     "guidelines/naming.md",
    #     "guidelines/performance.md",
    # ],
    "References" => ["API_public.md", "API_private.md", "references.md"],
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
        collapselevel = 1,
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