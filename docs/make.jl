push!(LOAD_PATH, "./../src/")
using SDDKleinianGroups

using Documenter
makedocs(
    modules = [SDDKleinianGroups],
    sitename = "SDDKleinianGroups Reference",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        warn_outdated = true,
        collapselevel=1,
        )
)
