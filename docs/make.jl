using Documenter
using Literate
using ShallowWaters

EXAMPLE = joinpath(@__DIR__, "assets", "swm_equations.jl")
OUTPUT = joinpath(@__DIR__, "src")

# Generate markdown
binder_badge = "# [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/milankl/ShallowWaters.jl/gh-pages?labpath=dev%2Fswm_equations.ipynb)"
function preprocess_docs(content)
  return string(binder_badge, "\n\n", content)
end

Literate.markdown(EXAMPLE, OUTPUT; preprocess=preprocess_docs, codefence="```julia" => "```")
Literate.notebook(EXAMPLE, OUTPUT)

pages = [
  "Introduction" => "index.md",
  "Tutorial" => "swm_equations.md",
  "Reference" => "reference.md",
]

makedocs(
  sitename = "ShallowWaters.jl",
  format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
  modules = [ShallowWaters],
  pages = pages,
)

deploydocs(
  repo = "github.com/milankl/ShallowWaters.jl.git",
  push_preview = true,
  devbranch = "main",
)
