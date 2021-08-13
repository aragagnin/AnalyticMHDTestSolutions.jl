using Documenter
using GadgetIO

makedocs(
    sitename = "AnalyticMHDTestSolutions.jl",
    format = Documenter.HTML(),
    modules = [GadgetIO],
    pages = [
            "Table of Contents" => "index.md",
            "Install"           => "install.md",
            "Sod Shocks"        => "sod_shocks.md",
            "Sedov-Taylor Blastwave"   => "sedov.md",
            "Contributing"      => "contribute.md",
            "API reference"     => "api.md"
            ]
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/AnalyticMHDTestSolutions.jl.git"
)
