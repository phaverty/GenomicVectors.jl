using Documenter, GenomicVectors

makedocs()

deploydocs(
           repo = "github.com/phaverty/GenomicVectors.jl.git",
           julia = "release",
           DOCUMENTER_DEBUG = true
           )
