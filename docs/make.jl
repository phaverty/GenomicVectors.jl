using Documenter, GenomicVectors

makedocs()

deploydocs(
           repo = "github.com/phaverty/GenomicVectors.jl.git",
           julia = "0.5",
           DOCUMENTER_DEBUG = true
           )
