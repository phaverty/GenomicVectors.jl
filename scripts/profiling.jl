using GenomicRanges
using JLD
using Base.Test

ds = load("/Users/phaverty/.julia/v0.4/GenomicRanges/test/testData/pinfo.jld")

positions = ds["pos"]
chromosomes = ds["chr"]
info = ds["info"]

gp = genopos( positions, chromosomes, info)
x = GenomicPositions( positions, chromosomes, info )


#genopos( x );
#chr( x );
#chrpos( x );
#chrpos2( x );
#chrpos4( x )
using ProfileView
Profile.clear(); @profile for i in 1:1e2 gp = genopos( positions, chromosomes, info); end; ProfileView.view()
Profile.clear(); @profile for i in 1:1e2 gp = chrpos(x); end; ProfileView.view()
Profile.clear(); @profile for i in 1:1e2 gp = chr(x); end; ProfileView.view()



#@time genopos( positions, chromosomes, info);
#@time chrpos( x );
#@time chrpos2( x );
#@time chrpos4( x );
#@time chr( x )
#

#@code_warntype chrpos2( x );
