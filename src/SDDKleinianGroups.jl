
"""
Drawings and calculations useful for research in the theory of Kleinian groups
and other groups of complex functions on the Riemann sphere.
"""
module SDDKleinianGroups

using Reexport

@reexport using SDDCore, SDDGeometry, SDDGraphics


include("traversegroup.jl")


include("orbits.jl")

export
    drawpointorbit,
    drawpointssetorbit,
    drawcircleorbit,
    drawcirclessetorbit


include("limitsets.jl")

export
    drawlimitset, drawΛ,
    drawfixedpointslimitset, drawfixedpointsΛ,
    drawchaosgamelimitset, drawchaosgameΛ,
#=    drawdensitylimitset, drawdensityΛ, =#
    drawtrappedpointslimitset, drawtrappedpointsΛ


#include("connectedlimitset.jl")

#=export
    drawconnectedlimitset, drawconnectedΛ,
=#


#include("fundamentalregions.jl")

#=export
=#

end
