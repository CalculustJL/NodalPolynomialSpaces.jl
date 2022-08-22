module NodalPolynomialSpaces

using Reexport
@reexport using AbstractPDEInterfaces

using LinearAlgebra
using SparseArrays
using Adapt

import FastGaussQuadrature: gausslobatto, gausslegendre, gausschebyshev

include("lagrange_matrices.jl")
include("type.jl")
include("converters.jl")
include("interface.jl")
include("operators.jl")

export NodalPolynomialSpace,
       GaussLobattoLegendreSpace, GaussLegendreSpace, GaussChebychevSpace

end # module
