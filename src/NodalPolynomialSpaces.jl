module NodalPolynomialSpaces

using Reexport
@reexport using CalculustCore
using CalculustCore.Spaces: AbstractSpace, AbstractDiscretization, TransformedSpace
using CalculustCore.Domains: AbstractDomain, ProductDomain

using SciMLOperators: AbstractSciMLOperator

using LinearAlgebra
using SparseArrays

import FastGaussQuadrature: gausslobatto, gausslegendre, gausschebyshev

include("lagrange_matrices.jl")
include("type.jl")
include("converters.jl")
include("interface.jl")
include("operators.jl")

export NodalPolynomialSpace,
       GaussLobattoLegendreSpace, GaussLegendreSpace, GaussChebyshevSpace,
       GLLSpace, GLSpace, GCSpace

end # module
