#
using NodalPolynomialSpaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(NodalPolynomialSpaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using LinearAlgebra, LinearSolve
using Plots, Test

N = 256
ν = 1e-4

space = GaussLobattoLegendreSpace(N)
discr = Galerkin()
(x,) = points(space)

v = @. 0*x + 1
f = @. 0*x + 1

A = diffusionOp(ν, space, discr)
C = advectionOp((v,), space, discr)

op = cache_operator(C-A, x)

bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => DirichletBC(),
          )

prob = BoundaryValueProblem(op, f, bcs, space, discr, abstol=1e-8, reltol=1e-8)
alg  = LinearBoundaryValueAlg(linalg=KrylovJL_CG())

@time sol = solve(prob, alg; verbose=false)
@show sol.resid
plt = plot(sol)
#
