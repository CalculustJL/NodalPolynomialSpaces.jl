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

N = 64

space = GaussLobattoLegendreSpace(N)
discr = Galerkin()
(x,)  = points(space)

op = -laplaceOp(space, discr)
f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => DirichletBC(),
          )

prob = BoundaryValueProblem(op, f, bcs, space, discr, abstol=1e-8, reltol=1e-8)
alg  = LinearBoundaryValueAlg(linalg=KrylovJL_CG())

@time sol = solve(prob, alg; verbose=false)
plt = plot(sol)
display(plt)
#savefig(plt, "bvp_dd")
@test sol.resid < 1e-6
#
