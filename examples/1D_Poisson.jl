#
using NodalPolynomialSpaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(NodalPolynomialSpaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using LinearAlgebra, LinearSolve, Plots

N = 32

space = GaussLobattoLegendreSpace(N)
discr = Galerkin()

(x,) = points(space)

op = laplaceOp(space, discr)
f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => DirichletBC(),
           :Upper1 => DirichletBC(),
          )

prob = BVPDEProblem(op, f, bcs, space, discr)
alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp_dd")
