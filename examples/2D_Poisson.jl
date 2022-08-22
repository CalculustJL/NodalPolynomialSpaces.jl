using NodalPolynomialSpaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(NodalPolynomialSpaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using LinearAlgebra, LinearSolve, OrdinaryDiffEq, Plots

N = 32

space = GaussLobattoLegendreSpace(N, N)
discr = Galerkin()

op = laplaceOp(space, discr)
(x, y,) = points(space)

f  = @. 0*x + 1
bcs = Dict(
           :Lower1 => NeumannBC(),
           :Upper1 => DirichletBC(),

           :Lower2 => DirichletBC(),
           :Upper2 => NeumannBC(),
          )

prob = BVPDEProblem(op, f, bcs, space, discr)
alg  = LinearBVPDEAlg(linalg=IterativeSolversJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp2d_dn")
