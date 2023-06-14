using NodalPolynomialSpaces
let
    # add dependencies to env stack
    pkgpath = dirname(dirname(pathof(NodalPolynomialSpaces)))
    tstpath = joinpath(pkgpath, "test")
    !(tstpath in LOAD_PATH) && push!(LOAD_PATH, tstpath)
    nothing
end

using LinearAlgebra, LinearSolve
using Test, Plots

N = 64
ν = 0.1

space = GaussLobattoLegendreSpace(N, N)
discr = Galerkin()

op = diffusionOp(ν, space, discr)
(x, y,) = points(space)

f = @. sin(x) * cos(y)
bcs = Dict(
           :D1_inf => NeumannBC(),
           :D1_sup => DirichletBC(),

           :D2_inf => DirichletBC(),
           :D2_sup => NeumannBC(),
          )

prob = BoundaryValueProblem(op, f, bcs, space, discr, abstol=1e-8, reltol=1e-8)
alg  = LinearBoundaryValueAlg(linalg=KrylovJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp2d_dn")
