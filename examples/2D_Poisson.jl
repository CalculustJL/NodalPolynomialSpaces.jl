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
           :Lower1 => NeumannBC(),
           :Upper1 => DirichletBC(),

           :Lower2 => DirichletBC(),
           :Upper2 => NeumannBC(),
          )

prob = BoundaryValueProblem(op, f, bcs, space, discr, abstol=1e-8, reltol=1e-8)
alg  = LinearBoundaryValueAlg(linalg=KrylovJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
plt = plot(sol)
savefig(plt, "bvp2d_dn")
