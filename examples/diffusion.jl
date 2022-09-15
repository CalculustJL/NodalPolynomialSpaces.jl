using AbstractPDEInterfaces, NodalPolynomialSpaces

N = 64
ν = 0.01

xdom = IntervalDomain(-1, 1, boundary_tags=(:left, :right))
ydom = IntervalDomain(-1, 1, boundary_tags=(:lower, :upper))
dom = xdom ⊗ ydom

space = GaussLobattoLegendreSpace(N, N; dom=dom)
(x, y,) = points(space)

bcs = Dict(
           :left => NeumannBC(),
           :right => DirichletBC(),
           :lower => DirichletBC(),
           :upper => NeumannBC(),
          )

discr = Galerkin()
lhsOp = diffusionOp(ν, space, discr)

f = @. 0*x + 1

prob = BoundaryValueProblem(lhsOp, f, bcs, space, discr,
                            abstol=1e-8, reltol=1e-8)

alg  = LinearBoundaryValueAlg(linalg=KrylovJL_CG())

@time sol = solve(prob, alg; verbose=false)
@test sol.resid < 1e-8
