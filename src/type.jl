#
"""
Lagrange polynomial spectral space
"""
struct NodalPolynomialSpace{T,
                            D,
                            Tpts<:NTuple{D},
                            Tdom<:Domains.AbstractDomain{T,D},
                            Tquad,
                            Tgrid,
                            Tmass,
                            Tderiv,
                            Tglo
                           } <: Spaces.AbstractSpace{T,D}
    """ size """
    npoints::Tpts
    """ Domain """
    dom::Tdom
    """ quadratures """
    quads::Tquad
    """ grid points """
    grid::Tgrid
    """ mass matrix """
    mass_mat::Tmass
    """ derivative matrices """
    deriv_mats::Tderiv
    """ global numbering """
    glo_num::Tglo
end

function NodalPolynomialSpace(n::Integer;
        domain::AbstractDomain{<:Any,1} = ChebyshevDomain(1),
        quadrature = gausschebyshev,
        T = Float64,
       )

    if !isa(domain, Domains.LogicallyRectangularDomain)
        msg = """Spectral polynomials work with logically rectangular
            domains. `domain` must either be a `Domains.IntervalDomain`,
            or product of interval domains created with `LinearAlgebra.×`.
            Optionally `domain` may be a `Domains.MappedDomain` generated
            as `Domains.deform(dom, mapping)`.
            """
        throw(ArgumentError(msg))
    end

    # check for deformation
    dom, mapping = if domain isa Domains.MappedDomain
        domain.domain, domain.mapping
    else
        domain, nothing
    end

    # put domain in a box
    if dom isa IntervalDomain
        dom = ProductDomain(dom)
    end

    @assert dom isa Domains.BoxedDomain

    # JUST MULTIPLY DIFF MATRIX AND MASS MATRIX with the right scaling

    #""" reset deformation to map from [-1,1]^D """
    #ref_dom = ChebyshevDomain(1)
    #dom = ref_domain # map_from_ref(domain, ref_dom) # TODO
    ## change dom eltype

    z, w = quadrature(n)

    z = T.(z)
    w = T.(w)

    D = lagrange_deriv_mat(z)

    npoints = (n,)
    dom     = T(dom)
    quads   = ((z, w),)
    grid    = vec.((z,))
    mass_mat = vec(w)
    deriv_mats = (D,)
    glo_num = reshape(1:prod(npoints), npoints)

    V = NodalPolynomialSpace(
                             npoints, dom, quads, grid,
                             mass_mat, deriv_mats, 
                             glo_num,
                            )

    isnothing(mapping) ? V : deform(V, mapping)
end

function NodalPolynomialSpace(nr::Integer, ns::Integer;
        domain::AbstractDomain{<:Number,2} = ChebyshevDomain(2),
        quadrature = gausslobatto,
        T = Float64,
       )

    if !isa(domain, Domains.LogicallyRectangularDomain)
        msg = """Trigonometric polynomials work with logically rectangular
            domains. `domain` must be a product of `Domains.IntervalDomain`
            created with `LinearAlgebra.×`. Optionally `domain` may be a
            `Domains.MappedDomain` generated as `Domains.deform(dom, map)`.
            """
        throw(ArgumentError(msg))
    end

    # check for deformation
    dom, mapping = if domain isa Domains.MappedDomain
        domain.domain, domain.mapping
    else
        domain, nothing
    end

    # put domain in a box
    @assert dom isa Domains.BoxedDomain

    #""" reset deformation to map from [-1,1]^D """
    #ref_dom = ChebyshevDomain(2)
    #dom = ref_dom # map_from_ref(dom, ref_dom) # TODO

    zr, wr = quadrature(nr)
    zs, ws = quadrature(ns)

    zr, wr = T.(zr), T.(wr)
    zs, ws = T.(zs), T.(ws)

    r, s = ndgrid(zr, zs)

    Dr = lagrange_deriv_mat(zr)
    Ds = lagrange_deriv_mat(zs)

    npoints = (nr, ns,)
    dom = T(dom)
    quads = ((zr, wr), (zs, ws),)
    grid = vec.((r, s,))
    mass_mat = vec(wr * ws')
    deriv_mats = (Dr, Ds,)
    glo_num = reshape(1:prod(npoints), npoints)

    V = NodalPolynomialSpace(
                             npoints, dom, quads, grid,
                             mass_mat, deriv_mats,
                             glo_num,
                            )

    isnothing(mapping) ? V : deform(V, mapping)
end

GaussLobattoLegendreSpace(args...; kwargs...) = NodalPolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendreSpace(args...; kwargs...) = NodalPolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebyshevSpace(args...; kwargs...) = NodalPolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)

const GLLSpace = GaussLobattoLegendreSpace
const GLSpace = GaussLegendreSpace
const GCSpace = GaussChebyshevSpace
#
