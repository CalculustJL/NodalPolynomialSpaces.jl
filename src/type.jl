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
                            Tloc,
#                           Tglo
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
    """ local numbering """
    loc_num::Tloc
#   """ global numbering """
#   glo_num::Tglo
end

function NodalPolynomialSpace(n::Integer;
        dom::Domains.AbstractDomain{<:Any,1}=ChebychevDomain(1),
        quadrature = gausschebychev,
        T = Float64,
       )

    if dom isa IntervalDomain
        dom = BoxDomain(dom)
    elseif !(dom isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_dom = ChebychevDomain(1)
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
    loc_num = reshape(1:prod(npoints), npoints)

    space = NodalPolynomialSpace(
                                 npoints, dom, quads, grid,
                                 mass_mat, deriv_mats, 
                                 loc_num,
                                )

    dom isa Domains.DeformedDomain ? deform(space, mapping) : space
end

function NodalPolynomialSpace(nr::Integer, ns::Integer;
        dom::Domains.AbstractDomain{<:Number,2}=ChebychevDomain(2),
        quadrature = gausslobatto,
        T = Float64,
       )

    if !(dom isa BoxDomain)
        @error "spectral polynomials work with logically rectangular domains"
    end

    #""" reset deformation to map from [-1,1]^D """
    #ref_dom = ChebychevDomain(2)
    #dom = ref_dom # map_from_ref(dom, ref_dom) # TODO

    zr, wr = quadrature(nr)
    zs, ws = quadrature(ns)

    zr, wr = T.(zr), T.(wr)
    zs, ws = T.(zs), T.(ws)

    r, s = ndgrid(zr,zs)

    Dr = lagrange_deriv_mat(zr)
    Ds = lagrange_deriv_mat(zs)

    npoints = (nr, ns,)
    dom = T(dom)
    quads = ((zr, wr), (zs, ws),)
    grid = vec.((r, s,))
    mass_mat = vec(wr * ws')
    deriv_mats = (Dr, Ds,)
    loc_num = reshape(1:prod(npoints), npoints)

    space = NodalPolynomialSpace(
                                 npoints, dom, quads, grid,
                                 mass_mat, deriv_mats,
                                 loc_num,
                                )

    dom isa Domains.DeformedDomain ? deform(space, mapping) : space
end

GaussLobattoLegendreSpace(args...; kwargs...) =
    NodalPolynomialSpace(args...; quadrature=gausslobatto, kwargs...)
GaussLegendreSpace(args...; kwargs...) =
    NodalPolynomialSpace(args...; quadrature=gausslegendre, kwargs...)
GaussChebychevSpace(args...; kwargs...) =
    NodalPolynomialSpace(args...; quadrature=gausschebyshev, kwargs...)
#
