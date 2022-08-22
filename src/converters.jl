#
function (::Type{T})(space::NodalPolynomialSpace) where{T<:Number}
    npoints = size(space)
    dom  = T(domain(space))

    quads = Tuple((T.(z), T.(w)) for (z, w) in quadratures(space))
    grid  = Tuple(T.(x) for x in points(space))

    mass_mat = T.(mass_matrix(space))
    deriv_mats = Tuple(T.(D) for D in space.deriv_mats)

    loc_num = local_numbering(space)

    NodalPolynomialSpace(
                         npoints, dom, quads, grid,
                         mass_mat, deriv_mats, loc_num,
                        )
end

function Adapt.adapt_structure(to, space::NodalPolynomialSpace)
    grid  = adapt_structure(to, points(space))
    quads = adapt_structure(to, quadratures(space))
    mass_mat = adapt_structure(to, mass_matrix(space))
    loc_num = adapt_structure(to, local_numbering(space))

    x = first(grid)
    T = eltype(x)

    npoints = size(space)
    dom     = T(domain(space))

    NodalPolynomialSpace(
                         npoints, dom, quads, grid,
                         mass_mat, deriv_mats, loc_num,
                        )
end
#
