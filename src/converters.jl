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
#
