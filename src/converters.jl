#
function (::Type{T})(V::NodalPolynomialSpace) where{T<:Number}
    npoints = size(V)
    dom  = T(domain(V))

    quads = Tuple((T.(z), T.(w)) for (z, w) in quadratures(V))
    grid  = Tuple(T.(x) for x in points(V))

    mass_mat = T.(V.mass_mat)
    deriv_mats = Tuple(T.(D) for D in V.deriv_mats)

    glo_num = global_numbering(V)

    NodalPolynomialSpace(
                         npoints, dom, quads, grid,
                         mass_mat, deriv_mats, glo_num,
                        )
end
#
