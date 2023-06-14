#
Spaces.size(V::NodalPolynomialSpace) = V.npoints

Spaces.domain(V::NodalPolynomialSpace) = V.dom

Spaces.points(V::NodalPolynomialSpace) = V.grid

Spaces.quadratures(V::NodalPolynomialSpace) = V.quads

Spaces.global_numbering(V::NodalPolynomialSpace) = V.glo_num

function Spaces.boundary_nodes(V::NodalPolynomialSpace)
    D = ndims(V)
    npoints = size(V)
    glo_num = global_numbering(V)

    indices = ()
    for i=1:D
        n = npoints[i]
        range_lower = ([1:npoints[j] for j=1:i-1]..., 1, [1:npoints[j] for j=i+1:D]...)
        range_upper = ([1:npoints[j] for j=1:i-1]..., n, [1:npoints[j] for j=i+1:D]...)
        indices = (indices..., glo_num[range_lower...])
        indices = (indices..., glo_num[range_upper...])
    end

    indices
end
#
