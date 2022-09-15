#
Spaces.size(space::NodalPolynomialSpace) = space.npoints

Spaces.domain(space::NodalPolynomialSpace) = space.dom

Spaces.points(space::NodalPolynomialSpace) = space.grid

Spaces.quadratures(space::NodalPolynomialSpace) = space.quads

Spaces.mass_matrix(space::NodalPolynomialSpace) = space.mass_mat

Spaces.global_numbering(space::NodalPolynomialSpace) = space.glo_num

function Spaces.boundary_nodes(space::NodalPolynomialSpace)
    D = dims(space)
    npoints = size(space)
    glo_num = global_numbering(space)

    indices = ()
    for i=1:D
        n = npoints[i]
        range_lower = ([1:npoints[j] for j=1:i-1]...,1,[1:npoints[j] for j=i+1:D]...)
        range_upper = ([1:npoints[j] for j=1:i-1]...,n,[1:npoints[j] for j=i+1:D]...)
        indices = (indices..., glo_num[range_lower...])
        indices = (indices..., glo_num[range_upper...])
    end

    indices
end
#
