#
Spaces.size(space::NodalPolynomialSpace) = space.npoints

Spaces.domain(space::NodalPolynomialSpace) = space.dom

Spaces.points(space::NodalPolynomialSpace) = space.grid

Spaces.quadratures(space::NodalPolynomialSpace) = space.quads

Spaces.mass_matrix(space::NodalPolynomialSpace) = space.mass_mat

Spaces.local_numbering(space::NodalPolynomialSpace) = space.loc_num

function Spaces.global_numbering(space::NodalPolynomialSpace)
    dom = domain(space)
    loc_num = local_numbering(space)
end

function Spaces.boundary_nodes(space::NodalPolynomialSpace)
    D = dims(space)
    npoints = size(space)
    loc_num = local_numbering(space)

    indices = ()
    for i=1:D
        n = npoints[i]
        range_lower = ([1:npoints[j] for j=1:i-1]...,1,[1:npoints[j] for j=i+1:D]...)
        range_upper = ([1:npoints[j] for j=1:i-1]...,n,[1:npoints[j] for j=i+1:D]...)
        indices = (indices..., loc_num[range_lower...])
        indices = (indices..., loc_num[range_upper...])
    end

    indices
end
#
