#
###
# vector calculus ops
###

# TODO - discretization option
function Spaces.massOp(space::NodalPolynomialSpace, ::Galerkin)
    mass_mat = mass_matrix(space)

    DiagonalOperator(mass_mat)
end

function Spaces.gradientOp(space::NodalPolynomialSpace{<:Number,1},
                           ::Spaces.AbstractDiscretization)
    (Dr,) = space.deriv_mats

    Dx = MatrixOperator(Dr)

    DD = [Dx,]
end

function Spaces.gradientOp(space::NodalPolynomialSpace{<:Number,2},
                          ::Spaces.AbstractDiscretization)
    (nr, ns) = space.npoints
    (Dr, Ds) = space.deriv_mats

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()

    Dx = ⊗(Is, Dr)
    Dy = ⊗(Ds, Ir)

    DD = [Dx
          Dy]
end

function Spaces.gradientOp(space::NodalPolynomialSpace{<:Number,3},
                           ::Spaces.AbstractDiscretization)
    (Dr, Ds, Dt) = space.deriv_mats
    (nr, ns, nt) = space.npoints

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()
    It = IdentityOperator{nt}()

    Dx = ⊗(It, Is, Dr)
    Dy = ⊗(It, Ds, It)
    Dz = ⊗(Dt, Is, It)

    DD = [Dx
          Dy
          Dz]
end

###
# interpolation operators
###

function Spaces.interpOp(space1::NodalPolynomialSpace{<:Number,1},
                         space2::NodalPolynomialSpace{<:Number,1})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    J = lagrange_interp_mat(r2, r1) # from 1 to 2

    MatrixOperator(J)
end

function Spaces.interpOp(space1::NodalPolynomialSpace{<:Number,2},
                         space2::NodalPolynomialSpace{<:Number,2})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    s1, _ = space1.quads[2]
    s2, _ = space2.quads[2]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)

    ⊗(Js, Jr)
end

function Spaces.interpOp(space1::NodalPolynomialSpace{<:Number,3},
                         space2::NodalPolynomialSpace{<:Number,3})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    s1, _ = space1.quads[2]
    s2, _ = space2.quads[2]

    t1, _ = space1.quads[3]
    t2, _ = space2.quads[3]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)
    Jt = lagrange_interp_mat(t2, t1)

    ⊗(Jt, Js, Jr)
end
#
