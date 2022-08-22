#
###
# vector calculus ops
###

function massOp(space::NodalPolynomialSpace, ::Galerkin)
    mass_mat = mass_matrix(space)

    DiagonalOperator(mass_mat)
end

function gradientOp(space::NodalPolynomialSpace{<:Number,1})
    (Dr,) = space.deriv_mats

    Dx = MatrixOperator(Dr)

    DD = [Dx,]
end

function gradientOp(space::NodalPolynomialSpace{<:Number,2})
    (nr, ns) = space.npoints
    (Dr, Ds) = space.deriv_mats

    Ir = IdentityOperator{nr}()
    Is = IdentityOperator{ns}()

    Dx = ⊗(Is, Dr)
    Dy = ⊗(Ds, Ir)

    DD = [Dx
          Dy]
end

function gradientOp(space::NodalPolynomialSpace{<:Number,3})
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

function interpOp(space1::NodalPolynomialSpace{<:Number,1},
                  space2::NodalPolynomialSpace{<:Number,1})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    J = lagrange_interp_mat(r2, r1) # from 1 to 2

    MatrixOperator(J)
end

function interpOp(space1::NodalPolynomialSpace{<:Number,2},
                  space2::NodalPolynomialSpace{<:Number,2})
    r1, _ = space1.quads[1]
    r2, _ = space2.quads[1]

    s1, _ = space1.quads[2]
    s2, _ = space2.quads[2]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)

    ⊗(Js, Jr)
end

function interpOp(space1::NodalPolynomialSpace{<:Number,3},
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
