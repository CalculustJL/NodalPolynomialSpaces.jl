#
###
# vector calculus ops
###

function Spaces.massOp(V::NodalPolynomialSpace, ::Galerkin)

    DiagonalOperator(V.mass_mat)
end

function Spaces.gradientOp(V::NodalPolynomialSpace{<:Number,1},
                           ::Spaces.AbstractDiscretization)
    (Dr,) = V.deriv_mats

    Dx = MatrixOperator(Dr)

    DD = AbstractSciMLOperator[Dx,]
end

function Spaces.gradientOp(V::NodalPolynomialSpace{<:Number,2},
                          ::Spaces.AbstractDiscretization)
    (nr, ns) = V.npoints
    (Dr, Ds) = V.deriv_mats

    Ir = IdentityOperator(nr)
    Is = IdentityOperator(ns)

    Dx = ⊗(Is, Dr)
    Dy = ⊗(Ds, Ir)

    DD = AbstractSciMLOperator[Dx
                               Dy]
end

function Spaces.gradientOp(V::NodalPolynomialSpace{<:Number,3},
                           ::Spaces.AbstractDiscretization)
    (Dr, Ds, Dt) = V.deriv_mats
    (nr, ns, nt) = V.npoints

    Ir = IdentityOperator(nr)
    Is = IdentityOperator(ns)
    It = IdentityOperator(nt)

    Dx = ⊗(It, Is, Dr)
    Dy = ⊗(It, Ds, It)
    Dz = ⊗(Dt, Is, It)

    DD = AbstractSciMLOperator[Dx
                               Dy
                               Dz]
end

###
# interpolation operators
###

function Spaces.interpOp(V1::NodalPolynomialSpace{<:Number,1},
                         V2::NodalPolynomialSpace{<:Number,1})
    r1, _ = V1.quads[1]
    r2, _ = V2.quads[1]

    J = lagrange_interp_mat(r2, r1) # from 1 to 2

    MatrixOperator(J)
end

function Spaces.interpOp(V1::NodalPolynomialSpace{<:Number,2},
                         V2::NodalPolynomialSpace{<:Number,2})
    r1, _ = V1.quads[1]
    r2, _ = V2.quads[1]

    s1, _ = V1.quads[2]
    s2, _ = V2.quads[2]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)

    ⊗(Js, Jr)
end

function Spaces.interpOp(V1::NodalPolynomialSpace{<:Number,3},
                         V2::NodalPolynomialSpace{<:Number,3})
    r1, _ = V1.quads[1]
    r2, _ = V2.quads[1]

    s1, _ = V1.quads[2]
    s2, _ = V2.quads[2]

    t1, _ = V1.quads[3]
    t2, _ = V2.quads[3]

    Jr = lagrange_interp_mat(r2, r1) # from 1 to 2
    Js = lagrange_interp_mat(s2, s1)
    Jt = lagrange_interp_mat(t2, t1)

    ⊗(Jt, Js, Jr)
end
#
