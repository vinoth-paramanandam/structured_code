module weno
    use constant
    use declaration
    
    implicit none
    contains

    subroutine weno5(q_stencil, ul, ur)
        !loop iters
        integer(i8) :: k
        !weno input and output variables
        real(prec),dimension(-2:3,4), intent(in) :: q_stencil
        real(prec), dimension(4), intent(out) :: ul, ur
        !weno reconstruction variables
        real(prec), dimension(4) :: g_0, g_1, g_2
        real(prec), dimension(4) :: is0, is1, is2
        real(prec), dimension(4) :: w0, w1, w2
        real(prec), dimension(4) :: alpha1, alpha2, alpha3
        !reconstruction for ul value
        do k = 1,4
            g_0(k) = sixth*(two*q_stencil(-2,k) - seven*q_stencil(-1, k) + eleven*q_stencil(0, k))
            g_1(k) = sixth*(-q_stencil(-1,k) + five*q_stencil(0,k) + two*q_stencil(1,k))
            g_2(k) = sixth*(two*q_stencil(0,k) + five*q_stencil(1,k) - q_stencil(2,k))

            is0(k) = k_is*(q_stencil(-2,k) - two*q_stencil(-1,k) + q_stencil(0,k))**2 +&
            & quarter*(q_stencil(-2,k) - four*q_stencil(-1,k) + three*q_stencil(0,k))**2

            is1(k) = k_is*(q_stencil(-1,k) - two*q_stencil(0,k) + q_stencil(1,k))**2 +&
            & quarter*(q_stencil(-1,k) - q_stencil(1,k))**2

            is2(k) = k_is*(q_stencil(0,k) - two*q_stencil(1,k) + q_stencil(2,k))**2 +&
            & quarter*(three*q_stencil(0,k) - four*q_stencil(1,k) + q_stencil(2,k))**2

            alpha1(k) = c0_js/(is0(k) + e_weno)**2
            alpha2(k) = c1_js/(is1(k) + e_weno)**2
            alpha3(k) = c2_js/(is2(k) + e_weno)**2

            w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
            w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
            w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
        end do
            ul = w0*g_0 + w1*g_1 + w2*g_2
        !reconstruction for ur value
        do k = 1,4
            g_0(k) = sixth*(two*q_stencil(3,k) - seven*q_stencil(2, k) + eleven*q_stencil(1, k))
            g_1(k) = sixth*(-q_stencil(2,k) + five*q_stencil(1,k) + two*q_stencil(0,k))
            g_2(k) = sixth*(two*q_stencil(1,k) + five*q_stencil(0,k) - q_stencil(-1,k))

            is0(k) = k_is*(q_stencil(3,k) - two*q_stencil(2,k) + q_stencil(1,k))**2 +&
            & quarter*(q_stencil(3,k) - four*q_stencil(2,k) + three*q_stencil(1,k))**2

            is1(k) = k_is*(q_stencil(2,k) - two*q_stencil(1,k) + q_stencil(0,k))**2 +&
            & quarter*(q_stencil(2,k) - q_stencil(0,k))**2

            is2(k) = k_is*(q_stencil(1,k) - two*q_stencil(0,k) + q_stencil(-1,k))**2 +&
            & quarter*(three*q_stencil(1,k) - four*q_stencil(0,k) + q_stencil(-1,k))**2

            alpha1(k) = c0_js/(is0(k) + e_weno)**2
            alpha2(k) = c1_js/(is1(k) + e_weno)**2
            alpha3(k) = c2_js/(is2(k) + e_weno)**2

            w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
            w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
            w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
        end do
            ur = w0*g_0 + w1*g_1 + w2*g_2
    end subroutine

    subroutine wenozq(q_stencil, ul, ur)
        implicit none
        !loop iters
        integer(i8) :: k
        !weno input and output variables
        real(prec),dimension(-2:3,4), intent(in) :: q_stencil
        real(prec), dimension(4), intent(out) :: ul, ur
        !weno reconstruction variables
        real(prec), dimension(4) :: g_0, g_1, g_2
        real(prec), dimension(4) :: is0, is1, is2
        real(prec), dimension(4) :: w0, w1, w2
        real(prec), dimension(4) :: alpha1, alpha2, alpha3
        real(prec), dimension(4) :: tau
        real(prec) :: a1, a2, a3, a4

    do k = 1,4
    a1 = -82*(q_stencil(-1, k) - q_stencil(1, k)) + 11*(q_stencil(-2, k) - q_stencil(2, k))
    a2 = 40*(q_stencil(-1, k) + q_stencil(1, k)) - three*(q_stencil(-2, k) + q_stencil(2, k)) - 74*q_stencil(0,k)
    a3 = two*(q_stencil(-1,k) - q_stencil(1,k)) - q_stencil(-2,k) + q_stencil(2, k)
    a4 = -four*(q_stencil(-1,k) + q_stencil(1, k)) + six*q_stencil(0, k) + q_stencil(-2, k) + q_stencil(2, k)

    g_0(k) = q_stencil(0, k) + (a1/120)*half + (a2/56)*(one/six) + (a3/12)*(one/20) + (a4/24)*(one/70)
    g_1(k) = q_stencil(0, k) + (q_stencil(0, k) - q_stencil(-1, k))*half
    g_2(k) = q_stencil(0, k) + (q_stencil(1, k) - q_stencil(0, k))*half

    is0(k) = (a1 + a3)*(a1 + a3)/14400 + (13/3)*(a2/56 + (123/455)*(a4/24))**2 + (781/20)*(a3/12)**2 + (1421461/2275)*(a4/24)**2
    is1(k) = (q_stencil(-1, k) - q_stencil(0, k))**2
    is2(k) = (q_stencil(0, k) - q_stencil(1, k))**2

    tau(k) = quarter*(abs(is0(k) - is1(k)) + abs(is0(k) - is2(k)))**2

    alpha1(k) = c0*(one + tau(k)/(is0(k) + e_weno))
    alpha2(k) = c1*(one + tau(k)/(is1(k) + e_weno))
    alpha3(k) = c2*(one + tau(k)/(is2(k) + e_weno))

    w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    end do
    ul = w0*((one/c0)*g_0 - (c1/c0)*g_1 - (c2/c0)*g_2) + w1*g_1 + w2*g_2

    do k = 1,4
    a1 = -82*(q_stencil(0, k) - q_stencil(2, k)) + 11*(q_stencil(-1, k) - q_stencil(3, k))
    a2 = 40*(q_stencil(0, k) + q_stencil(2, k)) - three*(q_stencil(-1, k) + q_stencil(3, k)) - 74*q_stencil(1,k)
    a3 = two*(q_stencil(0,k) - q_stencil(2,k)) - q_stencil(-1,k) + q_stencil(2, k)
    a4 = -four*(q_stencil(0,k) + q_stencil(2, k)) + six*q_stencil(1, k) + q_stencil(-1, k) + q_stencil(3, k)

    g_0(k) = q_stencil(1, k) - (a1/120)*half + (a2/56)*(one/six) - (a3/12)*(one/20) + (a4/24)*(one/70)
    g_1(k) = q_stencil(1, k) - (q_stencil(1, k) - q_stencil(0, k))*half
    g_2(k) = q_stencil(1, k) - (q_stencil(2, k) - q_stencil(1, k))*half

    is0(k) = (a1 + a3)*(a1 + a3)/14400 + (13/3)*(a2/56 + (123/455)*(a4/24))**2 + (781/20)*(a3/12)**2 + (1421461/2275)*(a4/24)**2
    is1(k) = (q_stencil(0, k) - q_stencil(1, k))**2
    is2(k) = (q_stencil(1, k) - q_stencil(2, k))**2

    tau(k) = quarter*(abs(is0(k) - is1(k)) + abs(is0(k) - is2(k)))**2

    alpha1(k) = c0*(one + tau(k)/(is0(k) + e_weno))
    alpha2(k) = c1*(one + tau(k)/(is1(k) + e_weno))
    alpha3(k) = c2*(one + tau(k)/(is2(k) + e_weno))

    w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    end do
    ur = w0*((one/c0)*g_0 - (c1/c0)*g_1 - (c2/c0)*g_2) + w1*g_1 + w2*g_2
        
    end subroutine wenozq

    subroutine wenozq_gen(q_stencil, ul, ur, xdist1, xdist2)
        implicit none
        !loop iters
        integer(i8) :: k
        !weno input and output variables
        real(prec),dimension(-2:3,4), intent(in) :: q_stencil
        real(prec), INTENT(IN) :: xdist1, xdist2
        real(prec), dimension(4), intent(out) :: ul, ur
        !weno reconstruction variables
        real(prec), dimension(4) :: g_0, g_1, g_2
        real(prec), dimension(4) :: is0, is1, is2
        real(prec), dimension(4) :: w0, w1, w2
        real(prec), dimension(4) :: alpha1, alpha2, alpha3
        real(prec), dimension(4) :: tau
        real(prec) :: a1, a2, a3, a4
        
    do k = 1,4
    a1 = -82.0d0*(q_stencil(-1, k) - q_stencil(1, k)) + 11.0d0*(q_stencil(-2, k) - q_stencil(2, k))
    a2 = 40.0d0*(q_stencil(-1, k) + q_stencil(1, k)) - three*(q_stencil(-2, k) + q_stencil(2, k)) - 74.0d0*q_stencil(0,k)
    a3 = two*(q_stencil(-1,k) - q_stencil(1,k)) - q_stencil(-2,k) + q_stencil(2, k)
    a4 = -four*(q_stencil(-1,k) + q_stencil(1, k)) + six*q_stencil(0, k) + q_stencil(-2, k) + q_stencil(2, k)

    g_0(k) = q_stencil(0, k) + (a1/120.0d0)*xdist1 + (a2/56.0d0)*(xdist1*xdist1 - a21) + &
            (a3/12.0d0)*(xdist1**3 - a31*(xdist1)) + (a4/24.0d0)*(xdist1**4 - a41*xdist1**2 + a42)
    g_1(k) = q_stencil(0, k) + (q_stencil(0, k) - q_stencil(-1, k))*xdist1
    g_2(k) = q_stencil(0, k) + (q_stencil(1, k) - q_stencil(0, k))*xdist1

    is0(k) = (a1 + a3)*(a1 + a3)/14400.0d0 + (13.0d0/3.0d0)*(a2/56.0d0 + (123.0d0/455.0d0)*(a4/24.0d0))**2 + &
             (781.0d0/20.0d0)*(a3/12.0d0)**2 + (1421461.0d0/2275.0d0)*(a4/24.0d0)**2
    is1(k) = (q_stencil(-1, k) - q_stencil(0, k))**2
    is2(k) = (q_stencil(0, k) - q_stencil(1, k))**2

    tau(k) = quarter*(abs(is0(k) - is1(k)) + abs(is0(k) - is2(k)))**2

    alpha1(k) = c0*(one + tau(k)/(is0(k) + e_weno))
    alpha2(k) = c1*(one + tau(k)/(is1(k) + e_weno))
    alpha3(k) = c2*(one + tau(k)/(is2(k) + e_weno))

    w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    end do
    ul = w0*((one/c0)*g_0 - (c1/c0)*g_1 - (c2/c0)*g_2) + w1*g_1 + w2*g_2

    do k = 1,4
    a1 = -82.0d0*(q_stencil(0, k) - q_stencil(2, k)) + 11.0d0*(q_stencil(-1, k) - q_stencil(3, k))
    a2 = 40.0d0*(q_stencil(0, k) + q_stencil(2, k)) - three*(q_stencil(-1, k) + q_stencil(3, k)) - 74.0d0*q_stencil(1,k)
    a3 = two*(q_stencil(0,k) - q_stencil(2,k)) - q_stencil(-1,k) + q_stencil(2, k)
    a4 = -four*(q_stencil(0,k) + q_stencil(2, k)) + six*q_stencil(1, k) + q_stencil(-1, k) + q_stencil(3, k)

    g_0(k) = q_stencil(1, k) + (a1/120.0d0)*xdist2 + (a2/56.0d0)*(xdist2*xdist2 - a21) + &
    (a3/12.0d0)*(xdist2**3 - a31*(xdist2)) + (a4/24.0d0)*(xdist2**4 - a41*xdist2**2 + a42)
    g_1(k) = q_stencil(1, k) + (q_stencil(1, k) - q_stencil(0, k))*xdist2
    g_2(k) = q_stencil(1, k) + (q_stencil(2, k) - q_stencil(1, k))*xdist2

    is0(k) = (a1 + a3)*(a1 + a3)/14400.0d0 + (13.0d0/3.0d0)*(a2/56.0d0 + (123.0d0/455.0d0)*(a4/24.0d0))**2 + &
    (781.0d0/20.0d0)*(a3/12.0d0)**2 + (1421461.0d0/2275.0d0)*(a4/24.0d0)**2
    is1(k) = (q_stencil(0, k) - q_stencil(1, k))**2
    is2(k) = (q_stencil(1, k) - q_stencil(2, k))**2

    tau(k) = quarter*(abs(is0(k) - is1(k)) + abs(is0(k) - is2(k)))**2

    alpha1(k) = c0*(one + tau(k)/(is0(k) + e_weno))
    alpha2(k) = c1*(one + tau(k)/(is1(k) + e_weno))
    alpha3(k) = c2*(one + tau(k)/(is2(k) + e_weno))

    w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    end do
    ur = w0*((one/c0)*g_0 - (c1/c0)*g_1 - (c2/c0)*g_2) + w1*g_1 + w2*g_2
    end subroutine wenozq_gen

    subroutine wenozq_gaussl(q_stencil, ul, ur, xdist1, xdist2)
        implicit none
        !loop iters
        integer(i8) :: k
        !weno input and output variables
        real(prec),dimension(-2:2,4), intent(in) :: q_stencil
        real(prec), INTENT(IN) :: xdist1, xdist2
        real(prec), dimension(4), intent(out) :: ul, ur
        !weno reconstruction variables
        real(prec), dimension(4) :: g_0, g_1, g_2
        real(prec), dimension(4) :: is0, is1, is2
        real(prec), dimension(4) :: w0, w1, w2
        real(prec), dimension(4) :: alpha1, alpha2, alpha3
        real(prec), dimension(4) :: tau
        real(prec) :: a1, a2, a3, a4
        
    do k = 1,4
    a1 = -82.0d0*(q_stencil(-1, k) - q_stencil(1, k)) + 11.0d0*(q_stencil(-2, k) - q_stencil(2, k))
    a2 = 40.0d0*(q_stencil(-1, k) + q_stencil(1, k)) - three*(q_stencil(-2, k) + q_stencil(2, k)) - 74.0d0*q_stencil(0,k)
    a3 = two*(q_stencil(-1,k) - q_stencil(1,k)) - q_stencil(-2,k) + q_stencil(2, k)
    a4 = -four*(q_stencil(-1,k) + q_stencil(1, k)) + six*q_stencil(0, k) + q_stencil(-2, k) + q_stencil(2, k)

    g_0(k) = q_stencil(0, k) + (a1/120.0d0)*xdist1 + (a2/56.0d0)*(xdist1*xdist1 - a21) + &
            (a3/12.0d0)*(xdist1**3 - a31*(xdist1)) + (a4/24.0d0)*(xdist1**4 - a41*xdist1**2 + a42)
    g_1(k) = q_stencil(0, k) + (q_stencil(0, k) - q_stencil(-1, k))*xdist1
    g_2(k) = q_stencil(0, k) + (q_stencil(1, k) - q_stencil(0, k))*xdist1

    is0(k) = (a1 + a3)*(a1 + a3)/14400.0d0 + (13.0d0/3.0d0)*(a2/56.0d0 + (123.0d0/455.0d0)*(a4/24.0d0))**2 + &
             (781.0d0/20.0d0)*(a3/12.0d0)**2 + (1421461.0d0/2275.0d0)*(a4/24.0d0)**2
    is1(k) = (q_stencil(-1, k) - q_stencil(0, k))**2
    is2(k) = (q_stencil(0, k) - q_stencil(1, k))**2

    tau(k) = quarter*(abs(is0(k) - is1(k)) + abs(is0(k) - is2(k)))**2

    alpha1(k) = c0*(one + tau(k)/(is0(k) + e_weno))
    alpha2(k) = c1*(one + tau(k)/(is1(k) + e_weno))
    alpha3(k) = c2*(one + tau(k)/(is2(k) + e_weno))

    w0(k) = alpha1(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w1(k) = alpha2(k)/(alpha1(k)+alpha2(k)+alpha3(k))
    w2(k) = alpha3(k)/(alpha1(k)+alpha2(k)+alpha3(k))

    ul(k) = w0(k)*((one/c0)*g_0(k) - (c1/c0)*g_1(k) - (c2/c0)*g_2(k)) + w1(k)*g_1(k) + w2(k)*g_2(k)
    
    g_0(k) = q_stencil(0, k) + (a1/120.0d0)*xdist2 + (a2/56.0d0)*(xdist2*xdist2 - a21) + &
            (a3/12.0d0)*(xdist2**3 - a31*(xdist2)) + (a4/24.0d0)*(xdist2**4 - a41*xdist2**2 + a42)
    g_1(k) = q_stencil(0, k) + (q_stencil(0, k) - q_stencil(-1, k))*xdist2
    g_2(k) = q_stencil(0, k) + (q_stencil(1, k) - q_stencil(0, k))*xdist2

    ur(k) = w0(k)*((one/c0)*g_0(k) - (c1/c0)*g_1(k) - (c2/c0)*g_2(k)) + w1(k)*g_1(k) + w2(k)*g_2(k)
    end do
end subroutine wenozq_gaussl

end module
