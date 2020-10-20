
module boundary
    use constant
    use declaration
    implicit none

    contains
    subroutine zero_extrapolation(istart, iend, jstart, jend, bside)
        implicit none
        integer(i32), intent(in) :: istart, iend, jstart, jend, bside
        integer(i32) :: i, j
        select case (bside)
        case(1)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(i, gc+1)
                    u(i, j) = u(i, gc+1)
                    v(i, j) = v(i, gc+1)
                    p(i, j) = p(i, gc+1)
                    t(i, j) = t(i, gc+1)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(2)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(nx-gc, j)
                    u(i, j) = abs(u(nx-gc, j))
                    v(i, j) = v(nx-gc, j)
                    p(i, j) = p(nx-gc, j)
                    t(i, j) = t(nx-gc, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(3)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(i, ny-gc)
                    u(i, j) = u(i, ny-gc)
                    v(i, j) = v(i, ny-gc)
                    p(i, j) = p(i, ny-gc)
                    t(i, j) = t(i, ny-gc)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(4)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(gc+1, j)
                    u(i, j) = u(gc+1, j)
                    v(i, j) = v(gc+1, j)
                    p(i, j) = p(gc+1, j)
                    t(i, j) = t(gc+1, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        end select
    end subroutine zero_extrapolation

    subroutine symmetry(istart, iend, jstart, jend, bside)
        implicit none
        integer(i32), intent(in) :: istart, iend, jstart, jend, bside
        real(prec) :: dxlen, dylen, dlen, normalx, normaly, un, ut
        integer(i32) :: i, j
        select case (bside)
        case(1)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i+1, gc+1) + x(i, gc+1))
                    dylen = (y(i+1, gc+1) - y(i, gc+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(i, gc+1)
                    u(i, j) = u(i, gc+1)
                    v(i, j) = v(i, gc+1)
                    p(i, j) = p(i, gc+1)
                    t(i, j) = t(i, gc+1)

                    un = u(i, j)*normalx + v(i, j)*normaly
                    ut = v(i, j)*normalx - u(i, j)*normaly

                    un = -un 
                    ut = ut

                    u(i, j) = un*normalx - ut*normaly
                    v(i, j) = ut*normalx + un*normaly

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(2)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(nx-gc, j) + x(nx-gc, j+1))
                    dylen = (y(nx-gc, j) - y(nx-gc, j+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(nx-gc, j)
                    u(i, j) = u(nx-gc, j)
                    v(i, j) = v(nx-gc, j)
                    p(i, j) = p(nx-gc, j)
                    t(i, j) = t(nx-gc, j)

                    un = u(i, j)*normalx + v(i, j)*normaly
                    ut = v(i, j)*normalx - u(i, j)*normaly

                    un = -un 
                    ut = ut

                    u(i, j) = un*normalx - ut*normaly
                    v(i, j) = ut*normalx + un*normaly

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(3)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i, ny-gc) + x(i+1, ny-gc))
                    dylen = (y(i, ny-gc) - y(i+1, ny-gc))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(i, ny-gc)
                    u(i, j) = u(i, ny-gc)
                    v(i, j) = v(i, ny-gc)
                    p(i, j) = p(i, ny-gc)
                    t(i, j) = t(i, ny-gc)

                    un = u(i, j)*normalx + v(i, j)*normaly
                    ut = v(i, j)*normalx - u(i, j)*normaly

                    un = -un 
                    ut = ut

                    u(i, j) = un*normalx - ut*normaly
                    v(i, j) = ut*normalx + un*normaly

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(4)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = -(x(gc+1, j) - x(gc+1, j+1))
                    dylen = y(gc+1, j) - y(gc+1, j+1)
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(gc+1, j)
                    u(i, j) = u(gc+1, j)
                    v(i, j) = v(gc+1, j)
                    p(i, j) = p(gc+1, j)
                    t(i, j) = t(gc+1, j)

                    un = u(i, j)*normalx + v(i, j)*normaly
                    ut = v(i, j)*normalx - u(i, j)*normaly

                    un = -un 
                    ut = ut

                    u(i, j) = un*normalx - ut*normaly
                    v(i, j) = ut*normalx + un*normaly

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        end select
    end subroutine symmetry

    subroutine wallboundary(istart, iend, jstart, jend, bside)
        implicit none
        integer(i32), intent(in) :: istart, iend, jstart, jend, bside
        real(prec) :: dxlen, dylen, dlen, normalx, normaly, un, ut
        integer(i32) :: i, j
        select case (bside)
        case(1)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i+1, gc+1) + x(i, gc+1))
                    dylen = (y(i+1, gc+1) - y(i, gc+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(i, gc+1)
                    u(i, j) = -u(i, gc+1)
                    v(i, j) = -v(i, gc+1)
                    p(i, j) = p(i, gc+1)
                    t(i, j) = two*288.16_prec - t(i, gc+1)
                    ! t(i, j) = t(i, gc+1)

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(2)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(nx-gc, j) + x(nx-gc, j+1))
                    dylen = (y(nx-gc, j) - y(nx-gc, j+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(nx-gc, j)
                    u(i, j) = -u(nx-gc, j)
                    v(i, j) = -v(nx-gc, j)
                    p(i, j) = p(nx-gc, j)
                    t(i, j) = t(nx-gc, j)

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(3)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i, ny-gc) + x(i+1, ny-gc))
                    dylen = (y(i, ny-gc) - y(i+1, ny-gc))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(i, ny-gc)
                    u(i, j) = -u(i, ny-gc)
                    v(i, j) = -v(i, ny-gc)
                    p(i, j) = p(i, ny-gc)
                    t(i, j) = t(i, ny-gc)

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(4)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = -(x(gc+1, j) - x(gc+1, j+1))
                    dylen = y(gc+1, j) - y(gc+1, j+1)
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(gc+1, j)
                    u(i, j) = -u(gc+1, j)
                    v(i, j) = -v(gc+1, j)
                    p(i, j) = p(gc+1, j)
                    t(i, j) = t(gc+1, j)

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        end select
    end subroutine wallboundary

    subroutine supersonic_inflow(istart, iend, jstart, jend, bside, rin, uin, vin, pin)
        implicit none
        integer(i32), intent(in) :: istart, iend, jstart, jend, bside
        real(prec), intent(in) :: rin, uin, vin, pin
        integer(i32) :: i, j

        do j = jstart, jend
            do i = istart, iend
                    rho(i, j) = rin
                    u(i, j) = uin
                    v(i, j) = vin
                    p(i, j) = pin
                    t(i, j) = pin/(r*rin)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            enddo
        enddo
    end subroutine supersonic_inflow

    subroutine pressure_oulet(istart, iend, jstart, jend, bside)
        implicit none
        integer(i32) :: istart, iend, jstart, jend, bside
        real(prec) :: dxlen, dylen, dlen, normalx, normaly
        integer(i32) :: i, j
        real(prec) :: pdom, rdom, udom, vdom, cref, rref, m, qn
        real(prec) :: pamb, uamb, vamb, ramb
        select case (bside)
        case(1)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i+1, gc+1) + x(i, gc+1))
                    dylen = (y(i+1, gc+1) - y(i, gc+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(i, gc+1)
                    u(i, j) = u(i, gc+1)
                    v(i, j) = v(i, gc+1)
                    p(i, j) = p(i, gc+1)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                    pdom = p(i, j)
                    rdom = rho(i, j)
                    vdom = v(i, j)
                    udom = u(i, j)

                    rref = rho(i, j)
                    cref = c(i,j)

                    m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                    qn = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen

                    if (m < one) then
                        if(qn < zero) then
                            p(i,j) = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                            rho(i,j) = ramb + (p(i,j) - pamb)/cref**2
                            u(i,j) = uamb - normalx*(pamb - p(i,j))/(rref*cref)
                            v(i,j) = vamb - normaly*(pamb - p(i,j))/(rref*cref)
                        else
                            p(i,j) = pamb
                            rho(i,j) = rdom + (p(i,j) - pdom)/cref**2
                            u(i,j) = udom + normalx*(pdom - pamb)/(rref*cref)
                            v(i,j) = vdom + normaly*(pdom - pamb)/(rref*cref)
                        endif
                    endif
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(2)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(nx-gc, j) + x(nx-gc, j+1))
                    dylen = (y(nx-gc, j) - y(nx-gc, j+1))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    pamb = 101325.0_prec
                    ramb = 1.225_prec
                    uamb = zero
                    vamb = zero

                    rho(i, j) = rho(nx-gc, j)
                    u(i, j) = u(nx-gc, j)
                    v(i, j) = v(nx-gc, j)
                    p(i, j) = p(nx-gc, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                    pdom = p(i, j)
                    rdom = rho(i, j)
                    vdom = v(i, j)
                    udom = u(i, j)

                    rref = rho(i, j)
                    cref = c(i,j)

                    m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                    qn = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen

                    if (m < one) then
                        if(qn < zero) then
                            p(i,j) = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                            rho(i,j) = ramb + (p(i,j) - pamb)/cref**2
                            u(i,j) = uamb - normalx*(pamb - p(i,j))/(rref*cref)
                            v(i,j) = vamb - normaly*(pamb - p(i,j))/(rref*cref)
                        else
                            p(i,j) = pamb
                            rho(i,j) = rdom + (p(i,j) - pdom)/cref**2
                            u(i,j) = udom + normalx*(pdom - pamb)/(rref*cref)
                            v(i,j) = vdom + normaly*(pdom - pamb)/(rref*cref)
                        endif
                    endif
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(3)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = (-x(i, ny-gc) + x(i+1, ny-gc))
                    dylen = (y(i, ny-gc) - y(i+1, ny-gc))
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    pamb = 101325.0_prec
                    ramb = 1.225_prec
                    uamb = zero
                    vamb = zero

                    rho(i, j) = rho(i, ny-gc)
                    u(i, j) = u(i, ny-gc)
                    v(i, j) = v(i, ny-gc)
                    p(i, j) = p(i, ny-gc)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                    pdom = p(i, j)
                    rdom = rho(i, j)
                    vdom = v(i, j)
                    udom = u(i, j)

                    rref = rho(i, j)
                    cref = c(i,j)

                    m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                    qn = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen

                    if (m < one) then
                        if(qn < zero) then
                            p(i,j) = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                            rho(i,j) = ramb + (p(i,j) - pamb)/cref**2
                            u(i,j) = uamb - normalx*(pamb - p(i,j))/(rref*cref)
                            v(i,j) = vamb - normaly*(pamb - p(i,j))/(rref*cref)
                        else
                            p(i,j) = pamb
                            rho(i,j) = rdom + (p(i,j) - pdom)/cref**2
                            u(i,j) = udom + normalx*(pdom - pamb)/(rref*cref)
                            v(i,j) = vdom + normaly*(pdom - pamb)/(rref*cref)
                        endif
                    endif
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(4)
            do j = jstart, jend
                do i = istart, iend
                    dxlen = -(x(gc+1, j) - x(gc+1, j+1))
                    dylen = y(gc+1, j) - y(gc+1, j+1)
                    dlen = dsqrt(dxlen**2 + dylen**2)
                    normalx = dylen/dlen
                    normaly = dxlen/dlen

                    rho(i, j) = rho(gc+1, j)
                    u(i, j) = u(gc+1, j)
                    v(i, j) = v(gc+1, j)
                    p(i, j) = p(gc+1, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                    pdom = p(i, j)
                    rdom = rho(i, j)
                    vdom = v(i, j)
                    udom = u(i, j)

                    rref = rho(i, j)
                    cref = c(i,j)

                    m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                    qn = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen

                    if (m < one) then
                        if(qn < zero) then
                            p(i,j) = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                            rho(i,j) = ramb + (p(i,j) - pamb)/cref**2
                            u(i,j) = uamb - normalx*(pamb - p(i,j))/(rref*cref)
                            v(i,j) = vamb - normaly*(pamb - p(i,j))/(rref*cref)
                        else
                            p(i,j) = pamb
                            rho(i,j) = rdom + (p(i,j) - pdom)/cref**2
                            u(i,j) = udom + normalx*(pdom - pamb)/(rref*cref)
                            v(i,j) = vdom + normaly*(pdom - pamb)/(rref*cref)
                        endif
                    endif
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        end select
    end subroutine pressure_oulet

    subroutine periodicboundary(istart, iend, jstart, jend, bside)
        implicit none
        integer(i32), intent(in) :: istart, iend, jstart, jend, bside
        integer(i32) :: i, j

        select case (bside)
        case(1)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(i, ny-2*gc+j)
                    u(i, j) = u(i, ny-2*gc+j)
                    v(i, j) = v(i, ny-2*gc+j)
                    p(i, j) = p(i, ny-2*gc+j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(2)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(i-nx+2*gc, j)
                    u(i, j) = u(i-nx+2*gc, j)
                    v(i, j) = v(i-nx+2*gc, j)
                    p(i, j) = p(i-nx+2*gc, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(3)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(i, j-ny+2*gc)
                    u(i, j) = u(i, j-ny+2*gc)
                    v(i, j) = v(i, j-ny+2*gc)
                    p(i, j) = p(i, j-ny+2*gc)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        case(4)
            do j = jstart, jend
                do i = istart, iend
                    rho(i, j) = rho(nx-2*gc+i, j)
                    u(i, j) = u(nx-2*gc+i, j)
                    v(i, j) = v(nx-2*gc+i, j)
                    p(i, j) = p(nx-2*gc+i, j)
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                    t(i, j) = p(i,j)/(r*rho(i,j))

                    q(1, i, j) = rho(i, j)
                    q(2, i, j) = rho(i, j)*u(i, j)
                    q(3, i, j) = rho(i, j)*v(i, j)
                    q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
                enddo
            enddo
        end select
    end subroutine periodicboundary

    subroutine vortexbc
        !left boundary condtions
        call periodicboundary(1, gc, gc+1, ny-gc, 4)
        !right boundary conditions
        call periodicboundary(nx-gc+1, nx, 1, ny, 2)
        !bottom boundary condition
        call periodicboundary(gc+1, nx-gc, 1, gc, 1)
        !top boundary condtions
        call periodicboundary(gc+1, nx-gc, ny-gc+1, ny, 3)
    end subroutine vortexbc

    subroutine nozzlebc
        real(prec) :: temp, tinf, dpres, pin, uin, vin, rin, cref, rref, qnorm, m, cdom
        real(prec) :: rb, ub, vb, pb, rdom, udom, vdom, pdom, pamb, uamb, vamb, ramb, un, ut, normalx, normaly, dlen, dxlen, dylen
        integer :: i, j
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if j = 1, jc = 6; j = 2, jc = 5; j = 3, jc = 4; 
        !!!!!!!!!!!!!!symmetric boundary!!!!!!!!!!!!!!!!!!           jc = 6;        jc = 5;        jc = 4;
        do j = 1, gc
            do i = gc+1, nx-gc
                dxlen = (-x(i+1, gc+1) + x(i, gc+1))
                dylen = (y(i+1, gc+1) - y(i, gc+1))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                rho(i, j) = rho(i, gc+1)
                u(i, j) = u(i, gc+1)
                v(i, j) = v(i, gc+1)
                p(i, j) = p(i, gc+1)
                t(i, j) = t(i, gc+1)

                un = u(i, j)*normalx + v(i, j)*normaly
                ut = v(i, j)*normalx - u(i, j)*normaly

                un = -un 
                ut = ut

                u(i, j) = un*normalx - ut*normaly
                v(i, j) = ut*normalx + un*normaly

                c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)

                q(1, i, j) = rho(i, j)
                q(2, i, j) = rho(i, j)*u(i, j)
                q(3, i, j) = rho(i, j)*v(i, j)
                q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if j = 1, jc = 6; j = 2, jc = 5; j = 3, jc = 4; 
        !!!!!!!!!!!!!!linear ex boundary!!!!!!!!!!!!!!!!!!           jc = 6;        jc = 5;        jc = 4;
        do j = gc+1, ny-gc
            do i = nx-gc+1, nx
                ! rho(i, j) = rho(nx-gc, j)!two*rho(i-1, j) - rho(i-2, j)!rho(nx-gc, j)
                ! u(i, j) = dabs(u(nx-gc, j))!two*dabs(u(i-1, j)) - dabs(u(i-2, j))!abs(u(nx-gc, j))
                ! v(i, j) = v(nx-gc, j)!two*v(i-1, j) - v(i-2, j)!v(nx-gc, j)
                ! p(i, j) = p(nx-gc, j)!two*p(i-1, j) - p(i-2, j)!p(nx-gc, j)
                
                ! ! rho(i, j) = two*rho(i-1, j) - rho(i-2, j)!rho(nx-gc, j)
                ! ! u(i, j) = two*dabs(u(i-1, j)) - dabs(u(i-2, j))!abs(u(nx-gc, j))
                ! ! v(i, j) = two*v(i-1, j) - v(i-2, j)!v(nx-gc, j)
                ! ! p(i, j) = two*p(i-1, j) - p(i-2, j)!p(nx-gc, j)
                
                ! t(i, j) = p(i,j)/(r*rho(i,j))!two*t(i-1, j) - t(i-2, j)!t(nx-gc, j)
                ! c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                ! h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                dxlen = (-x(nx-gc, j) + x(nx-gc, j+1))
                dylen = (y(nx-gc, j) - y(nx-gc, j+1))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                pamb = 101325.0_prec
                ramb = 1.225_prec
                uamb = zero
                vamb = zero

                rho(i, j) = rho(nx-gc, j)
                u(i, j) = dabs(u(nx-gc, j))
                v(i, j) = v(nx-gc, j)
                p(i, j) = p(nx-gc, j)
                c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                ! pdom = p(i, j)
                ! rdom = rho(i, j)
                ! vdom = v(i, j)
                ! udom = u(i, j)

                ! rref = rho(i, j)
                ! cref = c(i,j)

                ! m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                ! qnorm = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen
                ! if (m < one) then
                !     if(qnorm < zero) then
                !         pb = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                !         rb = ramb + (pb - pamb)/cref**2
                !         ub = uamb - normalx*(pamb - pb)/(rref*cref)
                !         vb = vamb - normaly*(pamb - pb)/(rref*cref)
                !     else
                !         pb = pamb
                !         rb = rdom + (pb - pdom)/cref**2
                !         ub = udom + normalx*(pdom - pamb)/(rref*cref)
                !         vb = vdom + normaly*(pdom - pamb)/(rref*cref)
                !     endif
                !     p(i, j) = pb!two*pb - pdom
                !     u(i, j) = ub!two*ub - udom
                !     v(i, j) = vb!two*vb - vdom
                !     rho(i, j) = rb!two*rb - rdom
                ! else
                !     p(i, j) =  pdom
                !     u(i, j) =  DABS(udom)
                !     v(i, j) =  vdom
                !     rho(i, j) =  rdom
                ! endif
                !     ! p(i, j) =  pdom
                !     ! u(i, j) =  DABS(udom)
                !     ! v(i, j) =  vdom
                !     ! rho(i, j) =  rdom
                h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                t(i, j) = p(i,j)/(r*rho(i,j))
                c(i, j) = dsqrt(g1*R*t(i, j))


                q(1, i, j) = rho(i, j)
                q(2, i, j) = rho(i, j)*u(i, j)
                q(3, i, j) = rho(i, j)*v(i, j)
                q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if j = 1, jc = 6; j = 2, jc = 5; j = 3, jc = 4; 
        !!!!!!!!!!!!!!far field boundary!!!!!!!!!!!!!!!!!!           jc = 6;        jc = 5;        jc = 4;
        do j = ny-gc+1, ny
            do i = gc+1, nx-gc
                dxlen = (-x(i, ny-gc) + x(i+1, ny-gc))
                dylen = (y(i, ny-gc) - y(i+1, ny-gc))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                pamb = 101325.0_prec
                ramb = 1.225_prec
                uamb = zero
                vamb = zero

                rho(i, j) = rho(i, ny-gc)
                u(i, j) = u(i, ny-gc)
                v(i, j) = v(i, ny-gc)
                p(i, j) = p(i, ny-gc)
                c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))

                pdom = p(i, j)
                rdom = rho(i, j)
                vdom = v(i, j)
                udom = u(i, j)

                rref = rho(i, j)
                cref = c(i,j)

                m = dsqrt(u(i,j)**2 + v(i,j)**2)/c(i,j)
                qnorm = u(i,j)*normalx*dlen + v(i,j)*normaly*dlen

                if (m < one) then
                    if(qnorm < zero) then
                        pb = half*(pdom + pamb - rref*cref*(normalx*(uamb -udom) + normaly*(vamb - vdom)))
                        rb = ramb + (pb - pamb)/cref**2
                        ub = uamb - normalx*(pamb - pb)/(rref*cref)
                        vb = vamb - normaly*(pamb - pb)/(rref*cref)
                    else
                        pb = pamb
                        rb = rdom + (pb - pdom)/cref**2
                        ub = udom + normalx*(pdom - pamb)/(rref*cref)
                        vb = vdom + normaly*(pdom - pamb)/(rref*cref)
                    endif
                    p(i, j) = pb!two*pb - pdom
                    u(i, j) = ub!two*ub - udom
                    v(i, j) = vb!two*vb - vdom
                    rho(i, j) = rb!two*rb - rdom
                endif

                p(i, j) = pamb!two*pb - pdom
                u(i, j) = uamb!two*ub - udom
                v(i, j) = vamb!two*vb - vdom
                rho(i, j) = ramb!two*rb - rdom
                h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                t(i, j) = p(i,j)/(r*rho(i,j))
                c(i, j) = dsqrt(g1*R*t(i, j))

                q(1, i, j) = rho(i, j)
                q(2, i, j) = rho(i, j)*u(i, j)
                q(3, i, j) = rho(i, j)*v(i, j)
                q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            enddo
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if j = 1, jc = 6; j = 2, jc = 5; j = 3, jc = 4; 
        !!!!!!!!!!!supersonic inlet and wall!!!!!!!!!!!!!!           jc = 6;        jc = 5;        jc = 4;
        temp = (1 + (g2)*half*mach_no**2)
        tinf = 300.0d0/temp
        dpres = temp**3.50d0
        pin = npr*101325.0d0/dpres
        rin = pin/(r*tinf)
        uin = mach_no*(dsqrt(g1*pin/rin))
        vin = zero
        ! print *, rin, uin, vin, tinf
        do j = gc+1, ny-gc
            do i = 1, gc
                dxlen = -(x(gc+1, j) - x(gc+1, j+1))
                dylen = y(gc+1, j) - y(gc+1, j+1)
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen
                if(yc(i, j) .le. 0.25d0) then
                    rho(i, j) = rin
                    ! if(counter .le. 500) then
                    !     if(yc(i, j) .le. 0.0625d0) then
                    !         rho(i, j) = rin - 0.35d0!0.25d0 + rin !0.25d0 + rin
                    !     endif
                    ! endif
                    u(i, j) = uin
                    v(i, j) = vin
                    p(i, j) = pin
                    t(i, j) = tinf
                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                else
                    rho(i, j) = rho(gc+1, j)
                    u(i, j) = u(gc+1, j)
                    v(i, j) = v(gc+1, j)
                    p(i, j) = p(gc+1, j)
                    t(i, j) = t(gc+1, j)

                    un = u(i, j)*normalx + v(i, j)*normaly
                    ut = v(i, j)*normalx - u(i, j)*normaly

                    un = -un 
                    ut = ut

                    u(i, j) = un*normalx - ut*normaly
                    v(i, j) = ut*normalx + un*normaly

                    c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                    h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                endif
                q(1, i, j) = rho(i, j)
                q(2, i, j) = rho(i, j)*u(i, j)
                q(3, i, j) = rho(i, j)*v(i, j)
                q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            enddo
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! !bottom boundary condition
        ! call symmetry(gc+1, nx-gc, 1, gc, 1)
        ! !right boundary condition
        ! call zero_extrapolation(nx-gc+1, nx, gc+1, ny-gc, 2)
        ! !top boundary condition
        ! call pressure_oulet(gc+1, nx-gc, ny-gc+1, ny, 3)
        ! !left boundary condition
        ! temp = (1 + (g2)*half*mach_no**2)
        ! tinf = 300.0d0/temp
        ! dpres = temp**3.50d0
        ! pin = npr*101325.0d0/dpres
        ! rin = pin/(r*tinf)
        ! uin = mach_no*(dsqrt(g1*pin/rin))
        ! vin = zero
        ! call supersonic_inflow(1, gc, gc+1, (ny-gc)/2, 4, rin, uin, vin, pin)
        ! call wallboundary(1, gc, (ny-gc)/2+1, ny-gc, 4)
        
        ! integer :: i, j
        ! !bottom boundary condition
        ! call wallboundary(gc+1, nx-gc, 1, gc, 1)
        ! !right boundary condition
        ! call zero_extrapolation(nx-gc+1, nx, gc+1, ny-gc, 2)
        ! !top boundary condition
        ! do j = ny-gc+1, ny
        !     do i = gc+1, nx-gc
        !         t(i,j) = 288.16_prec
        !         p(i,j) = 101325_prec
        !         rho(i,j) = p(i,j)/(r*t(i,j))
        !         v(i,j) = zero
        !         c(i, j) = sqrt(g1*p(i, j)/rho(i, j))
        !         u(i,j) = mach_no*c(i,j)
        !         h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
        !         q(1, i, j) = rho(i, j)
        !         q(2, i, j) = rho(i, j)*u(i, j)
        !         q(3, i, j) = rho(i, j)*v(i, j)
        !         q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
        !     enddo
        ! enddo
        ! ! call zero_extrapolation(gc+1, nx-gc, ny-gc+1, ny, 3)
        ! !left boundary condition
        ! do j = gc+1, ny-gc
        !     do i = 1, gc
        !         t(i,j) = 288.16_prec
        !         p(i,j) = 101325_prec
        !         rho(i,j) = p(i,j)/(r*t(i,j))
        !         v(i,j) = zero
        !         c(i, j) = sqrt(g1*p(i, j)/rho(i, j))
        !         u(i,j) = mach_no*c(i,j)
        !         h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
        !         q(1, i, j) = rho(i, j)
        !         q(2, i, j) = rho(i, j)*u(i, j)
        !         q(3, i, j) = rho(i, j)*v(i, j)
        !         q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
        !     enddo
        ! enddo
        ! ! call zero_extrapolation(1, gc, gc+1, ny-gc, 4)
        ! ! temp = (1 + (g2)*half*mach_no**2)
        ! ! tinf = 300.0d0/temp
        ! ! dpres = temp**3.5
        ! ! pin = npr*101325.0_prec/dpres
        ! ! rin = pin/(r*tinf)
        ! ! uin = mach_no*(dsqrt(g1*pin/rin))
        ! ! vin = zero
        ! ! call supersonic_inflow(1, gc, gc+1, (ny-gc)/4, 4, rin, uin, vin, pin)
        ! ! call wallboundary(1, gc, (ny-gc)/4+1, ny-gc, 4)
    end subroutine nozzlebc

    subroutine flatplatebc
        real(prec) :: temp, tinf, dpres, pin, uin, vin, rin
        !bottom boundary condition
        call symmetry(gc+1, 43, 1, gc, 1)
        call wallboundary(44, nx-gc, 1, gc, 1)
        !right boundary condition
        call zero_extrapolation(nx-gc+1, nx, gc+1, ny-gc, 2)
        !top boundary condition
        call zero_extrapolation(gc+1, nx-gc, ny-gc+1, ny, 3)
        !left boundary condition
        call zero_extrapolation(1, gc, gc+1, ny-gc, 4)
        ! temp = (1 + (g2)*half*mach_no**2)
        ! tinf = 300.0d0/temp
        ! dpres = temp**3.5
        ! pin = npr*101325.0_prec/dpres
        ! rin = pin/(r*tinf)
        ! uin = mach_no*(dsqrt(g1*pin/rin))
        ! vin = zero
        ! call supersonic_inflow(1, gc, gc+1, (ny-gc)/4, 4, rin, uin, vin, pin)
        ! call wallboundary(1, gc, (ny-gc)/4+1, ny-gc, 4)
    end subroutine flatplatebc

end module
