module flux
    use constant
    use declaration
    use misc
    use weno

    IMPLICIT NONE

    CONTAINS

    subroutine weno_reconstruction(dir)
    implicit none
    integer(i8), intent(in) :: dir
    !loop iters and dir
    integer(i32) :: i, j, k, m
    !normals for the grid
    real(prec) :: normalx, normaly, dlen, dxlen, dylen
    !roe average or arithmatic avg values
    real(prec) :: u_avg, v_avg, c_avg, h_avg
    !conservative and characterstic variables as vectors
    real(prec), dimension(4) :: ql, qr, qcons, vecl, vecr, vecm
    !stencil for the weno scheme
    real(prec), dimension(-2:3, 4) :: q_stencil
    !eigenvectors of left and right
    real(prec), dimension(4,4) :: leigenvec, reigenvec
    ! xdistance
    real(prec) :: ximh, xiph, xip3h, xi1, xi2
    real(prec) :: yimh, yiph, yip3h, yi1, yi2

    if(dir .eq. 1) then
        !loop for the x-dir flux calculation
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, xc, nx, ny, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc+1, ny-gc
            do i = gc, nx-gc
                !positive normals in the direction
                dxlen = -(x(i+1, j+1) - x(i+1, j))
                dylen = (y(i+1, j+1) - y(i+1, j))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                xiph = (x(i+1, j+1) + x(i+1, j))*half
                ximh = (x(i, j+1) + x(i, j))*half
                xip3h = (x(i+2, j+1) + x(i+2, j))*half

                xi1 = (xiph - xc(i, j))/(xiph - ximh)
                xi2 = (xiph - xc(i+1, j))/(xip3h - xiph)

                !  !calculation of roe or arithmetic averages
                u_avg = roeavg(rho(i,j), rho(i+1, j), u(i, j), u(i+1, j))
                v_avg = roeavg(rho(i,j), rho(i+1, j), v(i, j), v(i+1, j))
                h_avg = roeavg(rho(i,j), rho(i+1, j), h(i, j), h(i+1, j))
                c_avg = sqrt(g2*(h_avg - half*(u_avg**2 + v_avg**2)))

                call lefteigenvector(u_avg, v_avg, c_avg, normalx, normaly, leigenvec)

                !finding the relevant stencil and transforming it to characteristic variables
                do k = -2,3
                    do m = 1, 4
                        qcons(m) = q(m, i+k, j)
                    end do
                    vecm = matmul(leigenvec, qcons)
                    q_stencil(k, :) = vecm!matmul(qcons, leigenvec)
                end do

                !weno reconstruction procedure. the procedure can be found in toro 2005
                !other orders of reconstruction is also included in the module weno
                call weno5(q_stencil, ql, qr)
                ! call wenozq(q_stencil, ql, qr)
            
                !call wenozq_gen(q_stencil, ql, qr, xi1, xi2)

                call righteigenvector(u_avg, v_avg, c_avg, h_avg, normalx, normaly, reigenvec)

                vecl = matmul(reigenvec, ql)!q(:, i, j)!matmul(reigenvec, ql)
                vecr = matmul(reigenvec, qr)!q(:, i+1, j)!matmul(reigenvec, qr)

                uleft(:, i, j) = vecl
                uright(:, i, j) = vecr

                ! do k = 1, 4
                !     if(isnan(vecl(k))) then
                !         print *, i, j, q_stencil
                !         stop
                !     endif
                ! enddo
            end do
        end do
        !$OMP END PARALLEL DO
    elseif(dir .eq. 2) then
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, yc, nx, ny, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc, ny-gc
            do i = gc+1, nx-gc
                !positive normals in the direction
                dxlen = (-x(i, j+1) + x(i+1, j+1))
                dylen = (y(i, j+1) - y(i+1, j+1))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                
                yiph = (y(i+1, j+1) + y(i, j+1))*half
                yimh = (y(i+1, j) + y(i, j))*half
                yip3h = (y(i, j+2) + y(i+1, j+2))*half

                yi1 = (yiph - yc(i, j))/(yiph - yimh)
                yi2 = (yiph - yc(i, j+1))/(yip3h - yiph)

                ! !calculation of roe or arithmetic averages
                u_avg = roeavg(rho(i,j), rho(i, j+1), u(i, j), u(i, j+1))
                v_avg = roeavg(rho(i,j), rho(i, j+1), v(i, j), v(i, j+1))
                h_avg = roeavg(rho(i,j), rho(i, j+1), h(i, j), h(i, j+1))
                c_avg = sqrt(g2*(h_avg - half*(u_avg**2 + v_avg**2)))

                call lefteigenvector(u_avg, v_avg, c_avg, normalx, normaly, leigenvec)

                !finding the relevant stencil and transforming it to characteristic variables
                do k = -2,3
                    do m = 1, 4
                        qcons(m) = q(m, i, j+k)
                    end do
                    vecm = matmul(leigenvec, qcons)
                    q_stencil(k, :) = vecm!matmul(qcons, leigenvec)
                end do

                !weno reconstruction procedure. the procedure can be found in toro 2005
                !other orders of reconstruction is also included in the module weno
                call weno5(q_stencil, ql, qr)
                ! call wenozq(q_stencil, ql, qr)
                !call wenozq_gen(q_stencil, ql, qr, yi1, yi2)

                call righteigenvector(u_avg, v_avg, c_avg, h_avg, normalx, normaly, reigenvec)

                vecl = matmul(reigenvec, ql)!q(:, i, j)!matmul(reigenvec, ql)
                vecr = matmul(reigenvec, qr)!q(:, i+1, j)!matmul(reigenvec, qr)

                uleft(:, i, j) = vecl
                uright(:, i, j) = vecr
            end do
        end do
        !$OMP END PARALLEL DO
    endif
    end subroutine weno_reconstruction

    subroutine hllc_fluxfunction(ql, qr, normx, normy, dlen, fluxvalue)
        implicit none

        real(prec), DIMENSION(4), intent(in) :: ql, qr
        real(prec), intent(in) :: normx, normy, dlen
        real(prec), DIMENSION(4), intent(out) :: fluxvalue
        real(prec), dimension(4) :: fl, fr, cstar
        real(prec) :: rhol, ul, vl, pl, rhor, ur, vr, pr
        real(prec) :: qnl, qnr, asoundl, asoundr, qro, aro, sl, sr, sm, pstar


        !force flux calculation as described in toro2005
        rhol = ql(1)
        ul = ql(2)/ql(1)
        vl = ql(3)/ql(1)
        pl = g2*(ql(4) - half*rhol*(ul**2 + vl**2))

        rhor = qr(1)
        ur = qr(2)/qr(1)
        vr = qr(3)/qr(1)
        pr = g2*(qr(4) - half*rhor*(ur**2 + vr**2))

        ! left and right flux calculation
        call fluxcalc(rhol, ul, vl, pl, normx, normy, fl)
        call fluxcalc(rhor, ur, vr, pr, normx, normy, fr)
        !calculating wave speed estimates for hll

        qnl = ul*normx + vl*normy
        qnr = ur*normx + vr*normy

        asoundl = dsqrt(g1*pl/rhol)
        asoundr = dsqrt(g1*pr/rhor)

        qro = roeavg(rhol, rhor, qnl, qnr)
        aro = roeavg(rhol, rhor, asoundl, asoundr)

        sl = dmin1(qnl-asoundl , qro-aro, zero)
        sr = dmax1(qro+aro , qnr+asoundr, zero)
        sm = (rhor*qnr*(sr - qnr) - rhol*qnl*(sl - qnl) + pl - pr)/(rhor*(sr - qnr) - rhol*(sl - qnl))

        pstar = rhol*(qnl - sl)*(qnl - sm) + pl

        if(sl .ge. zero) then
            fluxvalue = fl*dlen
        else if(sl .le. zero .and. sm .gt. zero) then
            cstar(1) = ql(1)*(sl - qnl)/(sl -sm)
            cstar(2) = (ql(2)*(sl - qnl) + (pstar - pl)*normx)/(sl - sm)
            cstar(3) = (ql(3)*(sl - qnl) + (pstar - pl)*normy)/(sl - sm)
            cstar(4) = (ql(4)*(sl - qnl) - pl*qnl + pstar*sm)/(sl - sm)

            fluxvalue(1) = cstar(1)*sm*dlen
            fluxvalue(2) = (cstar(2)*sm + pstar*normx)*dlen
            fluxvalue(3) = (cstar(3)*sm + pstar*normy)*dlen
            fluxvalue(4) = (cstar(4) + pstar)*sm*dlen
        else if(sm .le. zero .and. sr .ge. zero) then
            cstar(1) = qr(1)*(sr - qnr)/(sr -sm)
            cstar(2) = (qr(2)*(sr - qnr) + (pstar - pr)*normx)/(sr - sm)
            cstar(3) = (qr(3)*(sr - qnr) + (pstar - pr)*normy)/(sr - sm)
            cstar(4) = (qr(4)*(sr - qnr) - pr*qnr + pstar*sm)/(sr - sm)

            fluxvalue(1) = cstar(1)*sm*dlen
            fluxvalue(2) = (cstar(2)*sm + pstar*normx)*dlen
            fluxvalue(3) = (cstar(3)*sm + pstar*normy)*dlen
            fluxvalue(4) = (cstar(4) + pstar)*sm*dlen
        else if (sr .gt. zero) then
            fluxvalue = fr*dlen
        else
            print *, sl, sr, sm
            stop 'impossible scenario stopping the programx'
        end if

    end subroutine hllc_fluxfunction


    subroutine lf(ql, qr, normx, normy, dlen, fluxvalue)
        implicit none

        real(prec), DIMENSION(4), intent(in) :: ql, qr
        real(prec), intent(in) :: normx, normy, dlen
        real(prec), DIMENSION(4), intent(out) :: fluxvalue
        real(prec), dimension(4) :: fl, fr, cstar
        real(prec) :: rhol, ul, vl, pl, rhor, ur, vr, pr
        real(prec) :: qnl, qnr, asoundl, asoundr, qro, aro, sl, sr, sm, pstar


        !force flux calculation as described in toro2005
        rhol = ql(1)
        ul = ql(2)/ql(1)
        vl = ql(3)/ql(1)
        pl = g2*(ql(4) - half*rhol*(ul**2 + vl**2))

        rhor = qr(1)
        ur = qr(2)/qr(1)
        vr = qr(3)/qr(1)
        pr = g2*(qr(4) - half*rhor*(ur**2 + vr**2))

        ! left and right flux calculation
        call fluxcalc(rhol, ul, vl, pl, normx, normy, fl)
        call fluxcalc(rhor, ur, vr, pr, normx, normy, fr)
        !calculating wave speed estimates for hll

        qnl = ul*normx + vl*normy
        qnr = ur*normx + vr*normy

        asoundl = dsqrt(g1*pl/rhol)
        asoundr = dsqrt(g1*pr/rhor)

        qro = roeavg(rhol, rhor, qnl, qnr)
        aro = roeavg(rhol, rhor, asoundl, asoundr)

        sr = dmax1(qro+aro , qnr+asoundr, zero)

        fluxvalue = half*(fl + fr - sr*(qr - ql))*dlen
    end subroutine lf

    ! HLL FLUX CALCULTATION
    !SIMPLE MUSTA FLUX CALCULATION
    subroutine ihllc
        !loop iters and dir
        integer(i32) :: i, j, k
        integer(i8) :: dir, m
        !normals for the grid
        real(prec) :: normalx, normaly, dlen, dxlen, dylen
        !roe average or arithmatic avg values
        real(prec), dimension(4) :: ql1, ql2, qr1, qr2, qcons, vecl, vecr
        !stencil for the weno scheme
        real(prec), dimension(-2:2, 4) :: q_stencil1
        real(prec), dimension(-2:2, 4) :: q_stencil2
        ! flux variables
        real(prec), dimension(4) :: fluxval1, fluxval2, fluxval
        ! xdistance
        real(prec) :: xi1, xi2
        ! ROe average variables
        real(prec) :: u_avg, v_avg, c_avg, h_avg
        real(prec), dimension(4, 4) :: leigenvec, reigenvec
        
        uleft = zero
        uright = zero

        dir = 1 ! defines the direction of the cartesian co-ordinates 1 - x, 2 - y
        call weno_reconstruction(dir)
        !loop for the x-dir flux calculation
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, xc, yc, nx, ny, l_res, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc+1, ny-gc
            do i = gc, nx-gc
                !positive normals in the direction
                dxlen = -(x(i+1, j+1) - x(i+1, j))
                dylen = (y(i+1, j+1) - y(i+1, j))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen
                
                vecl = uleft(:, i, j)
                vecr = uright(:, i, j)

                call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
                ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)

                ! if ( i .lt. gc+3 .or. i .gt. nx-2*gc .or. j .lt. gc+3 .or. j .gt. ny-2*gc) then
                !     vecl = q(:, i, j)
                !     vecr = q(:, i+1, j)
                !     ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)
                !     call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
                ! else
                !     vecl = uleft(:, i, j)
                !     vecr = uright(:, i, j)
                    
                !     ! xi1 = gn1
                !     ! xi2 = gn2

                !     ! do k = -2,2
                !     !     do m = 1, 4
                !     !         qcons(m) = uleft(m, i+k, j)
                !     !     end do
                !     !     q_stencil1(k, :) = qcons!matmul(qcons, leigenvec)!qcons!matmul(qcons, leigenvec)
                !     ! end do
                    
                !     ! do k = -2,2
                !     !     do m = 1, 4
                !     !         qcons(m) = uright(m, i+k, j)
                !     !     end do
                !     !     q_stencil2(k, :) = qcons!matmul(qcons, leigenvec)!qcons!matmul(qcons, leigenvec)
                !     ! end do

                !     ! call wenozq_gaussl(q_stencil1, ql1, ql2, xi1, xi2)
                !     ! call wenozq_gaussl(q_stencil2, qr1, qr2, xi1, xi2)

                !     ! call hllc_fluxfunction(ql1, qr1, normalx, normaly, dlen, fluxval1)
                !     ! call hllc_fluxfunction(ql2, qr2, normalx, normaly, dlen, fluxval2)
                    
                !     ! fluxval = half*(fluxval1 + fluxval2)
                !     ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)
                !     call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
                ! endif
                l_res(:, i, j) = l_res(:, i, j) + fluxval
                l_res(:, i+1, j) = l_res(:, i+1, j) - fluxval

            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine


    !construction of j flux values
    subroutine jhllc
        
        !loop iters and dir
        integer(i32) :: i, j, k
        integer(i8) :: dir, m
        !normals for the grid
        real(prec) :: normalx, normaly, dlen, dxlen, dylen
        !roe average or arithmatic avg values
        real(prec), dimension(4) :: ql1, ql2, qr1, qr2, qcons, vecl, vecr
        !stencil for the weno scheme
        real(prec), dimension(-2:2, 4) :: q_stencil1
        real(prec), dimension(-2:2, 4) :: q_stencil2
        !dummy variables for flux calculation
        real(prec) :: rhol, rhor, rhom, ul, ur, um, vl, vr, vm, pl, pr, pm, qnl, qnr, asoundl, asoundr, sl, sr, sm, qro, aro, pstar
        !left and right flux values
        real(prec), dimension(4) :: fl, fr, cstar, fluxval1, fluxval2, fluxval
        ! xdistance
        ! real(prec) :: ximh, xiph, xip3h, xi1, xi2
        real(prec) :: yimh, yiph, yip3h, yi1, yi2
        !roe average variables
        real(prec) :: u_avg, v_avg, c_avg, h_avg
        real(prec), dimension(4, 4) :: leigenvec, reigenvec



        dir = 2 ! defines the direction of the cartesian co-ordinates 1 - x, 2 - y
        uleft = zero
        uright = zero
        call weno_reconstruction(dir)
        !loop for the y-dir flux calculation
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, xc, yc, nx, ny, l_res, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc, ny-gc
            do i = gc+1, nx-gc
                !positive normals in the direction
                dxlen = (-x(i, j+1) + x(i+1, j+1))
                dylen = (y(i, j+1) - y(i+1, j+1))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                vecl = uleft(:, i, j)
                vecr = uright(:, i, j)

                call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
                ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)

                ! if ( i .lt. gc+3 .or. i .gt. nx-2*gc .or. j .lt. gc+3 .or. j .gt. ny-2*gc) then
                !     vecl = q(:, i, j)
                !     vecr = q(:, i, j+1)

                !     ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)
                !     call lf(vecl, vecr, normalx, normaly, dlen, fluxval)

                ! else
                !     vecl = uleft(:, i, j)
                !     vecr = uright(:, i, j)
                    
                !     ! yi1 = gn1
                !     ! yi2 = gn2

                !     ! do k = -2,2
                !     !     do m = 1, 4
                !     !         qcons(m) = uleft(m, i, j+k)
                !     !     end do
                !     !     q_stencil1(k, :) = qcons!matmul(qcons, leigenvec)!qcons!matmul(qcons, leigenvec)
                !     ! end do
                    
                !     ! do k = -2,2
                !     !     do m = 1, 4
                !     !         qcons(m) = uright(m, i, j+k)
                !     !     end do
                !     !     q_stencil2(k, :) = qcons!matmul(qcons, leigenvec)!qcons!matmul(qcons, leigenvec)
                !     ! end do

                !     ! call wenozq_gaussl(q_stencil1, ql1, ql2, yi1, yi2)
                !     ! call wenozq_gaussl(q_stencil2, qr1, qr2, yi1, yi2)

                !     ! call hllc_fluxfunction(ql1, qr1, normalx, normaly, dlen, fluxval1)
                !     ! call hllc_fluxfunction(ql2, qr2, normalx, normaly, dlen, fluxval2)
                    
                !     ! fluxval = half*(fluxval1 + fluxval2)
                !     ! call hllc_fluxfunction(vecl, vecr, normalx, normaly, dlen, fluxval)
                !     call lf(vecl, vecr, normalx, normaly, dlen, fluxval)

                ! endif
                ! vecl = uleft(:, i, j)
                ! vecr = uright(:, i, j)
    
                l_res(:, i, j) = l_res(:, i, j) + fluxval
                l_res(:, i, j+1) = l_res(:, i, j+1) - fluxval
                
                ! if (i .eq. 4 .and. (j .eq. 52 .or. j .eq. 51)) then
                !     print *, 'very', j, l_res(:, i, j)
                !     print *, 'fluxval', j, fluxval
                ! endif
                ! g(:, i, j) = half*(fl + fr - sr*(vecr - vecl))*dlen
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine
!end of hllc flux calculationfluxvalroutine ihll
        ! HLL FLUX CALCULTATION
    !SIMPLE MUSTA FLUX CALCULATION
    subroutine ijst
        !loop iters and dir
        integer(i32) :: i, j, k
        !normals for the grid
        real(prec) :: normalx, normaly, dlen, dxlen, dylen
        ! flux variables
        real(prec), dimension(4) :: fluxval, fval
        ! face average values
        real(prec), dimension(4) :: wavg, diss
        ! Pressure sensor
        real(prec) :: p_etai, p_etaip1, p_eta
        ! gama value
        real(prec) :: gamai, gamaip1, gama
        ! epsilon values
        real(prec) ::  epsiph2,  epsiph4
        ! primitive variable
        real(prec) :: rhol, ul, vl, pl


    
        !loop for the x-dir flux calculation
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, xc, yc, nx, ny, l_res, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc+1, ny-gc
            do i = gc, nx-gc
                !positive normals in the direction
                dxlen = -(x(i+1, j+1) - x(i+1, j))
                dylen = (y(i+1, j+1) - y(i+1, j))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen
                
                wavg = half*(q(:, i, j) + q(:, i+1, j))

                rhol = wavg(1)
                ul = wavg(2)/wavg(1)
                vl = wavg(3)/wavg(1)
                pl = g2*(wavg(4) - half*rhol*(ul**2 + vl**2))
                
                call fluxcalc(rhol, ul, vl, pl, normalx, normaly, fval)

                p_etai = dabs(p(i+1, j) - two*p(i, j) + p(i-1, j))/(p(i+1, j) + two*p(i, j) + p(i-1, j))
                p_etaip1 = dabs(p(i+2, j) - two*p(i+1, j) + p(i, j))/(p(i+2, j) + two*p(i+1, j) + p(i, j))
                
                p_eta = half*(p_etai + p_etaip1)

                gamai = (dabs(u(i, j)*normalx + v(i, j)*normaly) + c(i, j))*dlen
                gamaip1 = (dabs(u(i+1, j)*normalx + v(i+1, j)*normaly) + c(i+1, j))*dlen
                gama = half*(gamai + gamaip1)

                epsiph2 = k2*(dmax1(p_etai, p_etaip1))
                epsiph4 = dmax1(zero, (k4 - epsiph2))
                
                diss = gama*(epsiph2*(q(:,i+1,j) - q(:,i,j)) - epsiph4*(q(:,i+2,j) -three*q(:,i+1,j) + three*q(:,i,j) -q(:,i-1,j)))
                
                fluxval = fval*dlen - diss
                
                l_res(:, i, j) = l_res(:, i, j) + fluxval
                l_res(:, i+1, j) = l_res(:, i+1, j) - fluxval

            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine


    !construction of j flux values
    subroutine jjst
        
        !loop iters and dir
        integer(i32) :: i, j, k
        !normals for the grid
        real(prec) :: normalx, normaly, dlen, dxlen, dylen
        ! flux variables
        real(prec), dimension(4) :: fluxval, fval
        ! face average values
        real(prec), dimension(4) :: wavg, diss
        ! Pressure sensor
        real(prec) :: p_etai, p_etaip1, p_eta
        ! gama
        real(prec) :: gamai, gamaip1, gama
        ! epsilon values
        real(prec) ::  epsiph2,  epsiph4
        ! primitive variables
        real(prec) :: rhol, ul, vl, pl

        !loop for the y-dir flux calculation
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(PRIVATE) &
        !$OMP SHARED(rho, u, v, p, c, h, q, x, y, xc, yc, nx, ny, l_res, uleft, uright) &
        !$OMP SHARED(g2, r)
        do j = gc, ny-gc
            do i = gc+1, nx-gc
                !positive normals in the direction
                dxlen = (-x(i, j+1) + x(i+1, j+1))
                dylen = (y(i, j+1) - y(i+1, j+1))
                dlen = dsqrt(dxlen**2 + dylen**2)
                normalx = dylen/dlen
                normaly = dxlen/dlen

                wavg = half*(q(:, i, j) + q(:, i, j+1))

                rhol = wavg(1)
                ul = wavg(2)/wavg(1)
                vl = wavg(3)/wavg(1)
                pl = g2*(wavg(4) - half*rhol*(ul**2 + vl**2))

                call fluxcalc(rhol, ul, vl, pl, normalx, normaly, fval)

                p_etai = dabs(p(i, j+1) - two*p(i, j) + p(i, j-1))/(p(i, j+1) + two*p(i, j) + p(i, j-1))
                p_etaip1 = dabs(p(i, j+2) - two*p(i, j+1) + p(i, j))/(p(i, j+2) + two*p(i, j+1) + p(i, j))
                
                p_eta = half*(p_etai + p_etaip1)

                gamai = (dabs(u(i, j)*normalx + v(i, j)*normaly) + c(i, j))*dlen
                gamaip1 = (dabs(u(i, j+1)*normalx + v(i, j+1)*normaly) + c(i, j+1))*dlen
                gama = half*(gamai + gamaip1)

                epsiph2 = k2*(dmax1(p_etai, p_etaip1))
                epsiph4 = dmax1(zero, (k4 - epsiph2))
                
                diss = gama*(epsiph2*(q(:,i,j+1) - q(:,i,j)) - epsiph4*(q(:,i,j+2) -three*q(:,i,j+1) + three*q(:,i,j) -q(:,i,j-1)))
                
                fluxval = fval*dlen - diss
                
                l_res(:, i, j) = l_res(:, i, j) + fluxval
                l_res(:, i, j+1) = l_res(:, i, j+1) - fluxval
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine
END MODULE
