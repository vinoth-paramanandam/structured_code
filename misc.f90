module misc
    use constant
    use declaration
    use grid

    implicit none

    contains

    subroutine prm2cons
        integer(prec) :: i, j
        !converting the flow variables to conservative variables
        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i, j) &
        !$omp shared(rho, u, v, p, q, nx, ny, g2)
        do j = gc+1, ny-gc
            do i = gc+1, nx-gc
                q(1, i, j) = rho(i, j)
                q(2, i, j) = rho(i, j)*u(i, j)
                q(3, i, j) = rho(i, j)*v(i, j)
                q(4, i, j) = p(i, j)/g2 + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2)
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine cons2prm
        integer(prec) :: i, j
        !$omp parallel do &
        !$omp default(none) &
        !$omp private(i, j) &
        !$omp shared(rho, u, v, p, c, h, t, q, nx, ny, xc, yc, g2, g1, r, counter)
        do j = gc+1, ny-gc
            do i = gc+1, nx-gc
                rho(i, j) = q(1, i, j)
                u(i, j) = q(2, i, j)/rho(i, j)
                v(i, j) = q(3, i, j)/rho(i, j)
                p(i, j) = (q(4, i, j) - half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))*g2
                c(i, j) = dsqrt(g1*p(i, j)/rho(i, j))
                h(i, j) = (p(i, j)/g2 + p(i, j) + half*rho(i, j)*(u(i, j)**2 + v(i, j)**2))/rho(i,j)
                t(i,j) = p(i,j)/(R*rho(i,j))
                if(isnan(rho(i,j))) then
                    print *, counter, i, j, xc(i, j), yc(i,j)
                    print *, q(:, i, j)
                    STOP
                endif
                ! if(isnan(p(i,j)) .or. isnan(rho(i,j)) .or. isnan(c(i,j))) then
                !     print *, q(:,i,j), i, j, xc(i, j), yc(i, j), p(i,j), rho(i,j)
                !     stop
                ! endif
            end do
        end do
        !$omp end parallel do
    end subroutine

    !function to calculate the roe averages
    real(prec) function roeavg(r1, r2, u1, u2)
    real(prec), intent(in) :: r1, r2, u1, u2
    real(prec) :: d1, d2

    ! d = sqrt(r2/r1)
    ! roeavg = (u1 + d*u2)/(1 + d)
    d1 = dsqrt(r1)
    d2 = dsqrt(r2)
    roeavg = (u1*d1 + u2*d2)/(d1 + d2)
    end function

    !calculation of leigenvec given nx, ny, and roe averages
    subroutine lefteigenvector(u_, v_, c_, normalx, normaly, leigenvec)
        real(prec), intent(in) :: u_, v_, c_, normalx, normaly
        real(prec) :: un, ut, ek ! normal values of velocities will be used in calc
        real(prec), dimension(4, 4), intent(out) :: leigenvec

        !calculating the normal velocities
        un = u_*normalx + v_*normaly
        ut = v_*normalx - u_*normaly
        ek = half*(u_**2 + v_**2) ! kinetic energy

        leigenvec(1, 1) = half*(g2*ek + c_*un)/c_**2
        leigenvec(1, 2) = -half*(g2*u_ + c_*normalx)/c_**2
        leigenvec(1, 3) = -half*(g2*v_ + c_*normaly)/c_**2
        leigenvec(1, 4) = half*g2/c_**2

        leigenvec(2, 1) = (c_**2 - g2*ek)/c_**2
        leigenvec(2, 2) = (g2*u_)/c_**2
        leigenvec(2, 3) = (g2*v_)/c_**2
        leigenvec(2, 4) = -g2/c_**2

        leigenvec(3, 1) = ut
        leigenvec(3, 2) = normaly
        leigenvec(3, 3) = -normalx
        leigenvec(3, 4) = zero

        leigenvec(4, 1) = half*(g2*ek - c_*un)/c_**2
        leigenvec(4, 2) = half*(c_*normalx - g2*u_)/c_**2
        leigenvec(4, 3) = half*(c_*normaly - g2*v_)/c_**2
        leigenvec(4, 4) = half*g2/c_**2
    end subroutine

    subroutine righteigenvector(u_, v_, c_, h_, normalx, normaly, reigenvec)
        real(prec), intent(in) :: u_, v_, c_, h_, normalx, normaly
        real(prec) :: un, ut, h0, ek
        real(prec), dimension(4, 4) :: reigenvec

        !calculation of normal velocities
        un = u_*normalx + v_*normaly
        ut = v_*normalx - u_*normaly
        ek = half*(u_**2 + v_**2)
        h0 = c_**2/g2 + ek

        reigenvec(1, 1) = 1
        reigenvec(2, 1) = u_ - c_*normalx
        reigenvec(3, 1) = v_ - c_*normaly
        reigenvec(4, 1) = h0 - c_*un

        reigenvec(1, 2) = 1
        reigenvec(2, 2) = u_
        reigenvec(3, 2) = v_
        reigenvec(4, 2) = ek

        reigenvec(1, 3) = 0
        reigenvec(2, 3) = normaly
        reigenvec(3, 3) = -normalx
        reigenvec(4, 3) = -ut

        reigenvec(1, 4) = 1
        reigenvec(2, 4) = u_ + c_*normalx
        reigenvec(3, 4) = v_ + c_*normaly
        reigenvec(4, 4) = h0 + c_*un
    end subroutine

    subroutine fluxcalc(rhof, uf, vf, pf, normalx, normaly, f)
        real(prec), intent(in) :: rhof, uf, vf, pf, normalx, normaly
        real(prec) :: un
        real(prec), dimension(4), intent(out) :: f
        un = uf*normalx + vf*normaly
        f(1) = rhof*un
        f(2) = rhof*uf*un + pf*normalx
        f(3) = rhof*vf*un + pf*normaly
        f(4) = (pf/g2 + half*rhof*(uf**2 + vf**2) + pf)*un
    end subroutine

    subroutine eigenval(ur, vr, cr, normalx, normaly, el1, el2, el3, el4)
        real(prec), intent(in) :: ur, vr, cr, normalx, normaly
        real(prec) :: un
        real(prec), intent(out) :: el1, el2, el3, el4

        un = ur*normalx + vr*normaly
        el1 = un - cr
        el2 = un
        el3 = un
        el4 = un + cr
    end subroutine

    subroutine timestep
        integer :: i, j
        real(prec) :: qplus, asound
        if(cfl) then
            dtmin = 1.0d8
            !$omp parallel do reduction (min:dtmin) &
            !$omp default(none) &
            !$omp private(i, j, qplus, asound) &
            !$omp shared(u, v, p, rho, g1, cfl_no, dl, dt_cell, nx, ny)
            do j = gc+1, ny-gc
                do i = gc+1, nx-gc
                    asound = dsqrt(g1*p(i, j)/rho(i, j))
                    qplus = dsqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j)) + asound
                    dt_cell(i, j) = cfl_no*dl(i,j)/qplus
                    dtmin = dmin1(dtmin, dt_cell(i, j))
                enddo
            enddo
            !$omp end parallel do

            if(timeaccurate) then
                dt_cell = dtmin
            endif
        else
            dt_cell = dt
            dtmin = dt
        endif
       
    end subroutine timestep

    real(prec) function minmod(a, b)
        real(prec), intent(in) :: a, b
        !real(prec), intent(out) :: minmod

        if((a > zero .and. b > zero) .or. (a < zero .and. b < zero)) then
            if(abs(a) > abs(b)) then
                minmod = b
            else
                minmod = a
            end if
        else
            minmod = zero
        end if
    end function

    ! real(prec) function psi(x)
    !     real(prec), intent(in) :: x
    !     !real(prec), intent(out) :: psi

    !     psi = sqrt(dd + x**2)
    ! end function


    subroutine output()
        integer(prec) :: i, j

        open (unit = unit_id, status = 'replace', file = filename, form = 'formatted')
        ! !headers for tecplot file
        ! write (unit_id, *) "variables = x, y, rho, u ,v, p, m"
        ! write (unit_id, *) "zone j=", ny-2*gc, "i=",nx-2*gc, "f=point"

        ! !writing the data
        ! do j = gc+1, ny-gc
        !     do i = gc+1, nx-gc
        !         asound = dsqrt(g1*p(i, j)/rho(i, j))
        !         M_num = dsqrt(u(i, j)*u(i, j) + v(i, j)*v(i, j))/asound
        !         write(unit_id, '(6f15.5)') xc(i, j), yc(i, j), rho(i, j), u(i, j), v(i, j), p(i, j), M_num
        !     end do
        ! end do
        
        !headers for tecplot files
        write (unit_id, *) "VARIABLES = x, y, rho, u, v, p, t"
        write (unit_id, *) "ZONE I =", imax, ",J =", jmax, ", DATAPACKING=BLOCK, VARLOCATION=([3,4,5,6,7]=CELLCENTERED)"
        write (unit_id, *) "SOLUTIONTIME=", time
        
        do j = 1, jmax
            do i = 1, imax
                write(unit_id, *) x(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax
            do i = 1, imax
                write(unit_id, *) y(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax-1
            do i = 1, imax-1
                write(unit_id, *) rho(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax-1
            do i = 1, imax-1
                write(unit_id, *) u(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax-1
            do i = 1, imax-1
                write(unit_id, *) v(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax-1
            do i = 1, imax-1
                write(unit_id, *) p(gc+i, gc+j)
            enddo
        enddo
        do j = 1, jmax-1
            do i = 1, imax-1
                write(unit_id, *) t(gc+i, gc+j)
            enddo
        enddo
        ! write(unit_id, *) ((x(gc+i,gc+j), i = 1, imax), j= 1, jmax)
        ! write(unit_id, *) ((y(gc+i,gc+j), i = 1, imax), j= 1, jmax)
        ! write(unit_id, *) ((rho(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        ! write(unit_id, *) ((u(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        ! write(unit_id, *) ((v(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        ! write(unit_id, *) ((p(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        ! write(unit_id, *) ((t(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        ! ! write(unit_id, *) ((t(gc+i,gc+j), i = 1, imax-1), j= 1, jmax-1)
        
        close (unit_id)
    end subroutine output

    subroutine restart()
        integer :: i, j

        open(unit = 11, status='replace', file=restartfile, form='formatted')
        write(11, *) ((q(1,i,j), i = gc+1, nx-gc), j= gc+1, ny-gc)
        write(11, *) ((q(2,i,j), i = gc+1, nx-gc), j= gc+1, ny-gc)
        write(11, *) ((q(3,i,j), i = gc+1, nx-gc), j= gc+1, ny-gc)
        write(11, *) ((q(4,i,j), i = gc+1, nx-gc), j= gc+1, ny-gc)
        close(11)
    end subroutine restart
end module misc
