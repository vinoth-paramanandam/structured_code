module odesolver
    use constant
    use declaration

    implicit none

    contains

    subroutine rksolver(nstep, rkstep)
    implicit none
    integer(i8), intent(in) :: nstep, rkstep
    integer(i32) :: i, j
    integer(i8) :: k
    real(prec) :: dtmod

    if (nstep .eq. 4) then

        if(rkstep .eq. 1) then
            !$omp parallel do &
            !$omp default(none) &
            !$omp private(i, j, k, dtmod) &
            !$omp shared(nx, ny, a, dt_cell, q_rk, l_res, l_vres, qn, qi)
            do j = gc+1, ny-gc
                do i = gc+1, nx-gc
                    dtmod = dt_cell(i, j)/a(i,j)
                    ! print *, i, j, dt_cell(i,j)
                    do k = 1,4
                        q_rk(1, k, i, j) = dtmod*(-l_res(k, i, j) + l_vres(k,i,j))
                        qn(k, i, j) = qi(k, i, j) + half*q_rk(1, k, i, j)
                    end do
                end do
            end do
            !$omp end parallel do
        elseif(rkstep .eq. 2) then
            !$omp parallel do &
            !$omp default(none) &
            !$omp private(i, j, k, dtmod) &
            !$omp shared(nx, ny, a, dt_cell, l_res, l_vres, q_rk, qn, qi)
            do j = gc+1, ny-gc
                do i = gc+1, nx-gc
                    dtmod = dt_cell(i, j)/a(i,j)
                    do k = 1,4
                        q_rk(2, k, i, j) = dtmod*(-l_res(k, i, j) + l_vres(k,i,j))
                        qn(k, i, j) = qi(k, i, j) + half*q_rk(2, k, i, j)
                        ! q_rk(3, k, i, j) = half*q_rk(1, k, i, j) + half*q_rk(2, k, i, j) + half*dtmod*(l_res(k, i, j))
                        ! if( i .eq. 7 .and. j .eq. 103) then
                        !     print *, qn(k, i, j), qi(k, i, j), l_res(k, i, j)
                        ! endif 
                    end do
                end do
            end do
            !$omp end parallel do
        elseif(rkstep .eq. 3) then
            !$omp parallel do &
            !$omp default(none) &
            !$omp private(i, j, k, dtmod) &
            !$omp shared(nx, ny, a, q_rk, dt_cell, l_res, l_vres, qn, qi)
            do j = gc+1, ny-gc
                do i = gc+1, nx-gc
                    dtmod = dt_cell(i, j)/a(i,j)
                    do k = 1,4
                        q_rk(3, k, i, j) = dtmod*(-l_res(k, i, j) + l_vres(k,i,j))
                        qn(k, i, j) = qi(k, i, j) + q_rk(3, k, i, j)
                        ! q_rk(3, k, i, j) = half*q_rk(1, k, i, j) + half*q_rk(2, k, i, j) + half*dtmod*(l_res(k, i, j))
                    end do
                end do
            end do
            !$omp end parallel do
        elseif(rkstep .eq. 4) then
            !$omp parallel do &
            !$omp default(none) &
            !$omp private(i, j, k, dtmod) &
            !$omp shared(nx, ny, a, q_rk, dt_cell, l_res, l_vres, qn, qi)
            do j = gc+1, ny-gc
                do i = gc+1, nx-gc
                    dtmod = dt_cell(i, j)/a(i,j)
                    do k = 1,4
                        q_rk(4, k, i, j) = dtmod*(-l_res(k, i, j) + l_vres(k,i,j))
                        qn(k, i, j) = qi(k, i, j) + &
                        & (q_rk(1, k, i, j) + two*(q_rk(2, k, i, j) + q_rk(3, k, i, j)) + q_rk(4, k, i, j))/six
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    endif
    end subroutine rksolver

end module odesolver