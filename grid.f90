module grid
    use constant
    use declaration

    implicit none

    contains
    subroutine readmesh
        integer:: i, j

        integer:: kmax, nblocks
        real(prec):: x1, x2, x3, x4, y1, y2, y3, y4
        real(prec):: dl1, dl2, dl3, dl4
        logical:: file_status = .true.

            !Reading the plot 3d Mesh file
        ! INQUIRE(file = gridfile, exist = file_status)
        IF (file_status) THEN
            OPEN(unit = unit_id, file = gridfile, form = 'formatted')
            ! READ(unit_id, *) nblocks
            ! READ(unit_id, *) imax, jmax, kmax

            nblocks = 1
            imax = 201
            jmax = 101
            kmax = 1


            PRINT *, 'No of grid points in the domain'
            PRINT *, imax, jmax
            PRINT *, 'No of cells in the domain'
            nx = (imax-1)+2*gc
            ny = (jmax-1)+2*gc
            PRINT *, imax-1, jmax-1
            PRINT *, 'No of cells including the ghost point'
            PRINT *, nx, ny
            PRINT *, 'Allocation of variables to the domain'

            CALL allocate_variables(imax, jmax)

            ! READ(unit_id, *) ((x(i+gc, j+gc), i = 1, imax), j = 1, jmax)
            ! READ(unit_id, *) ((y(i+gc, j+gc), i = 1, imax), j = 1, jmax)
            do j = 1, jmax
                do i = 1, imax
                    x(i+gc, j+gc) = (i-1)*(1.0d0/imax)
                    y(i+gc, j+gc) = (j-1)*(0.5d0/jmax)
                enddo
            enddo
            CLOSE(unit_id)
        ELSE
            PRINT *, "File not found aborting the program"
            STOP
        ENDIF

        ! Creation of ghost cell parameters 
        DO j = gc+1, ny-gc
            ! left boundary cells
            x(gc, j) = two*x(gc+1, j) - x(gc+2, j) 
            x(gc-1, j) = two*x(gc+1-1, j) - x(gc+2-1, j) 
            x(gc-2, j) = two*x(gc+1-2, j) - x(gc+2-2, j)
            
            y(gc, j) = two*y(gc+1, j) - y(gc+2, j) 
            y(gc-1, j) = two*y(gc+1-1, j) - y(gc+2-1, j) 
            y(gc-2, j) = two*y(gc+1-2, j) - y(gc+2-2, j)

            ! right boundary cells
            x(nx-gc+1, j) = two*x(nx-gc, j) - x(nx-gc-1, j)
            x(nx-gc+2, j) = two*x(nx-gc+1, j) - x(nx-gc, j)
            x(nx-gc+3, j) = two*x(nx-gc+2, j) - x(nx-gc+1, j)

            y(nx-gc+1, j) = two*y(nx-gc, j) - y(nx-gc-1, j)
            y(nx-gc+2, j) = two*y(nx-gc+1, j) - y(nx-gc, j)
            y(nx-gc+3, j) = two*y(nx-gc+2, j) - y(nx-gc+1, j)
        ENDDO

        DO i = gc+1, nx-gc
            ! bottom boundary cells
            x(i, gc) = two*x(i, gc+1) - x(i, gc+2) 
            x(i, gc-1) = two*x(i, gc+1-1) - x(i, gc+2-1) 
            x(i, gc-2) = two*x(i, gc+1-2) - x(i, gc+2-2)
            
            y(i, gc) = two*y(i, gc+1) - y(i, gc+2) 
            y(i, gc-1) = two*y(i, gc+1-1) - y(i, gc+2-1) 
            y(i, gc-2) = two*y(i, gc+1-2) - y(i, gc+2-2)

            ! top boundary cells
            x(i, ny-gc+1) = two*x(i, ny-gc) - x(i, ny-gc-1)
            x(i, ny-gc+2) = two*x(i, ny-gc+1) - x(i, ny-gc)
            x(i, ny-gc+3) = two*x(i, ny-gc+2) - x(i, ny-gc+1)

            y(i, ny-gc+1) = two*y(i, ny-gc) - y(i, ny-gc-1)
            y(i, ny-gc+2) = two*y(i, ny-gc+1) - y(i, ny-gc)
            y(i, ny-gc+3) = two*y(i, ny-gc+2) - y(i, ny-gc+1)
        ENDDO

        !calculation of cell center
        ! $omp parallel do reduction(min : dl)
        DO j = 1, ny
            DO i = 1, nx
                ! xc(i, j) = quarter*(x(i, j) + x(i+1, j) + x(i+1, j+1) + x(i, j+1))
                ! yc(i, j) = quarter*(y(i, j) + y(i+1, j) + y(i+1, j+1) + y(i, j+1))

                ! dx(i, j) = x(i+1, j) - x(i, j)
                ! dy(i, j) = y(i, j+1) - y(i, j)

                x1 = x(i, j)
                x2 = x(i+1, j)
                x3 = x(i+1, j+1)
                x4 = x(i, j+1)

                y1 = y(i, j)
                y2 = y(i+1, j)
                y3 = y(i+1, j+1)
                y4 = y(i, j+1)

                xc(i, j) = quarter*(x1+x2+x3+x4)
                yc(i, j) = quarter*(y1+y2+y3+y4)

                dl1 = SQRT((x2-x1)**2 + (y2-y1)**2)
                dl2 = SQRT((x3-x2)**2 + (y3-y2)**2)
                dl3 = SQRT((x4-x3)**2 + (y4-y3)**2)
                dl4 = SQRT((x1-x4)**2 + (y1-y4)**2)

                dl(i, j) = DMIN1(dl1, dl2, dl3, dl4)

                a(i, j) = half*((x3-x1)*(y4-y2) - (x4-x2)*(y3-y1))
            ENDDO
        ENDDO

        ! DO j = gc+1, ny-gc
        !     xc(gc, j) = xc(gc+1, j) - (x(gc+2, j) - x(gc+1, j))  ! dx(gc+1, j)
        !     yc(gc, j) = yc(gc+1, j) - (y(gc+2, j) - y(gc+1, j))
        !     xc(nx-gc+1, j) = xc(nx-gc, j) - (x(nx-gc+1, j) - x(nx-gc, j))  ! dx(nx-gc, j)
        !     yc(nx-gc+1, j) = yc(nx-gc, j) - (y(nx-gc+1, j) - y(nx-gc, j)) 
        ! ENDDO

        ! DO i = gc+1, nx-gc
        !     xc(i, gc) = xc(i, gc+1)
        !     yc(i, gc) = yc(i, gc+1) - dy(i, gc+1)
        !     xc(i, ny-gc+1) = xc(i, ny-gc)
        !     yc(i, ny-gc+1) = yc(i, ny-gc) + dy(i, ny-gc)
        ! ENDDO
    end subroutine readmesh


    subroutine allocate_variables(ixmax, jymax)
    implicit none
    INTEGER, INTENT(IN):: ixmax, jymax
    INTEGER:: iostat

    ALLOCATE(x(ixmax+2*gc, jymax+2*gc), y(ixmax+2*gc, jymax+2*gc), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the points'

    ALLOCATE(xc(ixmax+2*gc-1, jymax+2*gc-1), yc(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the cell centers'

    ALLOCATE(a(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the areas'

    ALLOCATE(dl(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the areas'

    ALLOCATE(dt_cell(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the timestep'

    ALLOCATE(rho(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the density'

    ALLOCATE(u(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the x velocity'

    ALLOCATE(v(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the yvelocity'

    ALLOCATE(p(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the pressure'

    ALLOCATE(t(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the pressure'

    ALLOCATE(c(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the sound speed'

    ALLOCATE(h(ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the enthalpy'

    ALLOCATE(q(4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the Conservative variables'

    ALLOCATE(uleft(4, ixmax+2*gc-1, jymax+2*gc-1), uright(4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the weno reconstruction variables'

    ALLOCATE(l_res(4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the residual variables'
    
    ALLOCATE(l_vres(4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the residual variables'

    ALLOCATE(q_rk(4, 4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the timestepping variables'

    ALLOCATE(qn(4, ixmax+2*gc-1, jymax+2*gc-1), qi(4, ixmax+2*gc-1, jymax+2*gc-1), STAT = iostat)
    IF (iostat /= 0) STOP 'Error in allocating the source variables'

    !initialisation of variables
    x = zero
    y = zero
    xc = zero
    dl = zero
    a = zero
    q = zero
    q_rk = zero
    l_res = zero
    l_vres = zero
    ! f = zero
    ! g = zero

    rho = zero
    u = zero
    v = zero
    p = zero
    c = zero
    h = zero
    t = zero

    end subroutine allocate_variables

    subroutine deallocation
        implicit none
        
        DEALLOCATE(x, y, xc, yc, a, dl)
        DEALLOCATE(rho, u, v, p, t, c, h)
        DEALLOCATE(q, uleft, uright, l_res, l_vres, q_rk, qn, qi)
        
    end subroutine deallocation
    
end module grid


