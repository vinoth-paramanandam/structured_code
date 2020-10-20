program main
    use constant
    use declaration
    use grid
    use misc
    use boundary
    use flux
    use odesolver
    implicit none

    integer(i32):: i, j
    logical:: file_exist, maxcount = .false., l2norm = .false., norm
    real(prec):: sum1n, sum2n, temp, rss, diffro

    namelist/pbl_params/  cfl, ftime, timeaccurate, cfl_no, dt, final_time, max_counter, print_iter, print_val, norm
    namelist/const_params/  g1, amol, pr_lam, c1ref
    namelist/file_params/gridfile, restartfile, outputfile, restart_sim
    namelist/initial_params/mach_no, ismach, uinit, vinit, pinit, tinit
    namelist/specific_conds/problem_specific, npr
    namelist/rk_params/nsteps

    open(unit = 11, file = 'input.in')
    read(11, nml = pbl_params)
    read(11, nml = const_params)
    read(11, nml = file_params)
    read(11, nml = initial_params)
    read(11, nml = specific_conds)
    read(11, nml = rk_params)
    close(11)

    r = Runiv/amol
    g2 = g1-1.0d0
    cpval = g1*r/g2
    cvval = r/g2
    c23 = two/three

    if(ismach) then
        ainit = dsqrt(g1*r*tinit)
        uinit = mach_no*ainit
    endif

    rinit = pinit/(r*tinit)
    einit = cvval*tinit+half*(uinit*uinit+vinit*vinit)
    print *, pinit, einit, tinit

    ! Geomentry routine to read the grid file
    call readmesh

    ! Calling the initial condition routine
    if(restart_sim) then
        inquire(file = restartfile, exist = file_exist)
        if (file_exist) then
            print *, "Reading the restart file for loading initial conditions"
            open(unit = 10, file = restartfile, form = 'formatted')
            read(10, *) ((q(1, i, j), i = gc+1, nx-gc), j = gc+1, ny-gc)
            read(10, *) ((q(2, i, j), i = gc+1, nx-gc), j = gc+1, ny-gc)
            read(10, *) ((q(3, i, j), i = gc+1, nx-gc), j = gc+1, ny-gc)
            read(10, *) ((q(4, i, j), i = gc+1, nx-gc), j = gc+1, ny-gc)
            close(10)
        else
            stop 'Restart file not found. Exiting the program'
        endif
    else
        print *, 'Applying the initial copnditions'
        print *, rinit, einit
        do j = gc+1, ny-gc
            do i = gc+1, nx-gc
                q(1, i, j) = rinit
                q(2, i, j) = zero  ! rinit*uinit
                q(3, i, j) = zero  ! rinit*vinit
                q(4, i, j) = rinit*einit  ! pinit/g2+half*rinit*(uinit**2 + vinit**2)
                
            end do
        end do
    endif

    ! COnverting the conservative variable to primitive
    call cons2prm

    filename = 'output.dat'
    unit_id = 4
    call output

    open (unit = 3, status = 'replace', file = outputfile, form = 'formatted', position = 'append')
    write (3, *)"variables = counter, dt, l1norm, l2norm"
    write (3, *)"zone t = onlyzone i =", max_counter+1, "f = point"

    ! Apply boundary conditions
    call nozzlebc

    do
        ! Copying the solution to the new variable
        qi = q

        ! Calling the time step to calculate the dt
        call timestep
        ! time is calculated using dtmin from time step routine
        time = time+dtmin

        ! Calling the viscous routine if needed
        ! call laminar_viscous
        l_vres = zero
        ! RK timestepping for time marching
        do irkstep = 1, nsteps
            ! print *, irkstep
            ! if(irkstep .eq. 2) stop
            l_res = zero
           
            call ihllc
            call jhllc
            

            call rksolver(nsteps, irkstep)

            q = qn

            call cons2prm
            
            call nozzlebc
        enddo
    
    rss = 1.0d8
    diffro = zero
    sum1n = zero
    sum2n = zero
    
    if(mod(counter, print_val) .eq. 0) then
        if(norm) then
        !$omp parallel do reduction (+:sum1n) &
        !$omp default(none) &
        !$omp private(i, j) &
        !$omp shared(q, qi, nx, ny)
        do j = gc+1, ny-gc
            do i = gc+1, nx-gc
                sum1n = sum1n+dabs(q(1, i, j) - qi(1, i, j))
            enddo
        enddo

        !$omp parallel do reduction (+:sum2n) &
        !$omp default(none) &
        !$omp private(i, j, temp) &
        !$omp shared(q, qi, nx, ny)
        do j = gc+1, ny-gc
            do i = gc+1, nx-gc
                temp =  q(1, i, j) - qi(1, i, j)
                sum2n = sum2n+temp*temp
            enddo
        enddo

        rss = dsqrt(sum2n)

        open(unit = 3, status = 'old', file = outputfile, form = 'formatted', position = 'append')
        write(3, '(1i8, 5es14.7)') counter, time,  dtmin, dsqrt(sum1n), rss, npr
        write (*, '(1i8, 3es14.7)') counter, dsqrt(sum1n), rss, npr

        endif
   
    endif

    !if(mod(counter, print_iter) .eq. 0) npr = npr+0.05d0
    
    if (counter >= max_counter) maxcount = .true.
    if (rss .lt. 1.0d-8) l2norm = .true.

    if (maxcount .OR. ftime .OR. l2norm) then
        write(filename, '(a, i8.8, a)') "plot",counter, ".dat"
        unit_id = counter+100
        call output
        call restart
        exit
    else
        !condition for printing the file
        if (mod(counter, print_iter) .eq. 0) then
            write(filename, '(a, i8.8, a)') "plot",counter, ".dat"
            unit_id = 100+counter
            call output
            call restart
        end if
        !incrementing the counter
        counter = counter+1
    end if
    
    enddo

    ! Deallocating the allocated memory
    call deallocation
end program main
