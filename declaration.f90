!------------------------------------------------------------------------------
! Indian Institute of Technology - Madras, PhD scholar  
!------------------------------------------------------------------------------
!
! MODULE:  Declaration
!
!> @author
!> Vinoth P}
!
! DESCRIPTION: 
!>  Declaration of global variables for the program
!
! REVISION HISTORY:
! dd Mmm yyyy - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module declaration
    use constant
    implicit none
    
    ! Problem Initialisation Parameters
    logical :: cfl, ftime
    real(prec) :: cfl_no, dt, final_time
    integer(prec) :: max_counter, print_iter, print_val
    
    ! Gas dynamic constants
    real(prec) ::  g1, amol, pr_lam, c1ref
    real(prec) :: g2, r, c23, cpval, cvval

    ! Files for the problem
    character(len=30) :: gridfile, restartfile, outputfile
    logical :: restart_sim
    
    ! Counters initialised for the problem
    integer(i32) :: counter = 1
    integer(i32), parameter :: gc = 3

    ! variables for storing the grid
    integer(i32) :: imax, jmax
    integer(i32) :: nx, ny
    real(prec), allocatable, dimension(:, :) :: x, y
    real(prec), allocatable, dimension(:, :) :: xc, yc
    real(prec), allocatable, dimension(:, :) :: a, dl
    ! real(prec), allocatable, dimension(:, :) :: dx, dy
    
    ! Initila conditions for the problem
    real(prec) :: mach_no, uinit, vinit, pinit, tinit, rinit, einit, ainit
    logical :: ismach

    ! Problem specific initial conditions or boundary conditions
    real(prec) :: npr
    logical :: problem_specific
    
    !fluid properties
    real(prec), allocatable, dimension (:, :) :: rho, u, v, p, c, h, t

    !conservative and flux variables
    real(prec), allocatable, dimension(:, :, :) :: q

    !runge kutta integration variables
    real(prec), allocatable, dimension(:,:,:,:) :: q_rk
    real(prec), allocatable, dimension(:,:,:) :: l_res, l_vres

    !weno reconstruction variables
    real(prec),allocatable, dimension(:, :, :) :: uleft, uright

    !file pointers and file names
    character(len = 30) :: filename
    integer :: unit_id

    !time variable storage if needed
    real(prec), allocatable, dimension(:, :) :: dt_cell
    real(prec) :: time = zero, dtmin
    logical :: timeaccurate

    !source term
    real(prec), allocatable, dimension(:,:,:) :: qn, qi
    
    !Rk time marching step
    integer(i8) :: irkstep, nsteps
end module declaration