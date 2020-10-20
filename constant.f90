!------------------------------------------------------------------------------
! Indian Institute of Technology, PhD scholar
!------------------------------------------------------------------------------
!
! MODULE: Constant
!
!> @author
!> Vinoth P}
!
! DESCRIPTION: 
!> Constants used in the file and the precision is defined in this module
!
! REVISION HISTORY:
! dd Mmm yyyy - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module constant
    use iso_fortran_env
    implicit none
    
    ! sp - single precision 
    ! dp - double precision 
    ! precision for the real numbers in the program

    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

    ! precision for the integer number in the program
    integer, parameter :: i8 = int8
    integer, parameter :: i16 = int16
    integer, parameter :: i32 = int32
    integer, parameter :: i64 = int64

    ! selecting the precision for the real numbers
    integer, parameter :: prec = dp

    ! some constants for the program
    real(prec), parameter :: pi = 3.14159265358979311599796346854d0
    real(prec), parameter :: Runiv = 8.3144621d3

    ! some real number for the calculation
    real(prec), parameter :: one = 1.0d0
    real(prec), parameter :: two = 2.0d0
    real(prec), parameter :: three = 3.0d0
    real(prec), parameter :: four = 4.0d0
    real(prec), parameter :: five = 5.0d0
    real(prec), parameter :: six = 6.0d0
    real(prec), parameter :: seven = 7.0d0
    real(prec), parameter :: eight = 8.0d0
    real(prec), parameter :: nine = 9.0d0
    real(prec), parameter :: zero = 0.0d0

    real(prec), parameter :: ten = 10.0d0
    real(prec), parameter :: eleven = 11.0d0

    real(prec), parameter :: half = 0.5d0
    real(prec), parameter :: quarter = 0.25d0
    real(prec), parameter :: threefourths = 0.75d0
    real(prec), parameter :: sixth = one/six

    !some small numbers used in various schemes
    real(prec), parameter :: e_weno = 1.0d-6
    real(prec), parameter :: e_tvd = 1.0d-7
    real(prec), parameter :: d_tvd = 0.2d0
    !real(prec), parameter :: dd = 0.0625_prec

    !weights in weno schemes
    !wenojs weights
    real(prec), parameter :: c0_js = 0.1d0, c1_js = 0.6d0, c2_js = 0.3d0
    ! real(prec), parameter :: c0 = 0.2d0, c1 = 0.4d0, c2 = 0.4d0
    real(prec), parameter :: c0 = 0.98d0, c1 = 0.01d0, c2 = 0.01d0
    real(prec), parameter :: a21 = one/12.0d0, a31 = three/(20.0d0), a41 = three/14.0d0, a42 = three/(560.0d0)
    real(prec), parameter :: gn1 = half/dsqrt(three), gn2 = -half/dsqrt(three)
    real(prec), parameter :: k_is = 13.0_prec/12.0_prec
    real(prec), parameter :: k2 = half, k4 = 1.0d0/128.0d0
end module constant 
