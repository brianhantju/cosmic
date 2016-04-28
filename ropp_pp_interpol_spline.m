function [ f_int, fd_int, fd2_int ] = ropp_pp_interpol_spline( x, f, d2, x_int )
%ROPP_PP_INTERPOL_SPLINE
% !****s* PPSpline/ropp_pp_interpol_spline *
% !
% ! NAME
% !    ropp_pp_interpol_spline - Spline interpolation of a gridded function with
% !                              linear extrapolation outside grid extent
% !
% ! SYNOPSIS  
% !    call ropp_pp_interpol_spline(x, f, d2, x_int, f_int, fd_int, fd2_int)
% !
% ! DESCRIPTION
% !    Search for grid interval containing the interpolation point and summation
% !    of polynomail with given spline coefficients
% !
% !****
%   SUBROUTINE ropp_pp_interpol_spline(x, f, d2, x_int, f_int, fd_int, fd2_int)
    
%     ! 3.1 Declarations

%     REAL(wp), DIMENSION(:), INTENT(in)  :: x  ! Argument grid (monotonous)
%     REAL(wp), DIMENSION(:), INTENT(in)  :: f  ! Gridded function
%     REAL(wp), DIMENSION(:), INTENT(in)  :: d2 ! 2nd derivative of spline
%     REAL(wp),               INTENT(in)  :: x_int  ! Interpolation point
%     REAL(wp),               INTENT(out) :: f_int  ! Interpolated function value
%     REAL(wp), OPTIONAL,     INTENT(out) :: fd_int ! Interpolated 1st derivative
%     REAL(wp), OPTIONAL,     INTENT(out) :: fd2_int ! Interpolated 2nd deriv
    
%     REAL(wp) :: a1, a2, a3   ! Polynomial coefficients
%     REAL(wp) :: dx           ! Grid interval
%     REAL(wp) :: dx_t         ! Grid-point to interpolation-point distance
%     REAL(wp) :: x_t          ! Interpolation point projected to grid extent
%     REAL(wp) :: fd           ! Interpolated derivative
%     INTEGER  :: i            ! Array index
%     INTEGER  :: N            ! Number of data
%     INTEGER  :: i_int        ! Interpolation interval index

%     ! 3.2 Location of interpolation point inside grid

    N     = length(x);
    
    x_t   = min(max(x_int, min(x(1),x(N))), max(x(1),x(N)));
    
    i_int = ropp_pp_seek_index(x, x_t);
    i     = max(i_int, 1);
    
%     ! 2.3 Calculation of interpolation coefficients

    dx    = x(i+1) - x(i);
    a2    = d2(i)/2 ; % original 2.0_wp
    a3    = (d2(i+1) - d2(i))/(6 * dx) ; % original 6.0_wp
    a1    = (f(i+1) - f(i))/dx - dx * (a2 + dx * a3);
    
%     ! 3.4 Calculated interpolated value

    dx_t  = x_t - x(i);
    fd    = a1 + dx_t*(2*a2 + dx_t*3*a3);
    f_int = f(i) + dx_t*(a1 + dx_t*(a2 + dx_t*a3)) + fd*(x_int - x_t);
    
    if (exist( 'fd_int', 'var')) then % original, IF (PRESENT(fd_int)) THEN
       fd_int = fd;
    end
    
    if (exist( 'fd2_int', 'var' )) then % original, IF (PRESENT(fd2_int)) THEN
       fd2_int  = 2 * a2 + dx_t * 6 * a3;
    end
end