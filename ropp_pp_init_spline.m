function [ d2 ] = ropp_pp_init_spline( x, f )
%ROPP_PP_INIT_SPLINE
% !****s* PPSpline/ropp_pp_init_spline *
% !
% ! NAME
% !    ropp_pp_init_spline - Calculation of 2nd derivative of a spline
% !
% ! SYNOPSIS
% !    call ropp_pp_init_spline(x, f, d2)
% !
% ! DESCRIPTION
% !    Calculate 2nd derivative of a spline, drive-through solution
% !
% !****

%     ! 2.1 Declarations
    
%     REAL(wp), DIMENSION(:), INTENT(in)  :: x  ! Argument grid (monotonous)
%     REAL(wp), DIMENSION(:), INTENT(in)  :: f  ! Gridded function
%     REAL(wp), DIMENSION(:), INTENT(out) :: d2 ! 2nd derivative of spline

%     REAL(wp), DIMENSION(:), ALLOCATABLE :: d1
%     REAL(wp)                     :: dfl, dfr, df, a, b, c
%     INTEGER                      :: i, N

%     ! 2.2 Initialisation

    N =  length (x); % ALLOCATE(d1(N))
    d1 = zeros ( 1:N);
    d2(:) = 0;
    
%     ! 2.3 Drive-through calculation of spline coefficients
    
    d1(1) = 0;
    d2(N) = 0;
    d1(1) = 0;

    for i = 2: N-1

       dfl   = (f(i) - f(i-1))/(x(i) - x(i-1));
       dfr   = (f(i+1) - f(i))/(x(i+1) - x(i));
       df    = (dfr - dfl)/(x(i+1) - x(i-1));
       a     = (x(i) - x(i-1))/(2*(x(i+1) - x(i-1)));
       b     = (x(i+1) - x(i))/(2*(x(i+1) - x(i-1)));
       c     = 1 + a*d1(i-1);
       d1(i) = -b/c;
       d2(i) = (3*df - a*d2(i-1))/c;
    end
    
    for i = N-1:-1: 2 % original fortran is:  DO i = N-1, 2, -1
       d2(i) = d1(i)*d2(i+1) + d2(i);
    end
end