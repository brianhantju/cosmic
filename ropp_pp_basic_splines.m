function [ K ] = ropp_pp_basic_splines( X, Xs )
%ROPP_PP_BASIC_SPLINES Summary of this function goes here
% !****s* PPSpline/ropp_pp_basic_splines *
% !
% ! NAME
% !    ropp_pp_basic_splines - Generate a matrix of basic polynomials for
% !                            polynomial regression
% !
% ! SYNOPSIS
% !    call ropp_pp_basic_splines(X, Xs, K)
% !
% ! DESCRIPTION
% !    Generation of matrix of basic polynomials for polynomial regression
% !
% !****
% !
% ! 1.1 Declarations

%     REAL(wp), DIMENSION(:), INTENT(in)    :: X  ! Grid of argument Z
%     REAL(wp), DIMENSION(:), INTENT(in)    :: Xs ! X-grid of delta-splines
%     REAL(wp), DIMENSION(:,:), INTENT(out) :: K  ! Matrix of basic polynomials

%     INTEGER                       :: i    ! x index
%     INTEGER                       :: j    ! Delta-spline number
%     REAL(wp), DIMENSION(:), ALLOCATABLE :: S    ! Delta-spline
%     REAL(wp), DIMENSION(:), ALLOCATABLE :: D2S  ! Delta-spline 2nd derivative

%     ALLOCATE(S(SIZE(xs)))
%     ALLOCATE(D2S(SIZE(xs)))

S = zeros ( 1,length(Xs));
D2S = zeros (1,length (Xs));

% ! 1.2 Generate basic function matrix K

S(:) = 0;

for j=1:size ( K ,2)   %ori: j=1, SIZE(K,2) %% why us K before define it?
    S(j) = 1;
    % CALL ropp_pp_init_spline(Xs, S, D2S)
    
    [ D2S ] = ropp_pp_init_spline( Xs, S );
    
    for i=1:length(X)
        % CALL ropp_pp_interpol_spline(Xs, S, D2S, X(i), K(i,j))
        
        [ K(i,j) ] = ropp_pp_interpol_spline( Xs, S, D2S, X(i) );
    end
    S(j) = 0;
end  

end