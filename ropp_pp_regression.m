function [ a ] = ropp_pp_regression( K, y )
% !****s* PPUtils/ropp_pp_regression *
% !
% ! NAME
% !    ropp_pp_regression - Linear regression
% !
% ! SYNOPSIS
% !    call ropp_pp_regression(K, y, a)
% !
% ! DESCRIPTION
% !    Quasi-inversion of matrix of basic functions:
% !       || Sum_j a_j f_j(x_i) - y_i || = min
% !    Solution is a = Q y where Q is left inverse of K_ij = f_j(x_i) 
% !
% !****
%   SUBROUTINE ropp_pp_regression(K, y, a)

%     USE typesizes, ONLY: wp => EightByteReal
    
%     ! 6.1 Declarations
    
%     IMPLICIT NONE
    
%     REAL(wp), DIMENSION(:,:), INTENT(in)  :: K  ! matrix of functions
%     REAL(wp), DIMENSION(:),   INTENT(in)  :: y  ! regression data  (n*1)
%     REAL(wp), DIMENSION(:),   INTENT(out) :: a  ! regression coefficients
%     : The first number is the lower order of the polynomial
%     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Q ! quasi-inverse matrix
%     
%     ALLOCATE(Q(SIZE(K,2), SIZE(K,1)))


% The second input to ropp_pp_regression should be col vector
% The output of ropp_pp_regression is a column vector, with the first
% element being the lower order of the polynomial

% Convert y to column vector if it is not
if size (y,1) == 1 && size (y,2) > 1
    y = y';
end

    

%     ! 6.2 Invert matrix of functions
        
%     CALL ropp_pp_quasi_invert(K, Q)
    Q = ropp_pp_quasi_invert (K);

%     ! 6.3 Solve to find coefficients
    
    
%     a(:) = ropp_pp_matmul(Q, y(:));
      a = Q * y;
      % Dont know if it is correct

%     DEALLOCATE(Q)
    
%   END SUBROUTINE ropp_pp_regression

end

