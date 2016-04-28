function [ P, DP ] = ropp_pp_polynomial( c, x )
if 0 
% ****s* PPUtils/ropp_pp_polynomial *
% 
%  NAME
%     ropp_pp_polynomial - Calculation of a polynomial and its derivative
% 
%  SYNOPSIS
%     call ropp_pp_polynomial(c, x, P, DP)
% 
%  DESCRIPTION
%     Horner scheme
% 
% ****
% 
%   SUBROUTINE ropp_pp_polynomial(c, x, P, DP)
%
%     USE typesizes, ONLY: wp => EightByteReal
%
%      4.1 Declarations
%
%     IMPLICIT NONE
%
%     REAL(wp), DIMENSION(0:), INTENT(in)  :: c      polynomial coefficients
%     REAL(wp),                INTENT(in)  :: x      polynomial argument
%     REAL(wp),                INTENT(out) :: P      polynomial value
%     REAL(wp), OPTIONAL,      INTENT(out) :: DP     polynomial derivative
%
%     INTEGER :: n, i

%  ! 4.2 Calculate polynomial value
% 
%     n = UBOUND(c,1)
%     
%     P = c(n)
%     DO i = n-1, 0, -1
%        P = c(i) + x*P
%     ENDDO
%     
%     ! 4.3 Calculate polynomial derivative
%     
%     IF ( PRESENT(DP) ) THEN
%        DP = c(n)*REAL(n,wp)
%        DO i = n-1, 1, -1
%           DP = c(i)*REAL(i,wp) + x*DP
%        ENDDO
%     ENDIF
end


% ROPP_PP_POLYNOMIAL This function evaluate the polynomial value given coefficients
% Input
%   c : polynomial coefficients, 1*n , higher order to lower order, eg c =
%       [12,3,4] means polynomial of 12x^2 + 3x^1 + 4
%   x : polynomial argument, single value real number
% Output
%   P : polynomial value
%   DP: polynomial derivative
% Dependency Functions
%  
% Called in
%   C:\Users\bhan002\Dropbox\matlab\00 ROPP functions\ropp_pp_satellite_velocities.m
%
% Created by    Creation Date
%   HAN BO      Jul. 20, 2015
%
% Modified by  Modification Date
%   HAN BO      Jul. 20, 2015
%
% Reference
%   ROPP_PP module, the function inside the package 
% Copyright 2012-2016 HAN BO
% $ Revision: 1.0 $

c = c(end:-1:1);

% 4.2 Calculate polynomial value 
n = length(c);
P = c(end);
for i = (n-1):-1:1
    P = c(i) + x*P;
end

% 4.3 Calculate polynomial derivative
if nargout > 1
    DP = c(end)* (n-1);
    for i = n-1:-1:2
        DP = c(i)* (i-1) + x*DP;
    end
end
end
