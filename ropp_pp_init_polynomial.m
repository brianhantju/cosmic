function [ K ] = ropp_pp_init_polynomial( x, K )
% !****s* PPUtils/ropp_pp_init_polynomial *
% !
% ! NAME
% !    ropp_pp_init_polynomial - Generate matrix of basic polynomials
% !                              for regression
% !
% ! SYNOPSIS
% !    call ropp_pp_init_polynomial(x, K)
% !
% ! DESCRIPTION
% !    Compute basic polynomials as
% !        K_ij = f_j(x_i) = (x_i)^j     for j=0..UBound(K,2)
% !
% !****
% !
%   SUBROUTINE ropp_pp_init_polynomial(x, K)
%
%     USE typesizes, ONLY: wp => EightByteReal
%
%     ! 5.1 Declarations
%
%     IMPLICIT NONE

%     REAL(wp), DIMENSION(1:),    INTENT(in)  :: x  ! grid of argument x
%     REAL(wp), DIMENSION(1:,0:), INTENT(out) :: K  ! matrix of polynomials
%     INTEGER                                 :: i,j

%     ! 5.2 Generate polynomials


K (:,1) = 1;
for j = 2:size(K,2) 
    for i = 1:length(x)
        if (x(i) > 0)
            %           K(i, j) = x(i)**REAL(j,wp)
            K(i, j) = x(i)^ (j-1);
        else
            %           K(i, j) = x(i)**INT(REAL(j,wp))
%             K(i, j) = x(i)^int64(real(j-1)); % ?is this correct?
            
            % Added on Jul 22, 2015
            K(i, j) = x(i)^ (j-1);
        end
    end
end%   END SUBROUTINE ropp_pp_init_polynomial


end