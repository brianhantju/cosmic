function [ C ] = ropp_pp_matmul_sub( A, B )
% SUBROUTINE ropp_pp_matmul_sub(A, B, C)
%
%   IMPLICIT NONE
%
%   REAL(wp), DIMENSION(:,:), INTENT(in)  :: A
%   REAL(wp), DIMENSION(:,:), INTENT(in)  :: B
%   REAL(wp), DIMENSION(:,:), INTENT(out) :: C
%   INTEGER :: i, j, k
% DO i=1,SIZE(A,1)
%     DO j=1, SIZE(B,2)
%       C(i,j) = 0.0_wp
%       DO k=1, SIZE(B,1)
%         C(i,j) = C(i,j) + A(i,k)*B(k,j)
%       ENDDO
%     ENDDO
%   ENDDO
%
for i=1:size(A,1)
    for j=1: size(B,2)
        C(i,j) = 0;
        for k=1: size(B,1)
            C(i,j) = C(i,j) + A(i,k) * B(k,j);
        end
    end
end
% END SUBROUTINE ropp_pp_matmul_sub
end

