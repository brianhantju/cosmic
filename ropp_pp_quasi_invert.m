function [ Q ] = ropp_pp_quasi_invert( K )
% !****s* PPUtils/ropp_pp_quasi_invert *
% !
% ! NAME
% !    ropp_pp_quasi_invert - Quasi-inverse of a matrix K(dimY, dimX)
% !
% ! SYNOPSIS
% !    call ropp_pp_quasi_invert(K, Q)
% !
% ! DESCRIPTION
% !    1. dimY >= dimX:
% !       x = Qy - vector minimizing ||Kx - y||
% !       QK = E; Q is left inverse operator.
% !       Q = (K^T K)^-1 K^T
% !    2. dimY <= dimX:
% !       x = Qy - solution minimizing ||x||
% !       KQ = E; Q is right inverse operator.
% !       Q = K^T (KK^T)^-1
% !
% !****
% !
%   SUBROUTINE ropp_pp_quasi_invert(K, Q)

%     IMPLICIT NONE

%     ! 3.1 Declarations

%     REAL(wp), DIMENSION(:,:), INTENT(in)  :: K  ! Matrix to quasi-invert
%     REAL(wp), DIMENSION(:,:), INTENT(out) :: Q  ! Quasi-inverse
%
%     INTEGER :: DimX, DimY  ! Dimension of x and y spaces
%     INTEGER                               :: N  ! Work matrix dimension
%     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: W  ! Work arrays for inversion
%     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: WI ! Work arrays for inversion
%     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KT

%     ! 3.2 Initialization

DimX = size(K,2);
DimY = size(K,1);

%     if (ANY    (shape(Q) /= (/ DimX, DimY /))  )
%        return
%     end
%     Enter the loop if the dim dont agree with each other

if 0
    if ( size(Q,1) ~= DimX) || (size(Q,2) ~= DimY ) % why use Q before define
        return
    end
end


N = min ( DimX, DimY);
%     ALLOCATE (W(N,N))
%     ALLOCATE (WI(N,N))
%     ALLOCATE (KT(DimX,DimY))

KT = K';

%     ! 3.3 Quasi-Inversion

if (DimY >= DimX)
    
    %        ! 3.3.1 Left inversion   W=K^TK, WI=(K^TK)^-1, Q=(K^TK)^-1 K^T
    
    %        W = ropp_pp_matmul(KT, K)
    W = KT * K;
    %        CALL ropp_pp_invert_matrix(W, WI)
    WI = inv ( W );
    % !       Q = ropp_pp_matmul(WI, KT)
    %        CALL ropp_pp_matmul_sub(WI, KT, Q)
    Q = ropp_pp_matmul_sub ( WI, KT);
    
else
    
    %        ! 3.3.2. Right inversion W=KK^T, WI=(KK^T)^-1, Q=K^T(KK^T)^-1
    
    %        W = ropp_pp_matmul(K, KT)
    W = K * KT;
    %        CALL ropp_pp_invert_matrix(W, WI)
    WI = inv ( W);
    % !       Q = ropp_pp_matmul(KT, WI)
    %        CALL ropp_pp_matmul_sub(KT, WI, Q)
    Q = ropp_pp_matmul_sub ( KT, WI);
    
end

%     ! 3.4 Clean up
%
%     DEALLOCATE (W)
%     DEALLOCATE (WI)
%     DEALLOCATE (KT)

%   END SUBROUTINE ropp_pp_quasi_invert

end

