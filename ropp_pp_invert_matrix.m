function [ B ] = ropp_pp_invert_matrix( A )
%****s* PPUtils/ropp_pp_invert_matrix *
%
% NAME
%    ropp_pp_invert_matrix - Invert a matrix A(dimY, dimX)
%
% SYNOPSIS
%    call ropp_pp_invert_matrix(A, B)
%
% DESCRIPTION
%    Gauss elimination
%
%****
%
% SUBROUTINE ropp_pp_invert_matrix(A, B)

%   IMPLICIT NONE

% 2.1 Declarations

%   REAL(wp), DIMENSION(:,:), INTENT(inout) :: A     % Matrix to invert
%   REAL(wp), DIMENSION(:,:), INTENT(out)   :: B     % Inverted matrix

%   INTEGER  :: N             % Matrix dimension
%   INTEGER  :: i, k          % Matrix indeces
%   INTEGER  :: m(1)          % Search index
%   REAL(wp) :: Alpha         % Row combination coefficien
%   REAL(wp) :: Det           % Matrix determinant
%   REAL(wp) :: Amin, Amax    % Max and min diagonal elements
%   REAL(wp), DIMENSION(:), ALLOCATABLE :: F  % Exchange work space

% 2.2 Check array shape

%   IF (SIZE(A,1) /= SIZE(A,2)) THEN
%      RETURN
%   ENDIF
%   IF (ANY(SHAPE(B) /= SHAPE(B))) THEN
%      RETURN
%   ENDIF

if size ( A,1 ) ~= size ( A,2 )
    return
end
% Dont know how to translate this

%   ALLOCATE(F(SIZE(A,1)))

F = zeros ( size ( A,1) );

% 2.3 Initialization

%   N      = SIZE(A,1)
%   B(:,:) = 0.0_wp
%
%   DO i=1,N
%      B(i,i) = 1.0_wp
%   ENDDO


N = size ( A, 1);
% B = 0;

for i = 1:N
    B(i,i) = 1;
end


% 2.4 Elimination of sub-diagonal elements

%   Lower: DO k=1,N-1
%
%      IF (A(k,k) == 0.0_wp) THEN
%         m(:) = MAXLOC(A(k+1:N,k), A(k+1:N,k) /= 0.0_wp)
%         IF (m(1) == 0) THEN
%            EXIT Lower
%         END IF
%         i      = k + m(1)
%         F(:)   = A(k,:)
%         A(k,:) = A(i,:)
%         A(i,:) = F(:)
%         F(:)   = B(k,:)
%         B(k,:) = B(i,:)
%         B(i,:) = F(:)
%      END IF
%
%      DO i=k+1,N
%         Alpha  = A(i,k)/A(k,k)
%         A(i,:) = A(i,:) - Alpha*A(k,:)
%         B(i,:) = B(i,:) - Alpha*B(k,:)
%      END DO
%
%   END DO Lower


% Lower
for k = 1:N-1
    
    if ( A (k,k) == 0 )
        %         m (:) = find ( A(k+1:N, k) ~= 0 )
        
        newmat = A(k+1:N,k);
        newidx = find ( newmat ~= 0 ) ;
        newmat2 = newmat ( newidx, k);
        m  = max ( newmat2 );
        if isempty (m)
            m =0;
        end
        
        
        if ( m == 0 )
            break;
        end
        i = k + m(1);
        F(:) = A (k,:);
        A(k,:) = A (i,:);
        A(i,:) = F(:);
        F(:)   = B(k,:);
        B(k,:) = B(i,:);
        B(i,:) = F(:);
        
    end
    
    for i = k+1:N
        Alpha = A(i,k) / A (k,k);
        A(i,:) = A(i,:) - Alpha*A(k,:);
        B(i,:) = B(i,:) - Alpha*B(k,:);
    end
end



% 2.5 Checking for degenerated matrix
%
%   Det  = 1.0_wp
%   Amin = ABS (A(1,1))
%   Amax = ABS (A(1,1))
%   Diagonal: DO i=1,N
%      Amin = MIN(ABS (A(i,i)), Amin)
%      Amax = MAX(ABS (A(i,i)), Amax)
%      Det  = Det*A(i,i)
%   END DO Diagonal
%
%   IF (Det == 0.0_wp) THEN
%      RETURN
%   END IF


Det = 1;
Amin = abs ( A(1,1) );
Amax = abs ( A(1,1) );
for i = 1:N
    Amin = min ( abs ( A(i,i)), Amin );
    Amax = max ( abs ( A(i,i)), Amax );
    Det = Det * A ( i,i );
end

if ( Det == 0 )
    return
end


% 2.6 Elimination of super-diagonal elements

%   Upper: DO k=N,2,-1
%      DO i=k-1,1,-1
%         Alpha  = A(i,k)/A(k,k)
%         A(i,:) = A(i,:) - Alpha*A(k,:)
%         B(i,:) = B(i,:) - Alpha*B(k,:)
%      END DO
%   END DO Upper
%
%   DO i=1,N
%      B(i,:) = B(i,:)/A(i,i)
%   END DO

for k = N:-1:2
    for i = k-1:-1:1
        Alpha = A(i,k) / A(k,k);
        A(i,:) = A(i,:) - Alpha*A(k,:);
        B(i,:) = B(i,:) - Alpha*B(k,:);
    end
end

for i = 1:N
    B (i,:) = B(i,:)/A(i,i);
end




%   DEALLOCATE(F)
%
% END SUBROUTINE ropp_pp_invert_matrix

end

