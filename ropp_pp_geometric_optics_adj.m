function [impact, bangle, impact_dd, impact_dr, bangle_dd, bangle_dr] = ...
    ropp_pp_geometric_optics_adj ...
    ( r_leo, v_leo, r_gns, v_gns, doppler)
%****s* GeometricOptics/ropp_pp_geometric_optics_adj *
%
% NAME
%    ropp_pp_geometric_optics_adj - Calculate bending angle and impact
%                                   parameter from relative Doppler frequency
%                                   shift. ADJOINT VERSION.
%
% SYNOPSIS
%    call ropp_pp_geometric_optics_adj(r_leo, v_leo, r_gns, v_gns, doppler,
%                                      impact, bangle, impact_dd, impact_dr,
%                                      bangle_dd, bangle_dr)
%
% DESCRIPTION
%    This routine is the adjoint of ropp_pp_geometric_optics, which calculates
%    bending angle and impact parameter from relative Doppler frequency shift.
%    Iterative solution of the system of equations:
%               (c - (v_leo, U_leo))/(c - (v_gns, U_gns)) - 1 = doppler
%               [r_leo, U_leo] - [r_gns, U_gns] = 0
%               (U_leo, U_leo) = 1
%               (U_gns, U_gns) = 1
%   where U_leo and U_gns are the ray directions at the receiver and
%   transmitter respectively.
%
% INPUTS
%    real(wp), dimension(:) :: r_leo     % relative LEO position (ECI)
%    real(wp), dimension(:) :: v_leo     % LEO velocity (ECI)
%    real(wp), dimension(:) :: r_gns     % relative GPS position (ECI)
%    real(wp), dimension(:) :: v_gns     % GPS velocity (ECI)
%    real(wp)               :: doppler   % relative Doppler frequency shift
%
% OUTPUT
%    real(wp)               :: impact    % impact parameter (m)
%    real(wp)               :: bangle    % bending angle (rad)
%    real(wp)               :: impact_dd % d(IP)/d(d)
%    real(wp), dimension(:) :: impact_dr % d(IP)/d(r_leo,r_gns)
%    real(wp)               :: bangle_dd % d(BA)/d(d)
%    real(wp), dimension(:) :: bangle_dr % d(BA)/d(r_leo,r_gns)
%
% NOTES
%
% REFERENCES
%   Vorob'ev and Krasil'nikova 1994,
%   Estimation of the accuracy of the atmospheric refractive index recovery
%   from doppler Sshift measurements at frequencies used in the NAVSTAR system
%   Physics of the Atmosphere and Ocean (29) 602-609
%
%  Gorbunov, M.E. and Kornblueh, L. 2003,
%  Principles of variational assimilation of GNSS radio occultation data
%  Max PlancK Institute Report 350
%  http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/max_scirep_350.pdf
%
% AUTHOR
%   M Gorbunov, Russian Academy of Sciences, Russia.
%   Any comments on this software should be given via the ROM SAF
%   Helpdesk at http://www.romsaf.org
%
% COPYRIGHT
%   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
%   For further details please refer to the file COPYRIGHT
%   which you should have received as part of this distribution.
%
%****
% %-------------------------------------------------------------------------------
% % 1. Declarations
% %-------------------------------------------------------------------------------
%
%   USE typesizes, ONLY: wp => EightByteReal
%   USE ropp_pp, not_this => ropp_pp_geometric_optics_adj
%   USE ropp_utils
%
%   IMPLICIT NONE
%
%   REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo      % LEO position (m) [ECI]
%   REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo      % LEO velocity (m/s) [ECI]
%   REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns      % GPS position (m) [ECI]
%   REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns      % LEO velocity (m/s) [ECI]
%   REAL(wp),               INTENT(in)  :: doppler    % Relative Doppler shift
%   REAL(wp),               INTENT(out) :: impact     % Impact parameter (m)
%   REAL(wp),               INTENT(out) :: bangle     % Bending angle (rad)
%   REAL(wp),               INTENT(out) :: impact_dd  % d(IP)/d(d)
%   REAL(wp), DIMENSION(6), INTENT(out) :: impact_dr  % d(IP)/d(r_leo,r_gns)
%   REAL(wp),               INTENT(out) :: bangle_dd  % d(BA)/d(d)
%   REAL(wp), DIMENSION(6), INTENT(out) :: bangle_dr  % d(BA)/d(r_leo,r_gns)
%
%   REAL(wp), DIMENSION(:), ALLOCATABLE :: U_leo   % Ray direction at LEO
%   REAL(wp), DIMENSION(:), ALLOCATABLE :: U_gns   % Ray direction at GPS
%   REAL(wp), DIMENSION(:), ALLOCATABLE :: TC      % Temporary storage
%   REAL(wp), DIMENSION(6)              :: DZ      % Perturbation of (U_leo,U_gns)
%   REAL(wp), DIMENSION(6)              :: DF      % Discrepancy vector
%   REAL(wp), DIMENSION(6,6)            :: A       % Matrix system A*DZ = DF
%   REAL(wp), DIMENSION(6,6)            :: AI      % Inverted system matrix
%   REAL(wp), DIMENSION(6)              :: E_URT   % d(BA)/d(U_leo,U_gns)
%   REAL(wp), DIMENSION(3)              :: P_UR    % d(IP)/d(U_leo)
%   REAL(wp), DIMENSION(3)              :: P_XR    % d(IP)/d(r_leo)
%   REAL(wp), DIMENSION(6)              :: URT_D   % d(U_leo)/d(D)
%   REAL(wp), DIMENSION(3,3)            :: BXR     % Influence of X_leo
%   REAL(wp), DIMENSION(3,3)            :: BXT     % Influence of X_gns
%   REAL(wp), DIMENSION(6,6)            :: B       % Influence of X_leo and X_gns
%   REAL(wp), DIMENSION(6,6)            :: URT_XRT % d(U_leo,U_gns)/d(X_leo,X_gns)
%   INTEGER                             :: i, j, k % Counters
%
%   INTEGER,  PARAMETER :: nit = 10                % Number of iterations
%   REAL(wp), PARAMETER :: tensor(3,3,3) =   &     % Antisymmetrical tensor
%                              RESHAPE((/0,  0,  0,  0,  0,  1,  0, -1,  0,   &
%                                        0,  0, -1,  0,  0,  0,  1,  0,  0,   &
%                                        0,  1,  0, -1,  0,  0,  0,  0,  0/), &
%                                        Shape = (/3,3,3/), Order = (/3,2,1/))
%
% %-------------------------------------------------------------------------------
% % 2. Array allocation
% %-------------------------------------------------------------------------------
%
%   ALLOCATE(U_leo(SIZE(r_leo)))
%   ALLOCATE(U_gns(SIZE(r_gns)))
%   ALLOCATE(TC(SIZE(r_gns)))


Nit = 10; %number of iteration
% % tensor(3,3,3) = RESHAPE((0,  0,  0,  0,  0,  1,  0, -1,  0,    0,  0, -1,  0,  0,  0,  1,  0,  0,  ...
%                                         0,  1,  0, -1,  0,  0,  0,  0,  0); ...
%                                         Shape = (3,3,3), Order = (3,2,1));
tensor (1,:,1) = [0,0,0];
tensor (1,:,2) = [0,0,-1];
tensor (1,:,3) = [0,1,0];
tensor (2,:,1) = [0,0,1];
tensor (2,:,2) = [0,0,0];
tensor (2,:,3) = [-1,0,0];
tensor (3,:,1) = [0,-1,0];
tensor (3,:,2) = [1,0,0];
tensor (3,:,3) = [0,0,0];

U_leo = zeros (1,3);
U_gns = zeros (1,3);
A     = zeros (6,6);
lightvel = 299792458;
%-------------------------------------------------------------------------------
% 3. initial approximation -- straight line r_gns --> r_leo
%-------------------------------------------------------------------------------

U_leo(:) = (r_leo (1,:) - r_gns (1,:)) / sqrt ( sum ( (r_leo (:) - r_gns (:)) .^2 ));
U_gns (:) = U_leo (:);

%-------------------------------------------------------------------------------
% 4. iterative solution
%-------------------------------------------------------------------------------
for k = 1:Nit;
    % 4.1 Calculate discrepancy
    
    DF(1)   = doppler - ( dot(v_gns,U_gns) - dot(v_leo,U_leo)) / ...
        (lightvel - dot(v_gns,U_gns));
    DF(2:4) = cross(r_gns, U_gns) - cross(r_leo, U_leo);
    DF(5)   = 1 - dot(U_leo,U_leo);
    DF(6)   = 1 - dot(U_gns,U_gns);
    
    % 4.2 Calculate the matrix of linearized system
    A(:,:)   = 0;
    A(1,1:3) = -v_leo / (lightvel - dot(v_gns,U_gns));
    A(1,4:6) =  v_gns * ((lightvel - dot(v_leo,U_leo)) / ((lightvel - dot(v_gns,U_gns))^2));
    for i=1:3
        for j=1:3
            A(i+1,j)   =  sum(tensor(i,:,j)*r_leo(:)); % is the sign correct?
            A(i+1,j+3) = -sum(tensor(i,:,j)*r_gns(:));
        end
    end
    A(5,1:3) = 2 * U_leo;
    A(6,4:6) = 2 * U_gns;
    
    % 4.3 Solve liearized system
    
%     AI = inv (A);
    
    [ AI ] = ropp_pp_invert_matrix( A ) ;
    DZ = AI * DF'; % ori: DZ = MATMUL(AI, DF)
    DZ = DZ';
    % Due to fortran the matrix is not row - column
    
    % Calculation of next approximation
    
    TC = DZ(1:3);
    U_leo = U_leo + TC;
    TC = DZ(4:6);
    U_gns = U_gns + TC;
end
%-------------------------------------------------------------------------------
% 5. Calculate bending angle and impact parameter
%-------------------------------------------------------------------------------
% 5.1 Bending angle
%     bangle = vector_angle(U_gns, U_leo, vector_product(r_gns, r_leo));
% CALL vector_angle_adj(U_gns, U_leo, vector_product(r_gns, r_leo), bangle,   &
% E_URT(4:6), E_URT(1:3))

[bangle, E_URT(4:6), E_URT(1:3) ] = vector_angle_adj ( U_gns, U_leo, cross(r_gns,r_leo) );


% 5.2 impact = r_L sin(phi_l) = r_G sin(phi_G)

%  impact = SQRT(DOT_PRODUCT(r_leo, r_leo)) * SIN(vector_angle(r_leo, U_leo))

impact = norm( r_leo) * sin ( acos (dot (r_leo, U_leo) / (norm (r_leo)* norm (U_leo))));


%-------------------------------------------------------------------------------
% 6. Adjoint calculations
%-------------------------------------------------------------------------------

% 6.1 Derivatives of U_leo and U_gns

%   URT_D(:) = AI(:,1)

URT_D = AI(:,1);

%   DO i=1,3
%      DO j=1,3
%         BXR(i,j) = -SUM(tensor(i,j,:)*U_leo(:))
%         BXT(i,j) =  SUM(tensor(i,j,:)*U_gns(:))
%      ENDDO
%   ENDDO


for i = 1:3
    for j = 1:3
        tmp1 = tensor (i,j,1);
        tmp2 = tensor (i,j,2);
        tmp3 = tensor (i,j,3);
        tensortmp = [tmp1, tmp2, tmp3];
        %         tensortmp = [tensor (i,j,1),tensor(i,j,2),tensor(i,j,3)];
        BXR(i,j) = -sum ( tensortmp* U_leo(:));
        BXT(i,j) = sum(tensortmp*U_gns(:));
    end
end

%
%   B(:,:)       = 0
%   B(2:4,1:3)   = BXR(:,:)
%   B(2:4,4:6)   = BXT(:,:)
%   URT_XRT(:,:) = MATMUL(AI(:,:),B(:,:))




B(:,:)       = 0;
B(2:4,1:3)   = BXR(:,:);
B(2:4,4:6)   = BXT(:,:);
URT_XRT = AI * B';

% 6.2 Derivatives of bending angle

%
%   bangle_dd = SUM(E_URT(:)*URT_D(:))
%   bangle_dr = MATMUL(E_URT(:),URT_XRT(:,:))

bangle_dd = sum (E_URT*URT_D);
bangle_dr = E_URT * URT_XRT;

% 6.3 Derivatives of impact parameter

%   P_UR(:)    = -(vector_product(r_leo,(vector_product(r_leo,U_leo))))/impact
%   P_XR(:)    = -(vector_product(U_leo,(vector_product(U_leo,r_leo))))/impact
%
%   impact_dd  = SUM(P_UR(:)*URT_D(1:3))
%
%   impact_dr(:)   = MATMUL(P_UR(:),URT_XRT(1:3,:))
%   impact_dr(1:3) = impact_dr(1:3) + P_XR(:)


P_UR   = -(cross(r_leo,(cross(r_leo,U_leo))))/impact;
P_XR   = -(cross(U_leo,(cross(U_leo,r_leo))))/impact;

impact_dd  = sum (P_UR*URT_D(1:3));

impact_dr  = P_UR * URT_XRT(1:3,:);
impact_dr(1:3) = impact_dr(1:3) + P_XR;




    function [angle, AXY_X, AXY_Y ] = vector_angle_adj ( X, Y, A )
        
        % Function verified against ROPP on July 23, 2015
        
        %    SUBROUTINE vector_angle_adj(X, Y, A, angle,AXY_X, AXY_Y)
        
        %     USE typesizes,  ONLY: wp => EightByteReal
        %     USE ropp_utils, ONLY: vector_product
        
        %     REAL(wp), DIMENSION(3), INTENT(in)  :: X          % Cartesian vector
        %     REAL(wp), DIMENSION(3), INTENT(in)  :: Y          % Cartesian vector
        %     REAL(wp), DIMENSION(3), INTENT(in)  :: A          % Orientation axis
        %     REAL(wp),               INTENT(out) :: angle      % Angle between X and Y
        %     REAL(wp), DIMENSION(3), INTENT(out) :: AXY_X      % d(AXY)/d(X)
        %     REAL(wp), DIMENSION(3), INTENT(out) :: AXY_Y      % d(AXY)/d(Y)
        
        %     REAL(wp), DIMENSION(3) :: n, alpha, beta, gamma
        %     REAL(wp)               :: nn
        %     REAL(wp)               :: ag, bg
        
        %     nn = DOT_PRODUCT(A, A)
        nn = dot ( A, A);
        
        %     IF (nn == 0) THEN
        %        angle = 0.0_wp
        %        axy_x = 0.0_wp
        %        axy_y = 0.0_wp
        %     ELSE
        %        n = A/SQRT(nn)
        %        alpha = vector_product(n, X)
        %
        %        beta = X - DOT_PRODUCT(n, X) * n
        %        gamma = Y - DOT_PRODUCT(n, Y) * n
        %
        %        ag = DOT_PRODUCT(alpha,gamma)
        %        bg = DOT_PRODUCT(beta,gamma)
        %
        %        angle = ATAN2(ag, bg)
        %
        %        AXY_X = -(ag*gamma - bg*(vector_product(Y,n)))/(ag**2 + bg**2)
        %        AXY_Y = -(ag*beta  - bg*alpha) / (ag**2 + bg**2)
        %
        %     ENDIF
        
        if ( nn == 0)
            angle = 0;
            AXY_X = 0;
            AXY_Y = 0;
        else
            n = A/sqrt ( nn);
            alpha =  cross(n, X);
            beta = X - dot ( n,X ) * n;
            gamma = Y - dot ( n,Y) * n;
            ag = dot ( alpha, gamma);
            bg = dot ( beta , gamma);
            angle = atan2 ( ag, bg);
            
            AXY_X = -(ag*gamma - bg*( cross(Y,n) ))/(ag^2 + bg^2);
            AXY_Y = -(ag*beta  - bg*alpha) / (ag^2 + bg^2);
        end
        %   END SUBROUTINE vector_angle_adj
        
    end

end