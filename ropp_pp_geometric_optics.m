function [impact, bangle] = ropp_pp_geometric_optics ( r_leo, v_leo, r_gns, v_gns, doppler)
%!****s* GeometricOptics/ropp_pp_geometric_optics *
%!
%! NAME
%!    ropp_pp_geometric_optics - Calculate bending angle and impact parameter 
%!                               from relative Doppler frequency shift.
%!                   
%! SYNOPSIS
%!    call ropp_pp_geometric_optics(r_leo, v_leo, r_gns, v_gns, doppler,
%!                                  impact, bangle)
%! 
%! DESCRIPTION
%!    This routine calculates bending angle and impact parameter from relative
%!    Doppler frequency shift.
%!    Iterative solution of the system of equations:
%!               (c - (v_leo, U_leo))/(c - (v_gns, U_gns)) - 1 = doppler
%!               [r_leo, U_leo] - [r_gns, U_gns] = 0
%!               (U_leo, U_leo) = 1
%!               (U_gns, U_gns) = 1
%!   where U_leo and U_gns are the ray directions at the receiver and 
%!   transmitter respectively.
%!
%! INPUTS
%!    real(wp), dimension(:) :: r_leo     ! relative LEO position (m) [ECI]
%!    real(wp), dimension(:) :: v_leo     ! LEO velocity (m/s) [ECI]  
%!    real(wp), dimension(:) :: r_gns     ! relative GPS position (m) [ECI]
%!    real(wp), dimension(:) :: v_gns     ! GPS velocity (m/s) [ECI]  
%!    real(wp)               :: doppler   ! relative Doppler frequency shift
%!
%! OUTPUT
%!    real(wp)               :: impact    ! impact parameter (m)
%!    real(wp)               :: bangle    ! bending angle (rad)


% REAL(wp), DIMENSION(:), INTENT(in)  :: r_leo   ! LEO position (m) [ECI]
% REAL(wp), DIMENSION(:), INTENT(in)  :: v_leo   ! LEO velocity (m/s) [ECI]
% REAL(wp), DIMENSION(:), INTENT(in)  :: r_gns   ! GPS position (m) [ECI]
% REAL(wp), DIMENSION(:), INTENT(in)  :: v_gns   ! GPS velocity (m/s) [ECI]
% REAL(wp),               INTENT(in)  :: doppler ! Relative Doppler shift
% REAL(wp),               INTENT(out) :: impact  ! Impact parameter (m)
% REAL(wp),               INTENT(out) :: bangle  ! Bending angle (rad)
% REAL(wp), DIMENSION(:), ALLOCATABLE :: U_leo   ! Ray direction at receiver
% REAL(wp), DIMENSION(:), ALLOCATABLE :: U_gns   ! Ray direction at transmitter
% REAL(wp), DIMENSION(:), ALLOCATABLE :: TC      ! Temporary storage
% REAL(wp), DIMENSION(6)              :: DZ      ! Perturbation of (U_leo,U_gns)
% REAL(wp), DIMENSION(6)              :: DF      ! Discrepancy vector
% REAL(wp), DIMENSION(6,6)            :: A       ! Matrix system A*DZ = DF
% REAL(wp), DIMENSION(6,6)            :: AI      ! Inverted system matrix
% INTEGER                             :: i, j, k ! Counters
% 
% INTEGER,  PARAMETER :: Nit = 5                 ! Number of iterations
% REAL(wp), PARAMETER ::       &                 ! Antisymmetrical tensor
% tensor(3,3,3) = RESHAPE((/0,  0,  0,  0,  0,  1,  0, -1,  0,   &
% 0,  0, -1,  0,  0,  0,  1,  0,  0,   &
% 0,  1,  0, -1,  0,  0,  0,  0,  0/), &
% Shape = (/3,3,3/), Order = (/3,2,1/))


Nit = 5; %number of iteration
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
A = zeros (6,6);
lightvel = 299792458;
%-------------------------------------------------------------------------------
% initial approximation -- straight line r_gns --> r_leo
%-------------------------------------------------------------------------------

U_leo(:) = (r_leo (1,:) - r_gns (1,:)) / sqrt ( sum ( (r_leo (:) - r_gns (:)) .^2 ));
U_gns (:) = U_leo (:);

%-------------------------------------------------------------------------------
% iterative solution
%-------------------------------------------------------------------------------
for k = 1:Nit;
	%Calculate discrepancy
    
     DF(1)   = doppler - ( dot(v_gns,U_gns) - dot(v_leo,U_leo)) / ...
                             (lightvel - dot(v_gns,U_gns));
     DF(2:4) = cross(r_gns, U_gns) - cross(r_leo, U_leo);
     DF(5)   = 1 - dot(U_leo,U_leo); % is it 1?
     DF(6)   = 1 - dot(U_gns,U_gns);
	 
	 % Calculate the matrix of linearized system
	 A(:,:)   = 0;
     A(1,1:3) = -v_leo / (lightvel - dot(v_gns,U_gns));
     A(1,4:6) =  v_gns * (lightvel - dot(v_leo,U_leo)) /   ...
                             (lightvel - dot(v_gns,U_gns))^2;
     for i=1:3
        for j=1:3
           A(i+1,j)   =  sum(tensor(i,:,j)*r_leo(:)); % is the sign correct?
           A(i+1,j+3) = -sum(tensor(i,:,j)*r_gns(:));
        end
     end
     A(5,1:3) = 2 * U_leo;
     A(6,4:6) = 2 * U_gns; % is it 2?
	 
	 % Solve liearized system

     AI = inv (A);
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
% Calculate bending angle and impact parameter
%-------------------------------------------------------------------------------
    % Bending angle 
%     bangle = vector_angle(U_gns, U_leo, vector_product(r_gns, r_leo));
    
    n = cross (r_gns, r_leo) / sqrt ( dot ( cross(r_gns, r_leo) , cross(r_gns, r_leo) ));
    tmpa = cross (n, U_gns);
    tmpb = U_gns - dot (n,U_gns) * n;
    tmpc = U_leo - dot (n, U_leo) * n;
    bangle = atan2 ( dot(tmpa, tmpc) , dot (tmpb, tmpc ) );    
    
    % impact = r_L sin(phi_l) = r_G sin(phi_G)
    impact = norm( r_leo) * sin ( acos (dot (r_leo, U_leo) / (norm (r_leo)* norm (U_leo))));
end 