function [ xleo, vleo, xgns, vgns, abl, abg] = ...
    ropp_pp_satellite_velocities( time, r_leo, r_gns )
if 0
    %     ! $Id: ropp_pp_satellite_velocities.f90 2021 2009-01-16 10:49:04Z frhl $
    %
    % SUBROUTINE ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,   &
    %                                         xgns, vgns, abl, abg)
    %
    % !****s* Preprocessing/ropp_pp_satellite_velocities *
    % !
    % ! NAME
    % !    ropp_pp_satellite_velocities - Calculate satellite velocities by
    % !                                   polynomial regression
    % !
    % ! SYNOPSIS
    % !    call ropp_pp_satellite_velocities(time, r_leo, r_gns, xleo, vleo,    &
    % !                                      xgns, vgns, abl, abg)
    % !
    % ! DESCRIPTION
    % !    This routine calculates satellite velocities by polynomial regression
    % !
    % ! INPUTS
    % !    real(wp), dimension(:)   :: time    ! time of sample (s)
    % !    real(wp), dimension(:,:) :: r_leo   ! cartesian LEO position (ECI or ECF)
    % !    real(wp), dimension(:,:) :: r_gns   ! cartesian GPS position (ECI or ECF)
    % !
    % ! OUTPUTS
    % !    real(wp), dimension(:,:) :: xleo    ! LEO positions from regression
    % !    real(wp), dimension(:,:) :: vleo    ! LEO velocities from regression
    % !    real(wp), dimension(:,:) :: xgns    ! GPS positions from regression
    % !    real(wp), dimension(:,:) :: vgns    ! GPS velocities from regression
    % !    real(wp), dimension(:,:) :: abl     ! LEO regression coefficients (opt)
    % !    real(wp), dimension(:,:) :: abg     ! GPS regression coefficients (opt)
    % !
    % ! NOTES
    % !
    % ! REFERENCES
    % !
    % ! AUTHOR
    % !   M Gorbunov, Russian Academy of Sciences, Russia.
    % !   Any comments on this software should be given via the ROM SAF
    % !   Helpdesk at http://www.romsaf.org
    % !
    % ! COPYRIGHT
    % !   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
    % !   For further details please refer to the file COPYRIGHT
    % !   which you should have received as part of this distribution.
    % !
    % !****
    %
    % !-------------------------------------------------------------------------------
    % ! 1. Declarations
    % !-------------------------------------------------------------------------------
    %
    %   USE typesizes, ONLY: wp => EightByteReal
    %   USE ropp_pp, not_this => ropp_pp_satellite_velocities
    %
    %   IMPLICIT NONE
    %
    %   REAL(wp), DIMENSION(:), INTENT(in)    :: time    ! Time of samples (s)
    %   REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo   ! Cartesian LEO position (m)
    %   REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns   ! Cartesian GPS position (m)
    %   REAL(wp), DIMENSION(:,:), INTENT(out) :: xleo    ! LEO position vector (m)
    %   REAL(wp), DIMENSION(:,:), INTENT(out) :: vleo    ! LEO velocity vector (m)
    %   REAL(wp), DIMENSION(:,:), INTENT(out) :: xgns    ! GPS position vector (m)
    %   REAL(wp), DIMENSION(:,:), INTENT(out) :: vgns    ! GPS velocity vector (m)
    %   REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abl ! LEO regression coeffs
    %   REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(out) :: abg ! GPS regression coeffs
    %
    %   INTEGER, PARAMETER :: nv = 5  ! polynomial degree for calculation of velocity
    %   INTEGER :: n, i, m
    %   REAL(wp), DIMENSION(:),   ALLOCATABLE :: t_norm      ! Normalised time
    %   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: KV          ! Regress matrix velocity
    %   REAL(wp), DIMENSION(0:nv,3)           :: coeff_rleo  ! Coefficients for RLEO
    %   REAL(wp), DIMENSION(0:nv,3)           :: coeff_rgns  ! Coefficients for RGNS
    %
    % !-------------------------------------------------------------------------------
    % ! 2. Initialisation
    % !-------------------------------------------------------------------------------
    %
    %   n = SIZE(time)
    %
    %   ALLOCATE(t_norm(n))
    %   ALLOCATE(KV(n, 0:nv))
    %
    % !-------------------------------------------------------------------------------
    % ! 3. Re-normalise time
    % !-------------------------------------------------------------------------------
    %
    %   t_norm(:) = (time(:) - time(1))/(time(n) - time(1))
    %
    % !-------------------------------------------------------------------------------
    % ! 4. Calculate regression coefficients
    % !-------------------------------------------------------------------------------
    %
    %   CALL ropp_pp_init_polynomial(t_norm, KV)
    %
    %   DO m=1,3
    %      CALL ropp_pp_regression(KV, r_leo(:,m), coeff_rleo(:,m))
    %      CALL ropp_pp_regression(KV, r_gns(:,m), coeff_rgns(:,m))
    %   ENDDO
    %
    % !-------------------------------------------------------------------------------
    % ! 5. Regression
    % !-------------------------------------------------------------------------------
    %
    %   DO i=1,n
    %      DO m=1,3
    %         CALL ropp_pp_polynomial(coeff_rleo(:,m),t_norm(i),xleo(i,m),vleo(i,m))
    %         CALL ropp_pp_polynomial(coeff_rgns(:,m),t_norm(i),xgns(i,m),vgns(i,m))
    %      ENDDO
    %
    %      vleo(i,:) = vleo(i,:)/ (time(n) - time(1))
    %      vgns(i,:) = vgns(i,:)/ (time(n) - time(1))
    %
    %   ENDDO
    %
    %   IF (PRESENT(abl)) THEN
    %      abl(:,:) = coeff_rleo(:,:)
    %   ENDIF
    %
    %   IF (PRESENT(abg)) THEN
    %      abg(:,:) = coeff_rgns(:,:)
    %   ENDIF
    %
    % !-------------------------------------------------------------------------------
    % ! 6. Clean up
    % !-------------------------------------------------------------------------------
    %
    %   DEALLOCATE(t_norm)
    %   DEALLOCATE(KV)
    %
    % END SUBROUTINE ropp_pp_satellite_velocities
    
    
end

% ROPP_PP_SATELLITE_VELOCITIES This function generate satellite velocities
% Input
%   time   : 1*n
%   r_leo  : 3*n
%   r_gns  : 3*n
%
% Output
%   xleo   : 3*n
%   vleo   : 3*n
%   xgns   : 3*n
%   vgns   : 3*n
%   abl    : n*3, in each column, the first component being the highest
%            order
%   abg    : n*3
%
% Dependency Functions
%
% Called in
%   C:\Users\bhan002\Dropbox\matlab\7_3 CT2\ct2.m
%
% Created by    Creation Date
%   HAN BO      Jul. 01, 2015
%
% Modified by  Modification Date
%   HAN BO      Jul. 21, 2015
%
% Reference
%   ROPP-PP module
%   The input and output satisfying HB manner
%
% Copyright 2012-2016 HAN BO
% $ Revision: 1.0 $

time  = time' ; % Change to n*1
r_leo = r_leo'; % Change to n * 3
r_gns = r_gns';
nv    = 5;
coeff_rleo = zeros ( nv+1, 3);
coeff_rgns = zeros ( nv+1, 3);

% -------------------------------------------------------------------------------
%  2. Initialisation
% -------------------------------------------------------------------------------
n      = length(time);
t_norm = zeros ( n,1 );
% KV     = zeros ( n, nv+1);

% -------------------------------------------------------------------------------
%  3. Re-normalise time
% -------------------------------------------------------------------------------
for i = 1:n
    t_norm(i,1) = ( time (i) - time (1) )/ ( time(n) - time (1));
end

% -------------------------------------------------------------------------------
%  4. Calculate regression coefficients
% -------------------------------------------------------------------------------
for m=1:3
     [ coeff_mine_rleo(:,m) ] = poly_regress( t_norm, r_leo(:,m), nv );
     [ coeff_mine_rgns(:,m) ] = poly_regress( t_norm, r_gns(:,m), nv );
  
%     [ coeff_mine_rleo(:,m) ] = polyfit ( t_norm, r_leo(:,m), nv ); 
%     Cal coeffs use MATLAB intrinsic function
%     [ coeff_mine_rgns(:,m) ] = polyfit ( t_norm, r_gns(:,m), nv ); % Cal coeffs
%     
end

% KV = zeros ( n, nv + 1);
% [ KV ] = ropp_pp_init_polynomial ( t_norm, KV );
% for m = 1:3
%     [ coeff_mine_rleo(:,m) ] = ropp_pp_regression ( KV, r_leo(:,m) );
%     [ coeff_mine_rgns(:,m) ] = ropp_pp_regression ( KV, r_gns(:,m) );
% end
% % Change the first number to the highest order
% coeff_mine_rleo = coeff_mine_rleo ( end:-1:1, : );
% coeff_mine_rgns = coeff_mine_rgns ( end:-1:1, : );

% -------------------------------------------------------------------------------
%  5. Regression
% -------------------------------------------------------------------------------
for i=1:n
    for m=1:3
        [ xleo(i,m), vleo(i,m) ] = ...
            ropp_pp_polynomial ( coeff_mine_rleo(:,m)', t_norm(i));
        [ xgns(i,m), vgns(i,m) ] = ...
            ropp_pp_polynomial ( coeff_mine_rgns(:,m)', t_norm(i));
    end
    vleo(i,:) = vleo(i,:)/ (time(n) - time(1));
    vgns(i,:) = vgns(i,:)/ (time(n) - time(1));
end

if nargout > 4
    abl = coeff_mine_rleo;
    abg = coeff_mine_rgns;
end

xleo = xleo';
vleo = vleo';
xgns = xgns';
vgns = vgns';

% Evaluate the fitness

% Coeff
% SS_leo_resid = sum ( (xleo - r_leo').^2, 2 );
% SS_leo_total = ( length(xleo(1,:)) - 1) * var ( r_leo', 0, 2);
% rsq_leo      = 1 - SS_leo_resid./SS_leo_total;


end





