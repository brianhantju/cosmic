function [ FS, DS ] = ropp_pp_sliding_polynomial( T,S,W,np )
if 0
    % $Id: ropp_pp_sliding_polynomial.f90 2228 2009-09-01 15:36:10Z frhl $
    
    %****s* FFT/ropp_pp_sliding_polynomial *
    %
    % NAME
    %    ropp_pp_sliding_poly - Least-square fitting polynomial in sliding
    %                           windows
    %
    % SYNOPSIS
    %    call ropp_pp_sliding_polynomial(t, s, w, np, fs, ds)
    %
    % DESCRIPTION
    %    Least-square fitting polynomial in sliding windows
    %
    % INPUTS
    %    real(wp)             :: t       Time
    %    real(wp), dim([:],:) :: s       Signal samples ([channel],time)
    %    integer, [dim(:)],   :: w       Window width [npoints]
    %    integer              :: np      Polynomial degree
    %
    % OUTPUT
    %    real(wp), dim([:],:), optional  :: fs      Filtered signal
    %    real(wp), dim([:],:), optional  :: ds      Signal derivative
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
    %
    %-------------------------------------------------------------------------------
    % 1. 1-dimensional signal
    %-------------------------------------------------------------------------------
    
    % SUBROUTINE ropp_pp_sliding_poly_1d(t, s, w, np, fs, ds)
    %
    %     USE typesizes, ONLY: wp => EightByteReal
    %     USE ropp_pp_utils
    %
    %     % 7.1 Declarations
    %
    %     IMPLICIT NONE
    %
    %     REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: T  % Time
    %     REAL(wp), DIMENSION(:), TARGET, INTENT(in) :: S  % Signal samples (time)
    %     INTEGER,                  INTENT(in)  :: W  % Window width [samples]
    %     INTEGER,                  INTENT(in)  :: NP % Polynomial degree
    %     REAL(wp), DIMENSION(:),   INTENT(out) :: FS % Filtered signal
    %     REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: DS % Signal derivative
    %
    %     INTEGER     :: NS       % Number of signal samples
    %     INTEGER     :: WS       % Odd sliding window
    %     INTEGER     :: i        % Sample index
    %     INTEGER     :: Imin     % Lower limit of sliding window
    %     INTEGER     :: Imax     % Upper limit of sliding window
    %
    %     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K   % Basic polynomials
    %     REAL(wp), DIMENSION(:), ALLOCATABLE   :: A   % Polynomial coefficients
    %
    %     % 7.2 Initialization
    %
    %     % 7.2.1 Determination of dimensions
    %
    %     NS = SIZE(S)
    %
    %     % 7.2.2 Adjustment of odd filter width
    %
    %     WS = MIN(2*(W/2) + 1, NS)
    %
    %     % 7.2.3. Memory allocation
    %
    %     ALLOCATE(K(WS,0:NP))
    %     ALLOCATE(A(0:NP))
    %
    %     % 7.3 Sliding polynomial regression
    %
    %     DO i = 1, NS
    %
    %        % 7.3.1 Positioning sliding window
    %
    %        Imin = MAX(1,  i - WS/2)
    %        Imax = MIN(NS, i + WS/2)
    %        IF (Imin == 1) THEN
    %           Imax = WS
    %        END IF
    %        IF (Imax == NS) THEN
    %           Imin = NS - WS + 1
    %        END IF
    %
    %        % 7.3.2 Computation of basic polynomials
    %
    %        CALL ropp_pp_init_polynomial(T(imin:imax)-T(i), K)
    %
    %        % 7.3.3 Sliding polynomial regression
    %
    %        CALL ropp_pp_regression(K(:,:), S(imin:imax), A(:))
    %
    %        FS(i) = A(0)
    %
    %        IF (PRESENT(DS)) THEN
    %           DS(i) = A(1)
    %        END IF
    %
    %     END DO
    %
    %     % 7.4 Clean up
    %
    %     DEALLOCATE(K)
    %     DEALLOCATE(A)
    %
    %   END SUBROUTINE ropp_pp_sliding_poly_1d
end
%ROPP_PP_SLIDING_POLYNOMIAL This function do sliding polynomial
% Input
%   T   : time, 1*n
%   S   : signal, 1*n
%   W   : window width, integer
%   np  : polynomial degree, integer
%
% Output
%   FS  : filtered signal, 1*n
%   DS  : signal derivative, 1*n
%
% Dependency Functions
%   ropp_pp_init_polynomial
%   ropp_pp_regression
%
% Called in
%   C:\Users\bhan002\Dropbox\matlab\7_3 CT2\ct2.m
%
% Created by    Creation Date
%   HAN BO      Jul. 22, 2015 verified using ROPP original code in ubuntu
%
% Modified by  Modification Date
%   HAN BO      Jul. 22, 2015
%
% Reference
%   ROPP PP module
%
% Copyright 2012-2016 HAN BO
% $ Revision: 1.0 $

% 1-dimensional signal

if  size (T,1) == 1 && size (S,1) == 1 && size (W,1) == 1 && size (W,2) == 1
    
    
    % 7.2.1 Determination of dimensions
    NS = length (S);
    
    % 7.2.2 Adjustment of odd filter width
    WS = min ( 2 * floor(W/2) + 1, NS );
    
    % 7.2.3. Memory allocation
    K = nan ( WS, 1+np);
    A = nan ( 1, 1+np);
    
    % 7.3 Sliding polynomial regression
    for i = 1:NS
        
        % 7.3.1 Positioning sliding window
        Imin = max ( 1, i - floor(WS/2));
        Imax = min ( NS, i + floor(WS/2));
        
        if ( Imin == 1)
            Imax = WS;
        end
        
        if ( Imax == NS )
            Imin = NS - WS + 1;
        end
        
        % 7.3.2 Computation of basic polynomials
        [ K ] = ropp_pp_init_polynomial( T(Imin:Imax)-T(i), K);
        
        % 7.3.3 Sliding polynomial regression
        [ A ] = ropp_pp_regression( K, S(Imin:Imax) );
        
        
        % ++++++ Use Matlab function +++++
        % Use Matlab function
        %     [ polynomial_coeff ] = polyfit ( T(Imin:Imax)-T(i), S(Imin:Imax), np ); % Cal coeffs
        
        FS(i) = A(1);
        
        if nargout > 1
            DS (i) = A(2);
        end
    end
    
    
elseif size(T,1) ==1 && size(S,1) == 1 && size ( W,1) == 1 && size ( W,2 ) > 1
    
    % 1-dimensional signal ( vector window widths)
    % 9.2.1 Determination of dimensions
    NS = length (S);
    
    % 9.2.3. Memory allocation
    A = nan ( 1, 1+np );
    WS = nan ( 1,NS);
    
    % 9.2.4 Adjustment of odd filter width
    for i = 1:NS
        WS(i) = min ( 2 * floor(W(i)/2) + 1, NS );
    end
    
    % 9.3 Sliding polynomial regression
    for i = 1:NS
        
        % 9.3.1 Positioning sliding window
        Imin = max (1,  i - floor (WS(i)/2) );
        Imax = min ( NS, i + floor(WS(i)/2) );
        if (Imin == 1)
            Imax = WS(i);
        end
        if (Imax == NS)
            Imin = NS - WS(i) + 1;
        end
        
        % 7.3.2 Computation of basic polynomials
        K = nan (WS(i), 1+np );
        [ K ] = ropp_pp_init_polynomial ( T(Imin:Imax) - T(i), K );
        
        % 7.3.3 Sliding polynomial regression
        [ A ] = ropp_pp_regression ( K, S(Imin:Imax) );
        
        % Alternative - Use matlab function 
%         [ A_matlab ] = polyfit ( T(Imin:Imax) -T(i), S(Imin:Imax), np);
%         A_matlab = A_matlab ( end:-1:1);
        
        FS(i) = A(1);
        if nargout > 1
            DS(i) = A(2);
        end
    end
    
end

end
