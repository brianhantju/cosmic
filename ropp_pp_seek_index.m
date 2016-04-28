function [ ip ] = ropp_pp_seek_index( x, xp )
%ROPP_PP_SEEK_INDEX 
% !****s* PPUtils/ropp_pp_seek_index *
% !
% ! NAME
% !    ropp_pp_seek_index - Find grid interval containing given point
% !
% ! SYNOPSIS  
% !    ip = ropp_pp_seek_index(x, xp)
% !
% ! DESCRIPTION
% !    Search for grid interval containing the interpolation point.
% !    Combination of Newton-like and dichotomic iterative search.
% !     Result:
% !         >0   - index of grid interval containing xp (x(ip) <= xp <= x(ip+1))
% !          0   - xp is not inside grid
% !         -1   - iterations do not converge
%     REAL(wp), DIMENSION(:), INTENT(in) :: x   ! x-grid (homogeneous)
%     REAL(wp),               INTENT(in) :: xp  ! point to locate inside grid
%     INTEGER                            :: ip
    
%     INTEGER  :: N      ! Number of grid points
%     INTEGER  :: Imin   ! Upper estimate of index
%     INTEGER  :: Imax   ! Lower estimate of index
%     INTEGER  :: Dir    ! Direction of argument change
%     INTEGER  :: It     ! Iteration count
%     INTEGER  :: di     ! Index increment in iterations
%     INTEGER  :: is     ! Step direction count

%     ! 4.1 Grid size and direction calculation
   
    N    = length(x);
% original    Dir  = NINT( SIGN(1.0_wp, x(N)-x(1)) )
    if x(N) - x(1) >=0
        tmp = 1;
    elseif x(N) - x(1) < 0
        tmp = -1;
    end
    Dir = round (tmp);
        
%     ! 4.2 Checking if point is inside grid

% ori:    if ((Dir * xp < Dir * x(1)) .OR. (Dir*xp > Dir*x(N))) THEN
    if ((Dir * xp < Dir * x(1)) || (Dir*xp > Dir*x(N)))
       ip = 0;
%        RETURN
    end

%     ! 4.3 Initial approximation
    
    Imin = 1;
    Imax = N;
% ori:    ip   = Imin+FLOOR(REAL(Imax - Imin, wp)*(xp - x(Imin))/(x(Imax) - x(Imin)))
% is there any size issues?
% REAL: convert to real type
    ip   = Imin + floor (  (Imax - Imin) * (xp - x(Imin)) / (x(Imax) - x(Imin)));
    ip   = MAX(1, MIN(ip, N-1));
    
%     ! 4.4 Iterative index search

    It = 0;
    is = 0;
    
    Search: DO
       if ((Dir*x(ip) <= Dir*xp) .AND. (Dir*xp <= Dir*x(ip+1))) THEN
          EXIT Search
       end
       if (ABS(is) > 1)
          ip = (Imax + Imin)/2
          is = 0
       end
       if (Dir*x(ip+1) < Dir*xp)
          Imin = ip + 1
          di   = FLOOR(REAL(Imax - Imin, wp)*(xp - x(Imin))/(x(Imax) - x(Imin)))
          ip   = Imin + di
          is   = is + 1
       elseif (Dir*xp < Dir*x(ip))
           Imax = ip
           di   = FLOOR(REAL(Imin - Imax, wp)*(xp - x(Imax))/(x(Imin) - x(Imax)))
           ip   = Imax + di
           is   = is - 1
       end
       
       ip = MAX(1, MIN(ip, N-1))
       It = It + 1
       IF (It > N) THEN
       ip = -1;
       
       EXIT Search
       END IF
    END DO Search
    
end