A = [1,9,1,4;0,10,2,2;0,5,3,2;0,23,12,0];

% MAXLOC
clc

N = size ( A,1);
for k = 1:N-1
    newmat = A(k+1:N,k);
    newidx = find ( newmat ~= 0 ) ;
    newmat2 = newmat ( newidx, k);
    m  = max ( newmat2 );
    if isempty (m)
        m =0;
    end
    
    
    %  [ value, m ] = max (  ( A(k+1:N, k) ) );
    %
    %  if A(k+1:N,k) == 0
    %
    %      return
    %  end
    %
    %
    %  disp ('-----')
    %  disp (  A(k+1:N, k) )
    %
    %
    %  disp (m)
    %  disp ('-----')
end
