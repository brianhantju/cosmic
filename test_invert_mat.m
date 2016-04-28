

A = [0,9,1,4;0,10,2,2;0,5,3,2;0,23,12,0];
% 

[ B ] = ropp_pp_invert_matrix( A );

B1 = inv ( A);

B - B1
