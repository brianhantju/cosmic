clc
% c  = [ 3,2,3,5,6 ];
c =  [ 4,5,3];
c1 = [3,5,4];
x  = [ 4];


y = polyval(c,x);
yp = polyval ( c1,x);

k = polyder ( c1 );
yp2 = polyval ( k, x);

for x = 1:10
[ p, Dp ] = ropp_pp_polynomial( c, x );
disp ([num2str(p)])
disp ([num2str(Dp)])
end




% x = 4;
% % fx = 3 + 2*x + 3*x^2 + 5*x^3 + 6*x^4
% 
% fx = 3*x^2 + 5*x^1 + 4

% fx = 3*x^2+ 2*x +1