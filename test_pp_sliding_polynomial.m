clc
clear 

T = [-2,-1,3,4,6,7,20];
S = [12,3,42,32,1,13,43];
W = 3;
np = 2;


[ FS, DS ] = ropp_pp_sliding_polynomial( T,S,W,np );