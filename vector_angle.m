function [ angle ] = vector_angle( x, y )
%VECTOR_ANGLE This function return the vector angle of two vectors
% Input
%   x : 1*3 or 3*1 vector
%   y : the same shape as x 
%
% Output
%    angle : 1*1 real number
%
% Dependency Functions
%  
% Called in
%   C:\Users\bhan002\Dropbox\matlab\7_3 CT2\ct2.m
%
% Created by    Creation Date
%   HAN BO      Jul. 22, 2015
%
% Modified by  Modification Date
%   HAN BO      Jul. 22, 2015
%
% Reference
%  
% Copyright 2012-2016 HAN BO
% $ Revision: 1.0 $

angle = acos (  dot(x,y) / ( norm(x) * norm(y) ));

end

