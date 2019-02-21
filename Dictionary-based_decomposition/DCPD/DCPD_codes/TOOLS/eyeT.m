function [ T ] = eyeT( R )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

T = zeros(R,R,R);
for i=1:R
    T(i,i,i) = 1;
end
end

