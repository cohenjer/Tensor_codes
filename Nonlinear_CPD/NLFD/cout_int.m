function [f,e]=cout_int(X,A,B,C,S,beta);
% Computes the cost function ||X-X_e||^2

[Xe]=model_IF_int({A,B,C},1,S,beta);
M=numel(X);
Xd=reshape(X,M,1);
Xde=reshape(Xe,M,1);
e= Xde - Xd;
f=e'*e/2;
