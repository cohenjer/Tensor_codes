clear all
close all
clc


load('X.mat')
% 
% T=reshape(T,301*299*41,1);
% 
% temp = T(isnan(T));
% T(isnan(T))   =   0.01*rand(size(temp));
% 
% T   =   permute(reshape(T,299,301,41),[1,3,2]);
% 
% 
%          T=T(100:200,1:2:41,1:5:301);
% 
%          beta=2;
%          ex=[250:5*beta:450];
%          em=[300:5:600];
% 
%          S=(em(1)-ex(1))/(5*beta)+1;
% 
% 
%          [Lf,res,r]=LM_IF_int(T,3,S,beta,50,'b',ex,em);

F = 4;
close all
beta=2.5;
em=em(3:end);
T=T(:,:,3:end);
S=(em(1)-ex(1))/(2*beta)+1;
[Lf,res,r]=LM_IF_int(T,F,S,beta,50,'b',ex,em);