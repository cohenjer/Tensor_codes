% Select treshold to split the entries of x into two subsets
% 
% See Section 3.2 in 
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 

function [thres,delta,fobj] = fquad(x,s); 

if nargin == 1
    s = 0.01; % grid for the values of delta
end

[fdel,fdelp,delta,finter,gs] = fdelta(x,s); 
% fdel is the percentage of values smaller than delta
% finter is the number of points in a small interval around delta


warning('off'); 
fobj = -log( fdel.* (1-fdel) ) + exp(finter); 
% Can potentially use other objectives: 
%fobj = -log( fdel.* (1-fdel) ) + 2.^(finter); 
%fobj = ( 2*(fdel - 0.5) ).^2 + finter.^2; 
%fobj = -log( fdel.* (1-fdel) ) + finter.^2; 
%fobj = ( 2*(fdel - 0.5) ).^2 + finter.^2; 
warning('on');
[a,b] = min(fobj); 
thres = delta(b); 


% Evaluate the function fdel = sum( x_i <= delta)/n and its derivate 
% for all delta in interval [0,1] with step s
function [fdel,fdelp,delta,finter,gs] = fdelta(x,s); 

n = length(x); 
if nargin == 1
    s = 0.01;  
end
delta = 0:s:1; 
lD = length(delta); 

gs = 0.05; % Other values could be used, in [0,0.5]

for i = 1 : lD
    fdel(i) = sum(x <= delta(i))/n;
    if i == 2 % use only next point to evaluate fdelp(1)
        fdelp(1) = (fdel(2)-fdel(1))/s; 
    elseif i >= 2 % use next and previous point to evaluate fdelp(i)
        fdelp(i-1) = (fdel(i)-fdel(i-2))/2/s; 
        if i == lD % use only previous point to evaluate fdelp(lD)
            fdelp(lD) = (fdel(lD)-fdel(lD-1))/s; 
        end
    end 
    deltahigh = min(1,delta(i) + gs);
    deltalow = max(0,delta(i) - gs); 
    finter(i) = ( sum(x <= deltahigh) - sum(x < deltalow) )/n/(deltahigh-deltalow) ;
end