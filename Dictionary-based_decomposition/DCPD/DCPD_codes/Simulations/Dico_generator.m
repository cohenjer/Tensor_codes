function [ D,S,fun_sin ] = Dico_generator( type, l, k, nb, shift, mixing,k2 )
%-------------------------------------------------------------------------%
% [ output_args ] = Dico_generator( type, l, k, mixing)
%
% Generate a dictionnary for simulation purposes. A working condition is
% that nb<l.
%
% The inputs are:
% 
% - type        : Type of dictionnary ('sinusoids','splines','exponentials');
% - l           : (length) Size of the atoms;
% - k           : number of atoms;
% - nb          : number of points in the non-zero part of each atom;
% - shift       : shifts [shift1,shift2] for adjacent atoms and various frequencies, choose nb for no overlapping;
% - mixing      : 'mixed' for mixing the atoms of the dictionnary using
% - k2          : if mixing is used, sets the new number of atoms.
% sparse coefficients.
%
% The outputs are:
% 
% - D           : Dictionnary.
%
% List of updates                 -     10/10/2016  -  Jeremy E. Cohen
%                                       Creation of the file
%-------------------------------------------------------------------------%

switch(type)
    case('sinusoids')
    
% Square functions bank
fun_sin   =     zeros(l,k);
for i=1:k
    if mod(i*shift(1),l)==0
    fun       =     cos(pi/(1+floor(i*shift(1)/l))*(-1:2/nb:(1-2/nb)));    
    fun_sin(1:floor(nb)/2,i)    =     fun(nb-floor(nb/2)+1:end)';
    else
fun       =     cos(pi/(1+floor(i*shift(1)/l))*(-1:2/nb:(1-2/nb)));
fun_sin(mod(i*shift(1),l):min(mod(i*shift(1),l)+nb-1,l),i)    =     fun(1:min(nb,l-mod(i*shift(1),l)+1))';
    end
end

case('exponentials')
    
% Square functions bank
fun_sin   =     zeros(l,k);


fun       =    (1+0.1*floor(1*shift(1)/l))*exp(-((-3:6/nb:(3-6/nb)).^2)*(1+0.1*floor(1*shift(1)/l)^2));    
fun_sin(1:nb,1)    =     fun';   

for i=1:k-1
    shift2=floor(i*shift(1)/l)*shift(2);
    if mod(i*shift(1)+shift2,l)==0 
    fun       =    1/(1+0.1*floor(i*shift(1)/l))*exp(-((-3:6/nb:(3-6/nb)).^2)/(1+0.1*floor(i*shift(1)/l)^2));    
    fun_sin(min(shift2,l):(min(shift2+nb,l)-1),i+1)    =     fun(1:(min(shift2+nb,l)-min(shift2,l)))';
    else
fun       =     1/(1+0.1*floor(i*shift(1)/l))*exp(-((-3:6/nb:(3-6/nb)).^2)/(1+0.1*floor(i*shift(1)/l))^2);
fun_sin(mod(shift2+i*shift(1),l):min(mod(i*shift(1)+shift2,l)+nb-1,l),i+1)    =     fun(1:min(nb,l-mod(i*shift(1)+shift2,l)+1))';
    end
end

end

if strcmp(mixing,'mixed')==1
    S = max(randn(k,k2),0);
    D = fun_sin*S;
elseif strcmp(mixing,'correl')
    Gram   =   fun_sin'*fun_sin;
    [~,indices]   =   sort(Gram,1);
    S      =   zeros(k,k2);
    used   =   [];
    for j=1:floor(k2/2)
        n1 = 0;
        n2 = 1;
        while ismember(indices(k-n1,j),used)
            n1=n1+1;
            n2=n2+1;
        end
        while ismember(indices(k-n2,j),used)
            n2=n2+1;
        end
        S(indices(k-n1,j),2*j-1)     =     1;
        S(indices(k-n2,j),2*j)       =     1;
        used = [used,indices(k-n1,j),indices(k-n2,j)];
    end
    D   =   fun_sin*S;
else
    D = fun_sin;
    S = eye(k);
end

D         =     D.*repmat(1./sqrt(sum(D.^2)),l,1);

end

