function [ D ] = Dico_gen_group( L,K,c,NN )
%Generates a dictionnary with c classes of K atoms with linear baseline and 2
%features (sinus cardinal and triangular pulse).

x=0:(L-1);

for i=1:c

a = (2*rand(1)-1)/L; b = 2*rand(1)-1; 

for j=1:K

u_1 = (2*rand(1)-1)/4; u_2 = (2*rand(1)-1)/4; 
x_1 = floor(rand(1)*L)+1; x_2 = floor(rand(1)*L)+1;

if strcmp(NN,'NN')
D(:,j+(i-1)*K) = abs(a*x+b + u_1*sinc(pi/2/3*(x-x_1)) + ...
    u_2*(triangularPulse(x_2-4,x_2-2,x_2,x)-triangularPulse(x_2,x_2+2,x_2+4,x)));
else
D(:,j+(i-1)*K) = a*x+b + u_1*sinc(pi/2/3*(x-x_1)) + ...
    u_2*(triangularPulse(x_2-4,x_2-2,x_2,x)-triangularPulse(x_2,x_2+2,x_2+4,x));
end
end

end

D = col_norm(D);

end

