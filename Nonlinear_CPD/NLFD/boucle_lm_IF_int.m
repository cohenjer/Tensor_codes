function [L_newv,f_new,lda_new,nu,Hess,grad,Jac]=boucle_lm_IF_int(L,T,S,lda_old,nu,Hess,grad,Jac,beta)
% One step of the LM algorithm, with gradient update and computation of the
% cost function. Outputs the new variables, the Jacobian, Hessian and
% Gradient.

L_newv=L;
A=L{1};B=L{2};C=L{3};
N=size(A,2);
[I,J,K]=size(T);
p_old=[vec(A.');vec(B.');vec(C.')];
[f_old,e_old]=cout_int(T,A,B,C,S,beta);
% step;
dp= -(Hess + lda_old*eye((I+J+K)*N))\grad;
% update p-vector
p_new= p_old + dp;
p_new(p_new<0)=0;
% rebuild factor matrices
pA=p_new(1:I*N);
A=reshape(pA,N,I)';
pB=p_new(I*N+1:(J+I)*N);
B=reshape(pB,N,J)';
pC=p_new((J+I)*N+1:(J+I+K)*N);
C=reshape(pC,N,K)';
[f_new,e_new]=cout_int(T,A,B,C,S,beta);
Jp=(Jac*dp)';
f_new_approx = f_old + Jp*e_old + Jp*Jp'/2;
% calculate gain ratio
gain= ((f_old - f_new)/(f_old - f_new_approx));
% update lambda (regularization factor)
if gain >0
    lda_new=lda_old*max(1/3,1-(2*gain-1)^3);
    nu=2;
    L_newv={A,B,C};
    [Hess,Jac,grad]=construct_jacobian_IF_int(L_newv,T,S,beta);
else
    lda_new=lda_old*nu;
    nu=2*nu;
    f_new=f_old;
end


