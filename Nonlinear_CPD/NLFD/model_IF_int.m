function [T,Tl,Tif]=model_IF_int(L,nu,S,beta);
% Generate a tensor following the NLFD model from factors in L, with
% intensity nu, shift S and rate between spectral samplings beta.

A=L{1};
B=L{2};
C=L{3};

Tl=tenseur({A,B,C});
I=size(A,1);
J=size(B,1);
K=size(C,1);
R=size(A,2);
Tif=zeros(I,J,K);
for r=1:R
    for i=1:I
        for j=1:J
            for k=1:K 
                
                %if k<=(J-S)*beta+1
                 if k<=(J-S)*beta  
                    kp=floor((k-1)/beta)+S;
                    ks=min(kp+1,J);
                    alpha=(beta+beta*floor((k-1)/beta)+1-k)/beta;
                    Tif(i,j,k)=Tif(i,j,k)+A(i,r)*(B(j,r)+alpha*B(kp,r)+(1-alpha)*B(ks,r));
                else
                    Tif(i,j,k)=Tif(i,j,k)+A(i,r)*B(j,r);
                end
            end
        end
    end
end
Tif=exp(-nu*Tif);
T=Tl.*Tif;
%T=Tl;
