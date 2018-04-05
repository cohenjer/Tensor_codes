function [Hes,Jac,grad,f,e]=construct_jacobian_IF_int(L,T,S,beta)
% Builds the Jacobian matrix for the NLFD to run the LM algorithm, details
% in article CILS 2015-2016

A=L{1};
B=L{2};
C=L{3};
[K,R]=size(A);
L=size(B,1);
M=size(C,1);
I=K*L*M;
J=(K+L+M)*R;
Jac= zeros(I,J);
for i=1:I
    m=floor((i-1)/(K*L))+1;
    x=i-(m-1)*K*L;
    l=floor((x-1)/K)+1;
    k=x-(l-1)*K;   
    mp=floor((m-1)/beta)+S;
    ms=min(mp+1,L);
    alpha=(beta+beta*floor((m-1)/beta)+1-m)/beta;
    S1=sum(A(k,:).*B(l,:).*C(m,:));
    if m<=(L-S)*beta+1
        S2=sum(A(k,:).*(B(l,:)+alpha*B(mp,:)+(1-alpha)*B(ms,:)));  
        for r=1:R
            Jac(i,(k-1)*R+r)=(B(l,r)*C(m,r)-(B(l,r)+alpha*B(mp,r)+(1-alpha)*B(ms,r))*S1)*exp(-S2);
            Jac(i,K*R+(l-1)*R+r)=(A(k,r)*C(m,r)-A(k,r)*(1+(l==(mp))*alpha+(l==(ms))*(1-alpha))*S1)*exp(-S2);
            %if l~=mp
                Jac(i,K*R+(mp-1)*R+r)=-(l~=(mp))*A(k,r)*S1*alpha*exp(-S2);
            %end
            %if l~=ms
                Jac(i,K*R+(ms-1)*R+r)=-(l~=(ms))*A(k,r)*S1*(1-alpha)*exp(-S2);
            %end
            Jac(i,(K+L)*R+(m-1)*R+r)=(A(k,r)*B(l,r))*exp(-S2);
        end
    else
        S3=sum(A(k,:).*B(l,:));
        for r=1:R
            Jac(i,(k-1)*R+r)=(B(l,r)*C(m,r)-(B(l,r))*S1)*exp(-S3);
            Jac(i,K*R+(l-1)*R+r)=(A(k,r)*C(m,r)-A(k,r)*S1)*exp(-S3);
            Jac(i,(K+L)*R+(m-1)*R+r)=(A(k,r)*B(l,r))*exp(-S3);
        end
    end
end
Hes=Jac'*Jac;
[f,e]=cout_int(T,A,B,C,S,beta);
grad=Jac'*e;



