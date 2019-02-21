function [Le,res,errA,it,Sv]=dec_diagN(T,K,P,thresh,itmax,Lv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Computes the rank K truncated Canonical Polyadic Decomposition (CPD or
%   CANDECOMP/PARAFAC) of a real tensor T (of any order), using the DIAG algorithm (Luciani,
%   Abera). This version relies on the JDTM algorithm (Luciani, Abera) to
%   perform the Joint EigenValue Decomposition (JEVD) process.
% 
%     DIAG algorithm has been described in:
%
%     Xavier Luciani, Laurent Albera, Canonical Polyadic Decomposition based on joint eigenvalue decomposition,
%     Chemometrics and Intelligent Laboratory Systems, Volume 132, 15 March 2014, Pages 152-167, ISSN 0169-7439,
%     http://dx.doi.org/10.1016/j.chemolab.2013.12.009.
%     (http://www.sciencedirect.com/science/article/pii/S0169743913002396)
%     
%     Please cite this reference article if you use the DIAG algorithm for your work
%     and see it for further details about the method.
%
%   Other input parameters:
%
%   P is an integer value (0 < P < Q-1, where Q is the tensor order), it controls the choice 
%   of the unfolding matrix used to compute the SVD (first step of the algorithm).
%   Indeed DIAG algorithm imposes some restrictions on the tensor dimensions to
%   work:
%
%   In order to make the algorithm work, you must choose a permutation of the tensor dimensions and a
%   value of P so that the CPD rank is not greater than: 
% 
%   1) the product of the P first tensor dimensions and 
%   2) the product of the dimensions P+1 to Q-1 (for a Q-order tensor).
% 
%       For instance: 
% 
%   Third order tensors, Q = 3.
%   Here we have necessarily P = 1 hence at least two of the tensor dimensions must be greater or equal
%   to the CPD rank. 
% 
%   Fourth order tensors, Q = 4. 
%   Here we can choose either P = 1 or P = 2 but the condition remains the same in both cases and is simply:
%   at least one tensor dimension is greater than the CPD rank and 
%   at least one product of two of the remaining dimensions is also greater than the CPD rank. 
% 
%   Fifth order tensors, Q = 5. 
%   Here 1 ≤ P ≤ 3:
%   if we choose P = 1 or P = 3 then  at least one tensor dimension is greater than the CPD rank and 
%   at least one product of three of the remaining dimensions is also greater than the CPD rank.
%   if we choose P = 2 then at least one product between two tensor dimensions and
%   another product between two of the remaining dimensions are greater than the CPD rank.
%
%   Choose P=0 (not always optimal but recommended) to make the algorithm choose automaticaly a suitable
%   permutation and value of P.
%
%   thresh is a threshold value used to stop the JEVD
%   process, usually comprized between 1e-3 and 1e-6.
%
%   itmax is the maximal number of iterations allowed to perform the JEVD
%   process (useful if the threshold is not reach). With the JDTM algorithm used
%   in this version, usually there is no need to choose itmax greater than
%   20.
%
%   Lv is an optional parameter. If they are known, actual loading matrices can be given in the cell Lv and
%   used in order to solve permutation and scaling ambiguity. Usefull to test
%   the algorithm on known data sets.
%
%
%   Output parameters:
%
%   Le is a cell array which stores the estimated loading matrices.
%
%   res is the relative root squared reconstruction error term (root value of the objective
%   function to be minimized).
%
%   If Lv is given then vector errA contains the relative root squared error,
%   computed between estimated and actual loading matrices.
%
%   it is the actual number of iterations used during the JEVD process
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I=size(T);
Q=length(I);
ind=1:Q;

if P~=0
M=prod(I(1:P));
N=prod(I(P+1:Q-1));
if K > min([M N]);
  mess='\n In order to make the algorithm work, you must choose a permutation of the tensor dimensions and a value of P \n so that  the CPD rank is not greater than: \n \n 1) the product of the P first tensor dimensions and \n 2) the product of the dimensions P+1 to Q-1 (for a Q-order tensor). \n \n  For instance: \n \n Third order tensors, Q = 3.\n Here we have necessarily P = 1 hence at least two of the tensor dimensions must be greater or equal to the CPD rank. \n \n Fourth order tensors, Q = 4. \n Here we can choose either P = 1 or P = 2 but the condition remains the same in both cases and is simply:\n at least one tensor dimension is greater than the CPD rank and \n at least one product of two of the remaining dimensions is also greater than the CPD rank. \n \n Fifth order tensors, Q = 5. \n Here 1 ≤ P ≤ 3:\n if we choose P = 1 or P = 3 then  at least one tensor dimension is greater than the CPD rank and \n at least one product of three of the remaining dimensions is also greater than the CPD rank.\n if we choose P = 2 then at least one product between two tensor dimensions and \n another product between two of the remaining dimensions are greater than the CPD rank.';
  fprintf(mess)
  return
end
else 
    [~,po]=min(I);
    po=po(1);
    ind(Q)=po;
    ind(po)=Q;
    indc=ind(1:end-1);  
    Tindc=perms(indc);
    mref=0;
    for i=1:size(Tindc,1)
        for p=1:Q-2
            Ip=I([Tindc(i,:) po]);
            m=min([prod(Ip(1:p)) prod(Ip(p+1:Q-1))]);
            if m>mref
                P=p;
                indc=Tindc(i,:);
                mref=m;
            end
        end
    end
ind(1:end-1)=indc;
I=I(ind);
T=permute(T,ind);
M=prod(I(1:P));
N=prod(I(P+1:Q-1));
end



Td=reshape(T,[M N*I(Q)]);
[U,S,V]=svds(Td,K);
diag(S)
Sv=S;
X=S*V';
TM=zeros((I(end)-1)*I(end)/2,K,K);
k=1;
for k1=1:I(end)-1
    for k2=k1+1:I(end)
        G1=X(:,(k1-1)*N+1:k1*N);
        G2=X(:,(k2-1)*N+1:k2*N);
        TM(k,:,:)=pinv(G1')*G2';
        k=k+1;
    end
end
[R,it]=JDTM(TM,thresh,itmax);
E=U/(R');
F=V*S*R;
L=cell(P);
for i=1:K
    if P~=1
        a=E(:,i);
        Z=reshape(a,I(1:P));
        [S,W]=comp_hosvd(Z,1);
        k=1;
        for j=1:P;
            L{j}(:,i)=W{k};
            k=k+1;
        end
        L{j}(:,i)=L{j}(:,i)*S;
    else
        L{1}=E;
    end
    if N~=length(I)-1;
        cb=F(:,i);
        Z=reshape(cb,I(P+1:end));
        [S,W]=comp_hosvd(Z,1);
        k=1;
        for j=P+1:Q;
            L{j}(:,i)=W{k};
            k=k+1;
        end
        L{j}(:,i)=L{j}(:,i)*S;
    else
        L{P+1}=F;
    end
end
Ee=L{P};
for i=P-1:-1:1
    Ee=pkr(Ee,L{i});
end
Fe=L{Q};
for i=Q-1:-1:P+1
    Fe=pkr(Fe,L{i});
end

T1=tenseur(L);
res1=(T(:)-T1(:));
res1=res1'*res1./(T(:)'*T(:));
%res=sqrt(res1);
res=res1
Le=cell(Q);
for i=1:Q
    Le{ind(i)}=L{(i)};
end

errA=[];
if nargin==nargin('dec_diagN');
    Le=permscal3(Le,Lv);
    errA=zeros(1,Q);
    for i=1:length(I)
        A=Le{i};
        Av=Lv{i};
        err=(A(:)-Av(:));
        errA(i)=norm(err)/norm(Av(:));
    end
end



function T=tenseur(L)
% create tensor T of dimensions dim_T and rank rang_T from the matrices
% stored in the cell array L
ordre_T=length(L);
rang_T=size(L{1},2);
dim_T=zeros(1,ordre_T);
for i=1:ordre_T
    dim_T(i)=size(L{i},1);
end
T=zeros(dim_T);
for f=1:rang_T
    TT=ones(dim_T);
    for d=1:ordre_T
        v1=ones(1,ordre_T);
        v1(d)=dim_T(d);
        v2=dim_T;
        v2(d)=1;
        TT=repmat(reshape(L{d}(:,f),v1),v2).*TT;
    end
    T=TT+T;
end


function [Ln,Lv]=permscal3(L,Lv)
% Removes the permutation and the scaling indeterminancy in the loading
% matrices stored in L from the real matrices stored in Lv 
N=length(L);
F=size(L{1},2);
X=[];
Xv=[];
for n=1:N
    A=L{n};
    An=zeros(size(A));
    Avn=zeros(size(A));    
    Av=Lv{n};
    for i=1:F
        s=sum(A(:,i));
        An(:,i)=A(:,i)/s;
        s=sum(Av(:,i));
         Avn(:,i)=Av(:,i)/s;
    end
    X=[X;An];
    Xv=[Xv;Avn];
    clear An Avn
end
X=abs(X);
Xv=abs(Xv);
P=perms(1:F);
r=zeros(1,size(P,1));
for i=1:size(P,1)
    r(i)=norm(Xv-X(:,P(i,:)));
end
[~,pos]=min(r);
p=P(pos,:);
Ln=L;
for n=1:N
    A=L{n}(:,p);
    Av=Lv{n};
    Ae=zeros(size(A));
    for i=1:F
        Ae(:,i)=A(:,i)*A(:,i)'*Av(:,i)/(A(:,i)'*A(:,i));
    end
    Ln{n}=Ae;
    clear Ae
end


function [S,W]=comp_hosvd(X,rg)
% S=comp_hosvd(X) computes the High-Order Singular Value Decomposition
% (HOSVD) of an N-th order tensor X of dimensions I_1 x I_2 x ... x I_N.
%
% The outputs parameters: 
% S: is an N-th array of dimensions R_1 x R_2 x ... x R_N, 
% where R_1, R_2, ..., R_N are the mode-1, mode-2, ..., mode-N ranks of the 
% unfolded matrix representations X1, X2, ..., X3 of the tensor X, respectively.
% W: is a cell-array of N elements, the n-th element is the mode-n
% eigenmatrix. 
%
% Principle: From S and W{1},...,W{N} we can rebuild the original tensor X.

I= size(X);
N=length(I);
% hosvd
S=X;
W=cell(N);
for n=1:N 
    % n-th mode unfolded matrix
    if I(1)>=rg
    Sd= reshape(S,I(1),prod(I(2:end))); 
    % SVD of Sd
    [W{n},~]=svds(Sd,rg);
    Sd= W{n}'*Sd;
    S= reshape(Sd,[rg I(2:end)]);
    % circularly permute dimensions of tensor S
    S=shiftdim(S,1);
    I=size(S);
    else
        W{n}=[];
    end
end

function Z=pkr(A,B)
F=size(A,2);
Z=zeros(size(A,1)*size(B,1),F);
for f=1:F
    Z(:,f)=kron(A(:,f),B(:,f));
end


