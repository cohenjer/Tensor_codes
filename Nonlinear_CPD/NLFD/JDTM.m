function [A,it,r1,TD]=JDTM(TM,thresh,itmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   JDTM algorithm (Luciani, Abera) that performs the Joint EigenValue Decomposition (JEVD)
%   of a set of real square matrices stored in the three way array TM.
%   
%   JDTM algorithm has been described in:
%
%   Xavier Luciani, Laurent Albera, Canonical Polyadic Decomposition based on joint eigenvalue decomposition,
%   Chemometrics and Intelligent Laboratory Systems, Volume 132, 15 March 2014, Pages 152-167, ISSN 0169-7439,
%   http://dx.doi.org/10.1016/j.chemolab.2013.12.009.
%   (http://www.sciencedirect.com/science/article/pii/S0169743913002396)
%     
%   Please cite this reference article if you use the JDTM algorithm for your work
%   and see it for further details about the method.
%
%
%   Other input parameters:
%
%   thresh is a threshold value used to stop the JEVD
%   process, usually comprized between 1e-3 and 1e-6.
%
%   itmax is the maximal number of iterations allowed to perform the JEVD
%   process (useful if the threshold is not reached), usually there is no need to choose itmax greater than
%   20.
%
%
%   Onput parameters:
%
%   A is the estimation of the eigenvertor matrix.
% 
%   it is the actual number of iterations used during the JEVD process.
% 
%   r1 stores the value of the cost function to be minimized during the
%   JEVD process at each iteration.
% 
%   TD is a three way array of the same size of TM that stores the
%   diagonilized matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=size(TM,1);
N=size(TM,2);
A=eye(N);
it=1;
TD=TM;
m=0;
for k=1:K
    m=m+(sum(sum((squeeze(TM(k,:,:))-diag(diag(squeeze(TM(k,:,:))))).^2)));%/((sum(sum(squeeze(TM(k,:,:))))).^2);
end
r1(1)=m;
encore=1;
while it<=itmax && encore
    for i=1:N-1
        for j=i+1:N
            [G]=findGopt(TD,i,j);
            for k=1:K
                TD(k,:,:)=inv(G)*squeeze(TD(k,:,:))*G;
            end
            A=A*(G);
            [H]=findHopt(TD,i,j);
            for k=1:K
                TD(k,:,:)=inv(H)*squeeze(TD(k,:,:))*H;
            end
            A=A*(H);
        end
    end
    it=it+1;
    m=0;
    for k=1:K
        m=m+(sum(sum((squeeze(TD(k,:,:))-diag(diag(squeeze(TD(k,:,:))))).^2)));%/((sum(sum(squeeze(TD(k,:,:))))).^2);
    end
    r1(it)=m;
    encore=abs(r1(it-1)-r1(it))/r1(it)>thresh ;
end


function [G,s]=findGopt(TM,i,j)
K=size(TM,1);
N=size(TM,2);
C=zeros(K,2);
for k=1:K
    C(k,:)=[TM(k,i,i)-TM(k,j,j) (TM(k,i,j)+TM(k,j,i))];
end
R=C'*C;
[V,D]=eig(R);
[p]=find(diag(D)==max(diag(D)));
v=V(:,p);
teta2=atan(v(2)/v(1));
c=cos(teta2/2);
s=sin(teta2/2);
G=eye(N);
G(i,i)=c;
G(i,j)=-s;
G(j,i)=s;
G(j,j)=c;

function [H,s]=findHopt(TM,i,j)
K=size(TM,1);
N=size(TM,2);
C=zeros(K,2);
for k=1:K
    C(k,:)=[(TM(k,i,i)-TM(k,j,j)) (TM(k,i,j)-TM(k,j,i))];
end
J=diag([-1 1]);
R=J*(C'*C);
[V,D]=eig(R);
[p]=find(diag(D)==max(diag(D)));
v=V(:,p);
phi2=atanh(v(1)/v(2));
c=cosh(phi2/2);
s=sinh(phi2/2);
H=eye(N);
H(i,i)=c;
H(i,j)=s;
H(j,i)=s;
H(j,j)=c;



