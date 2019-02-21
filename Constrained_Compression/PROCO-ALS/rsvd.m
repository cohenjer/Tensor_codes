function  [U,S,V]     =     rsvd(A,k,i)
%-------------------------------------------------------------------------%
% [U,S,V]     =     rsvd(A, k, i)
%
% Randomized SVD. This algorithm implements the randomized SVD presented in
% 
% Halko, N.; Martinsson, P. G.; & Tropp, J. A.; Finding structure with ran-
% domness: Probabilistic algorithms for constructing approximate matrix de-
% compositions. Arxiv preprint arXiv:0909.4061. 2009.
%
% The inputs are:
% 
% - A     : input matrix.
% - k     : number of singular vectors.
% - i     : number of power iterates. (0 for no power iterates)
%
% The outputs are:
% 
% - U     : first factor.
% - S     : second factor.
% - V     : third factor.
%
% List of updates                 -     02/09/2014  -     R C. Farias
%                                       Creation of the file 
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
[m,n]     =     size(A);
transp    =     0; 
% Use conjugate of A if m<n
if  m     <     n
    A     =     A';
    transp=     1;
    [m,n] =     size(A);
end
% Number of columns of the random matrix (can be increased)
l   =     k+2;
%-------------------------------------------------------------------------%

%----------------------------------RSVD-----------------------------------%
% Random matrix
G   =     randn(n,l);
% Matrix to be QR decomposed
H   =     zeros(m,(i+1)*l);
H(:,1:l)        =     A*G;
if  i     >     0
    for   j=1:i
          H(:,j*l+1:(j+1)*l)=     A*(A'*H(:,(j-1)*l+1:(j)*l));     
    end
end
% QR decomposition of H
[Q,~]     =     qr(H,0);
clear     H;
% Compressed matrix to evaluate the SVD
T = A'*Q;
% SVD of T
[V_tilde,S_tilde,W] =     svd(T,'econ');
% Left and right singular vectors 
U_tilde   =     Q*W(:,1:k);
V_tilde   =     V_tilde(:,1:k);
% Singular values
S   =     S_tilde(1:k,1:k);
% Exchange U and V if m<n
if  transp==    1
    V     =     U_tilde;
    U     =     V_tilde;     
    else
    U     =     U_tilde;
    V     =     V_tilde;
end
%-------------------------------------------------------------------------%

