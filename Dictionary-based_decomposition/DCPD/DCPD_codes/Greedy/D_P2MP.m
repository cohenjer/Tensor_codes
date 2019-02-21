function  [A_aux,B_aux,C_aux,K_aux,P,Y_est,err_P2MP] = ...
          D_P2MP(Y,iter,A_0,B_0,C_0,D)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,err_als_vec] = ...
%           parafac2(Y,iter,A_0,S_0,C_0,D)
%
% Preprocessed MPANLS algorithm for Coupled CP decomposition for N matrices. Factors Bk 
% are coupled through estimated factors Pk and without coupling error; the latent variable
% B* is taken from columns of a dictionnary D.
%
% Inputs:
% 
% - Y           : cell of data tensors;
% - iter        : maximum number of iterates;
% - P_0         : cell of initital transformation matrices for the factors
%                 C;
% 
%                 
% - sigma_n     : vector with noise standard deviations;
% - A_0,...     : cell of initial factors;
% - D           : Dictionnary of columns of C*;
%
% Outputs:
% 
% - A,B,C       : cell of estimated factors;
% - P           : cell of estimated coupling matrices;
% - map_als_vec : objective function.
%
% External functions:
%
% - kr.m        : Khatri-Rao product;
%
% List of updates                 -     05/11/2015   J.E.Cohen
%                                       Creation of the file
%                                 -     10/11/2015   J.E.Cohen
%                                       Handling missing data
%                                 -     07/11/2016   J.E.Cohen
%                                       Dictionnary version
%                                 -     01/03/2017   J.E.Cohen
%                                       removed S, introduced K
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%

fprintf('\n Parafac2 as in Bro[98] \n\n')

%-------------------------------Parameters--------------------------------%
dim       =     size(Y);
J         =     prod(dim)./dim;
N         =     dim(3); 
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     100;


err_P2MP     =     zeros(iter,1);

R         =     size(A_0,2);

% Storing initial condition
A_aux     =     A_0;
B_aux     =     B_0; % Coupled mode, in the preprocessed space
C_aux     =     C_0; % Stacked mode, scaling ratios
C_aux     =     C_aux.*repmat(sqrt(sum(B_aux.^2)),size(C_aux,1),1);
B_aux     =     col_norm(B_aux);

% Setting unknown values to 0
vecy      =     Y(:);
M         =     isnan(vecy);
vecy(M)   =     min(vecy);
Y         =     reshape(vecy,dim);


%---- Outer Loop for global iterations ---%

for k=1:iter

%-----------------------------------------------------%


% Slices definition
for i=1:N

slice{i}  =     Y(:,:,i); 

end
 

%-------------------------------------------------------------------------%

%-----------------Coupling Estimation and Preprocessing-------------------%

P{1}      =    eye(dim(2));
Y_pre(:,:,1)   =    Y(:,:,1);

for i=2:N

[u,~,v]   =     svd(B_aux*diag(C_aux(i,:))*A_aux'*slice{i},'econ');
P{i}      =     v(:,1:R)*u(:,1:R)';%v*u';

slice_pre{i}  =     slice{i}*P{i}';
Y_pre(:,:,i)  =     slice_pre{i};

end

Ynorm     =     sum(Y_pre(:).^2);

%--------------------------ALS Steps--------------------------------------%

% Unfoldings1
for i=1:N
Y_mode1     =     reshape(permute(Y_pre,[1,2,3]),dim(1),dim(2)*dim(3));  
Y_mode2     =     reshape(permute(Y_pre,[2,1,3]),dim(2),dim(1)*dim(3));
Y_mode3     =     reshape(permute(Y_pre,[3,1,2]),dim(3),dim(1)*dim(2));
end



%--------------------------NALS one iteration-----------------------------%
    
    % A factors update
    temp_BC     =     (B_aux'*B_aux).*(C_aux'*C_aux);
    A_aux       =     (Y_mode1*kr(C_aux,B_aux))/temp_BC; 
    A_aux(A_aux<0)    =     0;
    
    % Normalization
    A_aux       =     A_aux.*repmat(1./sqrt(sum(A_aux.^2)),dim(1),1);
    
    % B factors update
    temp_AC     =     (A_aux'*A_aux).*(C_aux'*C_aux);
    B_aux       =     (Y_mode2*kr(C_aux,A_aux))/temp_AC; 
    B_aux(B_aux<0)    =     0;
    % Atom selection via Matching Pursuit
    K_aux       =     MatchP(B_aux,D);
    B_aux       =     D(:,K_aux);
    
    % C factors update
    temp_AB     =     (A_aux'*A_aux).*(B_aux'*B_aux);
    Y_AB        =     Y_mode3*kr(B_aux,A_aux);
    C_aux       =     Y_AB/temp_AB;
    C_aux(C_aux<0)    =     0;
    
    if mod(k,iter_err)==0
    
    % Estimation error
    err_P2MP(k)    =     Ynorm-2*sum(sum(Y_AB.*C_aux))+...
                             sum(sum((C_aux'*C_aux).*temp_AB));
    end
    
    % Estimated tensor
    
    Y_est     =     [];

    for i=1:N
        slice_est       =       A_aux*diag(C_aux(i,:))*B_aux'*P{i}';   
        Y_est(:,:,i)    =       slice_est;
    end
        
    % Estimation of the missing data

    vecy_est  =     Y_est(:);
    vecy(M)   =     vecy_est(M);
    Y         =     reshape(vecy,dim);
    
    % Handling NANs
    if sum(sum(isnan(A_aux))) || sum(sum(isnan(B_aux))) || sum(sum(isnan(C_aux)))
    A_aux(isnan(A_aux))     =     1;
    B_aux(isnan(B_aux))     =     1;
    C_aux(isnan(C_aux))     =     1;
    break
    end
    
    if mod(k,iter_err)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_P2MP(k),k)
    end
    
end
    
end