function [A_aux,B_aux,C_aux,K_aux,err_anls_vec]=     D_fastgrad(Y,iter,A_0,B_0,C_0,D,opts)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,K_aux,err_als_vec]=      D_fastgrad(Y,iter,A_0,B_0,C_0,D)
%
% Nesterov's Fast Gradient for 3-way CP nonnegative
% decomposition with known dictionnary for second mode. The projection can
% be either alternated for l1 and l2 norms or simultaneous.
% Factors A is normalized with the L_2 norm at each iterate.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - R           : number of components;
% - A_0,C_0     : initial parameters for the factors;
% - B_0         : initial dictionarried factor;
% - D           : known dictionary for factor B (normalized).
% - opts        : options. 
%       opts{1} = 'NN' for non-negative updates
%       opts{2} = 'true' for initializing nnls with previous estimate (NN case only).
%
% Outputs:
% 
% - A_aux       : first factor;
% - B_aux       : second factor;
% - K_aux       : atoms indexes;
% - C_aux       : third factor;
% - err_als_vec : reconstruction errors
%
% External functions:
%
% - kr.m        : Khatri-Rao product;
% - col_norm    : l2 normalisation of columns;
% - simplex_proj: projection on the probability simplex;
% - MatchP      : finds the best atoms in the dictionnary.
%
% List of updates                 -     04/01/2017   J.E.Cohen
%                                       Creation of the file 
%                                 -     08/02/2017   J.E.Cohen
%                                       Added initialisation options
%                                 -     01/03/2017   J.E.Cohen
%                                       introduced K as output variable
%-------------------------------------------------------------------------%

%------------------ l1 Regularization Parameter --------------------------%
lambda_max    =     0.7;
%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     100;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     10^-4;%10^-4
stop_err  =     eps;
% Size of the tensor
dim       =     size(Y);
% Rank of the tensor
R         =     size(A_0,2);
% Index for n-mode matricization
J         =     prod(dim)./dim;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
% Iterates for NNLS
inner     =     10;
% Number of inner iterations
iter_inner=     3; 
%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%
% Three different unfoldings of the tensor
Y_mode    =     cell(3);
% Data tensor unfolding
Y_mode{1} =     reshape(permute(Y,[1,2,3]),dim(1),J(1));  
Y_mode{2} =     reshape(permute(Y,[2,1,3]),dim(2),J(2));
Y_mode{3} =     reshape(permute(Y,[3,1,2]),dim(3),J(3));
% Storing auxiliary variables
A_aux     =     A_0;
B_aux     =     B_0;
C_aux     =     C_0;
C_aux     =     C_aux.*repmat(sqrt(sum(B_aux.^2)),size(C_aux,1),1);
B_aux     =     col_norm(B_aux);
K_aux     =     MatchP(B_aux,D);
S_aux     =     0.2*rand(size(D,2),R);
for j=1:R
S_aux(K_aux(j),j)     =     1;
end
%S_aux     =     rand(size(D,2),R);

% Filter for the relative error decrease
diff_err        =     10^10;
% Initial squared norm of the error
res             =     Y_mode{1}-(A_aux)*transpose(kr(C_aux,B_aux));
err_anls_old    =     (res(:)')*res(:);
% Squared norm of Y
snorm_Y         =     Y_mode{1}(:)'*Y_mode{1}(:);
% Vector of squared norms of the error
%err_anls_vec    =     zeros(1,iter+1);
% First error
err_anls_vec(1) =     err_anls_old;
fprintf('\n Alternating Fast Gradient \n')
fprintf('----\n')

% Accelerated computations

Gram = D'*D;
[~,stepD,~]=svds(Gram,1);
DT   = D'*Y_mode{2};
%stepS_old = 0.01;

% Others
alpha_A_new = 0.5;
alpha_S_new = 0.5;
alpha_C_new = 0.5;

%-------------------------------------------------------------------------%

%--------------------------Algorithm iterate------------------------------%
t         =     0;

while     (diff_err>stop_diff &&  t<iter ...
          &&  err_anls_old>stop_err   )
    % Iteration counter (counts/iter_err)
    t=    t+1;



    % Gradient iterates
    for   i           =     1:iter_err
    
lambda = lambda_max*((t-1)*iter_err+i)/(iter*iter_err);       

%------------------------- S gradient ------------------------------------%
        
        temp_CA        =  (C_aux'*C_aux).*(A_aux'*A_aux);
        stepB = svds(temp_CA,1); 
        stepS =  stepD*stepB;       
        
        temp_DT       =  DT*kr(C_aux,A_aux);
        
        for j=1:iter_inner
            
        S_old    =  S_aux;
        GradS    = -temp_DT+Gram*S_aux*temp_CA+lambda*ones(size(S_aux));
    
            % Direction and acceleration
            
                % Step safety choice
                GradS_NN = GradS; 
                GradS_NN(GradS_NN<=0) = -10^10; % will not matter in next step
                scores = S_aux./GradS_NN;
                step_ast = min(max(scores));
                step = min(1/stepS, step_ast-10^-12);
                
        S_aux    =  max(S_aux - step*GradS,0);
        S_aux    =  col_norm(S_aux);
          
        
        alpha_S_old = alpha_S_new;
        alpha_S_new = 1/2*(-alpha_S_old^2+sqrt(alpha_S_old^4+4*alpha_S_old^2));
        beta_S      = alpha_S_old*(1-alpha_S_old)/(alpha_S_old^2+alpha_S_new);
        S_aux       = S_aux + beta_S*(S_aux-S_old);       
        
        
        end
        
        %[~,K_aux]  =  max(S_aux);
        B_aux      =  D*S_aux;
%        --------- Visualisation -------
%        pause(0.05)
%        spy(S_aux)

%------------------------ A gradient -------------------------------------%
    
         Y_BC          =  Y_mode{1}*kr(C_aux,B_aux);
         temp_BC       =  (C_aux'*C_aux).*(B_aux'*B_aux);


if strcmp(opts{1},'NN')
    if strcmp(opts{2},'true')
        A_aux             =     nnlsHALSupdt(Y_mode{1}',kr(C_aux,B_aux),A_aux',inner)';
    else
        A_aux             =     nnlsHALSupdt(Y_mode{1}',kr(C_aux,B_aux))';
    end
else
        A_aux             =     Y_BC/temp_BC;
end
        A_aux             =     col_norm(A_aux);
        
 %       end  
        



    % C gradient
        
%         temp_AC       =  (C_aux'*C_aux).*(A_aux'*A_aux);
%         [~,stepC,~]   =  svds(temp_AC,1);
%         
%         for j=1:min(i,iter_inner)
%         
%         C_old = C_aux;
%         GradC    = -Y_mode{3}*kr(B_aux,A_aux)+C_aux*((B_aux'*B_aux).*(A_aux'*A_aux));
%     
%         % Direction and acceleration
%     
%         C_aux    =  max(C_aux - 1/stepC*GradC,0);
%         
%         alpha_C_old = alpha_C_new;
%         alpha_C_new = 1/2*(-alpha_C_old^2+sqrt(alpha_C_old^4+4*alpha_C_old));
%         beta_C      = alpha_C_old*(1-alpha_C_old)/(alpha_C_old^2+alpha_C_new);
%         C_aux       = C_aux + beta_C*(C_aux-C_old); 
%         
%         end
%       

if strcmp(opts{1},'NN')
    if strcmp(opts{2},'true')
    C_aux             =     nnlsHALSupdt(Y_mode{3}',kr(B_aux,A_aux),C_aux',inner)';
    else
    C_aux             =     nnlsHALSupdt(Y_mode{3}',kr(B_aux,A_aux))';
    end
else
    Y_AB          =  Y_mode{3}*kr(B_aux,A_aux);
    temp_AB       =  (B_aux'*B_aux).*(A_aux'*A_aux);
    C_aux         =  Y_AB/temp_AB;
end
    
    end

    
    % Squared norm of the error (of the previous iterate)
    err_anls          =     snorm_Y-2*sum(sum(Y_BC.*A_aux))+...
                            sum(sum((A_aux'*A_aux).*temp_BC));
    
    % Evaluation of the relative error decrease
    diff_err          =     abs((err_anls-err_anls_old)/err_anls);
    
    % Storing the last error
    err_anls_old      =     err_anls;
    err_anls_vec(t+1) =     err_anls;
    
    if mod(t*iter_err,iter_err)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_anls,t*iter_err)
    end
    
    if t==(iter-1)
        inner = 50;
        iter_inner=1;
        fprintf('final acceleration\n')
    end
    
end

    
   [~,K_aux]  =  max(S_aux);
   B_aux      =  D(:,K_aux);

%-------------------------------------------------------------------------%
end