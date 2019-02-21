function  [A_aux,B_aux,C_aux,A_p_aux,B_p_aux,C_p_aux,map_als_vec] = ...
          ccp_als(Y_1,Y_2,iter,H,H_p,sigma_n,sigma_c,type...
                                            ,A_0,B_0,C_0,A_p_0,B_p_0,C_p_0)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,A_p_aux,B_p_aux,C_p_aux,map_als_vec] = ...
%           ccp_als(Y_1,Y_2,iter,H,H_p,sigma_n,sigma_c,type...
%                                           ,A_0,B_0,C_0,A_p_0,B_p_0,C_p_0)
%
% ALS algorithm for Coupled CP decomposition. Factors C and C_p are coupled
% through transformations H and H_p. Factors A, B, A_p and are normalized 
% with the L_2 norm at each iterate.
%
% Inputs:
% 
% - Y_1         : first data matrix;
% - Y_2         : second data matrix;
% - iter        : maximum number of iterates;
% - H_1,H_2     : transformation matrices for the factors.
%                 
%                 For type='joint' the transformation is for vec(C)
%                 If the transformations are given for transformations on 
%                 the matrix factors do 
%                     H=kron(eye(R),H_mat) and H_p=kron(eye(R),H_p_mat);
%
%                 For type='alter' the transformation are on the matrix
%                 factors.
%
%                 For type='direct' the transformations are ignored but an
%                 input is still required.
%
% - sigma_n     : vector with noise standard deviations;
% - sigma_c     : coupling strenght;
% - type        : type of algorithm to be used
%                 'joint'   -     joint coupled factors update;
%                 'alter'   -     alternating coupled factors update;
%                 'direct'  -     joint update with H=H_p=Identity(ignore
%                                 given H and H_p matrices);
% - A_0,...     : initial factors for the first model;
% - A_p_0,...   : initial factors for the second model.  
%
% Outputs:
% 
% - A,...       : factors for the first model;
% - A_p,...     : factors for the second model;
% - map_als_vec : iterations (Row 1) x MAP objective (Row 2).
%
% External functions:
%
% - kr.m        : Khatri-Rao product;
% - sylvester.m : Sylvester equation solver, necessary if type='alter'. 
%                 Built-in function since Matlab 2014.
%
% List of updates                 -     23/02/2015  -     J. E. Cohen and
%                                                         R. C. Farias 
%                                       Creation of the file
%                                 -     22/07/2015  -     R. C. Farias and
%                                                         J. E. Cohen
%                                       Updated with fast right inverse and
%                                       added options to find C factors
%                                       ('type' input)
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the MAP objective
iter_map  =     10;
% Size of the filter for the relative error decrease
n_diff_map=     3;
% Stopping criteria parameters
stop_diff =     eps;
stop_map  =     eps;
% Maximum number of iterations (/iter_map)
iter      =     ceil(iter/iter_map);
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Size of the tensor
dim_1     =     size(Y_1);
dim_2     =     size(Y_2);
R         =     size(A_0,2);
% Index for n-mode matricization
J_1       =     prod(dim_1)./dim_1;
J_2       =     prod(dim_2)./dim_2;

% Three different unfoldings of the tensor
Y_1_mode  =     cell(3);
Y_2_mode  =     cell(3);
% Data tensors unfolding
% First tensor
Y_1_mode{1}     =     reshape(permute(Y_1,[1,2,3]),dim_1(1),J_1(1));  
Y_1_mode{2}     =     reshape(permute(Y_1,[2,1,3]),dim_1(2),J_1(2));
Y_1_mode{3}     =     reshape(permute(Y_1,[3,1,2]),dim_1(3),J_1(3));
% Second tensor
Y_2_mode{1}     =     reshape(permute(Y_2,[1,2,3]),dim_2(1),J_2(1));  
Y_2_mode{2}     =     reshape(permute(Y_2,[2,1,3]),dim_2(2),J_2(2));
Y_2_mode{3}     =     reshape(permute(Y_2,[3,1,2]),dim_2(3),J_2(3));

% Squared norm of the data
snorm_Y_1       =     Y_1_mode{1}(:)'*Y_1_mode{1}(:);
snorm_Y_2       =     Y_2_mode{1}(:)'*Y_2_mode{1}(:);

% Precisions and coupling strength
prec      =     1./(sigma_n.^2);
coupl     =     1/(sigma_c^2);

% Normalization of the MAP objective function function
norm      =     mean([prec coupl]);
% Normalized constants in the MAP objective function
prec      =     prec/norm;
coupl     =     coupl/norm;

% Coupling terms
coupl_11  =     coupl*(H')*H;
coupl_12  =     coupl*(H')*H_p;
coupl_21  =     (coupl_12');
coupl_22  =     coupl*(H_p')*H_p;
%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%
% Filter for the relative error decrease
diff_map_filt   =     1000*ones(1,n_diff_map);
% Storing initial condition on the compressed space
A_aux     =     A_0;
B_aux     =     B_0;
C_aux     =     C_0;
A_p_aux   =     A_p_0;
B_p_aux   =     B_p_0;
C_p_aux   =     C_p_0;
% Initial MAP objective value
res_1     =     Y_1_mode{1}-A_aux*kr(C_aux,B_aux)';
err_1     =     res_1(:)'*res_1(:);
res_2     =     Y_2_mode{1}-A_p_aux*kr(C_p_aux,B_p_aux)';
err_2     =     res_2(:)'*res_2(:);
if              strcmp(type,'alter')
res_c_m   =     H*C_aux-H_p*C_p_aux;
res_c     =     res_c_m(:);
elseif          strcmp(type,'direct')
res_c     =     C_aux(:)-C_p_aux(:);
else
res_c     =     H*C_aux(:)-H_p*C_p_aux(:);
end
err_c     =     res_c'*res_c;
map_als_old     =     prec(1)*err_1+prec(2)*err_2+coupl*err_c;
% Vector of squared norms of the error
map_als_vec     =     zeros(1,iter+1);
map_als_vec(1)  =     map_als_old;     
%-------------------------------------------------------------------------%

fprintf('\n Coupled CP ALS - CCP ALS \n\n')

%--------------------------Algorithm iterate------------------------------%
t         =     0;
while     (mean(diff_map_filt)>stop_diff &&  t<iter ...
          &&  map_als_old>stop_map   )
    % Iteration counter (counts/iter_map)
    t=    t+1;
    % Coupled ALS (10 iterates)
    for   i           =     1:iter_map
    
    % A factors update
    temp_CB           =     (C_aux'*C_aux).*(B_aux'*B_aux);
    A_aux             =     (Y_1_mode{1}*kr(C_aux,B_aux))/temp_CB;
    
    temp_CB           =     (C_p_aux'*C_p_aux).*(B_p_aux'*B_p_aux);
    A_p_aux           =     (Y_2_mode{1}*kr(C_p_aux,B_p_aux))/temp_CB;
    
    % Normalization
    A_aux       =     A_aux.*repmat(1./sqrt(sum(A_aux.^2)),dim_1(1),1);
    A_p_aux     =     A_p_aux.*repmat(1./sqrt(sum(A_p_aux.^2)),dim_2(1),1);
    
    % B factors update
    temp_CA     =     (C_aux'*C_aux).*(A_aux'*A_aux);
    B_aux       =     (Y_1_mode{2}*kr(C_aux,A_aux))/temp_CA;
    
    temp_CA     =     (C_p_aux'*C_p_aux).*(A_p_aux'*A_p_aux);
    B_p_aux     =     (Y_2_mode{2}*kr(C_p_aux,A_p_aux))/temp_CA;
    
    % Normalization
    B_aux       =     B_aux.*repmat(1./sqrt(sum(B_aux.^2)),dim_1(2),1);
    B_p_aux     =     ...
                      B_p_aux.*repmat(1./sqrt(sum(B_p_aux.^2)),dim_1(2),1);
    
    %---------------------------Factor C update---------------------------%
    switch(type)
        
    %----------------------Joint C factors update-------------------------%
          case('joint')
    % Intermediate quantities which will be also used to evaluate the MAP
    % objective
    temp_1_BA   =     (B_aux'*B_aux).*(A_aux'*A_aux);
    temp_2_BA   =     (B_p_aux'*B_p_aux).*(A_p_aux'*A_p_aux);
    Y_1_BA=     Y_1_mode{3}*kr(B_aux,A_aux);
    Y_2_BA=     Y_2_mode{3}*kr(B_p_aux,A_p_aux);
    
    % Vectorized system of equations
    D_1   =     prec(1)*(kron(temp_1_BA,eye(dim_1(3))))+...
                coupl_11;
    D_2   =     prec(2)*(kron(temp_2_BA,eye(dim_2(3))))+...
                coupl_22;
    D     =     [D_1 -coupl_12;-coupl_21 D_2];
    
    % Right-hand side of the vectorized system of equations
    E_1   =     prec(1)*Y_1_BA;
    E_2   =     prec(2)*Y_2_BA;
    E     =     [E_1(:);E_2(:)];
    % Vectorized factors         
    vec_C =     D\E;
    % Unvectorized factors
    C_aux       =      reshape(vec_C(1:dim_1(3)*R),[dim_1(3) R]);
    C_p_aux     =      reshape(vec_C(dim_1(3)*R+1:end),[dim_2(3) R]);
    %---------------------------------------------------------------------%
    
    %-Alternating C factors update - Least squares with Sylvester equation%
          case('alter')
    % Intermediate quantities which will be also used to evaluate the MAP
    % objective
    temp_1_BA   =     (B_aux'*B_aux).*(A_aux'*A_aux);
    temp_2_BA   =     (B_p_aux'*B_p_aux).*(A_p_aux'*A_p_aux);
    Y_1_BA=     Y_1_mode{3}*kr(B_aux,A_aux);
    Y_2_BA=     Y_2_mode{3}*kr(B_p_aux,A_p_aux);
          
    % C factor update
    D     =     coupl_11;
    E     =     prec(1)*temp_1_BA;
    F     =     prec(1)*Y_1_BA+coupl_12*C_p_aux;
    C_aux =     sylvester(D,E,F);
    
    % C_p factor update
    D           =     coupl_22;
    E           =     prec(2)*temp_2_BA;
    F           =     prec(2)*Y_2_BA+coupl_21*C_aux;
    C_p_aux     =     sylvester(D,E,F);
    %---------------------------------------------------------------------%
    
    %-------------Joint C factors update assuming direct coupling---------%
          case('direct')
    % Intermediate quantities which will be also used to evaluate the MAP
    % objective
    temp_1_BA   =     (B_aux'*B_aux).*(A_aux'*A_aux);
    temp_2_BA   =     (B_p_aux'*B_p_aux).*(A_p_aux'*A_p_aux);
    Y_1_BA=      Y_1_mode{3}*kr(B_aux,A_aux);
    Y_2_BA=      Y_2_mode{3}*kr(B_p_aux,A_p_aux);
    
     % Joint system of equations
    D     =     [prec(1)*temp_1_BA+coupl*eye(R) -coupl*eye(R);...
                -coupl*eye(R) prec(2)*temp_2_BA+coupl*eye(R)];
    
    % Right-hand side of the joint system of equations
    E_1   =     prec(1)*Y_1_BA;
    E_2   =     prec(2)*Y_2_BA;
    E     =     [E_1 E_2];
    % Stacked (row) C factors       
    stack_C     =     E/D;
    % Unvectorized sources
    C_aux       =      stack_C(:,1:R);
    C_p_aux     =      stack_C(:,R+1:end);
    %---------------------------------------------------------------------%
          otherwise
          disp('Unknown method. Choose between: joint, alter or direct.')
    end
    %---------------------------------------------------------------------%
    end
    
    % MAP objective value
    err_1       =     snorm_Y_1-2*sum(sum(Y_1_BA.*C_aux))+...
                      sum(sum((C_aux'*C_aux).*temp_1_BA));
    err_2       =     snorm_Y_2-2*sum(sum(Y_2_BA.*C_p_aux))+...
                      sum(sum((C_p_aux'*C_p_aux).*temp_2_BA));
    if                strcmp(type,'alter')
    res_c_m     =     H*C_aux-H_p*C_p_aux;
    res_c       =     res_c_m(:);
    elseif            strcmp(type,'direct')
    res_c       =     C_aux(:)-C_p_aux(:);
    else
    res_c       =     H*C_aux(:)-H_p*C_p_aux(:);
    end
    err_c       =     res_c'*res_c;
    map_als     =     prec(1)*err_1+prec(2)*err_2+coupl*err_c;
    
    % Evaluation of the relative error decrease
    diff_map_filt(1:end-1)  =     diff_map_filt(2:end);
    diff_map_filt(end)      =     abs((map_als-map_als_old)/map_als);
    
    % Storing the last error
    map_als_old        =    map_als;
    map_als_vec(t+1)   =    map_als;
    
    if mod(t*iter_map,100)  ==    0
    % Printing the iterate MAP value and step size
    fprintf('MAP objective value : %g  \tit : %d\n',map_als,t*iter_map)
    end
end
% Deleting zero error trailing at the end
map_als_vec(map_als_vec==0)=[];
map_als_vec=[0 iter_map:iter_map:iter_map*t;map_als_vec];
%-------------------------------------------------------------------------%
end