%------------------------------Hard_tests--------------------------------%
function out = Hard_tests( input )
% This file tests the Dictionnary CPD model in difficult scenarios.
%
% How to use: call the following
%
%           out = Hard_tests(input)
%
% where
% input : {N, [R,Re], [d,c], sigma, pos, corr}
% where - N is the number of simulation points
%       - R is the rank, Re the estimated rank
%       - [d,c] is [number of atoms, number of classes]
%       - sigma is noise std
%       - pos is '' for normal, 'NN' for nonnegative factors 
%       - corr is the coefficient ro in C2 = ro*C2+(1-ro)C1. 0 is complete
%       colinearity of two first columns of C and 1 is stochastic
%       independance.
%
% A dictionnary with d=1000, c=20 and fixed K has been preset and is used
% by default. To override it and use another dictionary, the code has to be
% modified.
%
% External functions:
%
% - kr.m              : Khatri-Rao product;
% - D_anls.m          : ANLS algorithm;
% - amb_correct.m     : Matching Permutation for two CP models.
%
% List of updates                 -     23/02/2016  -     J. E. Cohen 
%                                       Creation of the file
%                                 -     19/01/2017  -     J. E. Cohen 
%                                       Modified tests to match TSP
%                                       Submission
%                                 -     08/03/2017  -     J. E. Cohen
%                                       Added SMPALS, time, unconstrained
%                                       case.
%                                 -     24/09/2017  -     J. E. Cohen
%                                       Made the script into a function,
%                                       added variable correlation of
%                                       columns of C.
%-------------------------------------------------------------------------%

% ---- How to ----
% For gridding on delta, use delta*N in functions flexMP and DADMM and fix
% the initialization. Improvement from 0 to 87*0.005 in delta
% For correlation, set S = [eye(R),[1,0...];...]
% For missing values, uncomment according lines.

%--------------Cleaning the workspace and loading functions---------------%
%clc;
%clear     all;
%close     all;
%-------------------------------------------------------------------------%

%-----------------------------Model parameters----------------------------%
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS
addpath D:\Travail\Travail\Toolbox\Other_Toolboxes\nway331
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\Continuous
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\Codes_Nicolas
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\Data_Dico
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\Greedy
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\Matrix
addpath D:\Travail\Travail\Toolbox\Post-doc-Perso\Dictionnary_ALS\TOOLS

% 3-way block dimensions
%dim       =     [60 50 7]; % large
%dim       =     [20 50 20]; normal
dim       =     [20 50 7]; % small
% Number of trials
N         =     input{1}; %100
% Number of components
R         =     input{2}(1); %10
% Number of atoms and classes
d         =     input{3}(1); c=input{3}(2); % 1000 20
% Standard deviations of the measurements noise
sigma_n   =     input{4}; %0.01
% Standard deviation of the distance to dictionnary
sigma_d   =     0;
% Number of iterations
iter      =     1000;
%-------------------------------------------------------------------------%
%Stock     =     zeros(4,N,3);

% Recollecting stocked Dictionary and index set K (for reproducibility)
Dmem       =     load('Sim_Unconstrained_TSP','D','Ktrue');%Dico_gen_group(dim(2),d/c,c,'NN');
%Dmem       =     load('Sim_Unconstrained_TSP_Uncorr','D','Ktrue');
D          =     Dmem.D;
% For setting a random dictionary
%D          =     Dico_gen_group(dim(2),d/c,c,'NN');

for i=1:N




%for i=1:N

%-------- For correlating group atoms, comment if normal use --------%
% alpha(i)=0.05*(1-exp(-3*(-1+i/N))/exp(3));
% %Corr_matr    = ones(K,K);
% D = col_norm(D0*(eye(K)-alpha(i)*ones(K,K)));



%------------------------------Tensor generation--------------------------%
% Random factors generation
if strcmp(input{5},'NN')
A         =     abs(randn(dim(1),R)); % for non-negative
C         =     abs(randn(dim(3),R)); %for non-negative
else
A         =     randn(dim(1),R);
C         =     randn(dim(3),R);
end
% Correlating
%C(:,2)    =     input{6}*C(:,2) + (1-input{6})*C(:,1); % change one column
C         =     C*(input{6}*eye(R) + (1-input{6})*ones(R,R)/R);  

% Selecting the true atoms
Ktrue     =     50*linspace(0,R-1,R)+1; %or d/c*...
%Ktrue     =     10*linspace(0,R-1,R)+1; %or d/c*...
B         =     D(:,Ktrue)+sigma_d*randn(dim(2),R);

% For other selection method of simulated picked atoms
% Ktrue     =     randperm(d,R);
% Ktrue     =     linspace(1,R^2-R+1,R);
% Ktrue(2)=Ktrue(1);
% Ktrue     =     Dmem.Ktrue(1:R);     

% Normalization
A         =     A.*repmat(1./sqrt(sum(A.^2)),dim(1),1);
C         =     C.*repmat(1./sqrt(sum(C.^2)),dim(3),1);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%Noisy measurements
% Noisy tensor
Y           =     reshape(A*transpose(kr(C,B)),dim)+...
                  sigma_n*randn(dim);
              
SNR(i)         =     10*log10(sqFnorm(reshape(A*transpose(kr(C,B)),dim))/sqFnorm(Y-reshape(A*transpose(kr(C,B)),dim)));

% Missing values (not implemented as of 19/01/2017)
% percentage = 0;
% % mask of missing data
% M    =    false(prod(dim),1); % No missing data
% M(randperm(prod(dim),floor(prod(dim)*percentage/100))) =    true; % Random missing data
% 
% Y(M) =    nan;
% 
% Y    =    reshape(Y,dim);

% for i=1:N

fprintf('\n--------- Iteration number %d ----------\n',i)    


%-----------------------------Initialization------------------------------%

Re          =     input{2}(2); %R+2

% Non-negative case
if strcmp(input{5},'NN')
A_0         =     abs(randn(dim(1),Re));
K_0         =     randperm(d,Re);
A_0         =     A_0.*repmat(1./sqrt(sum(A_0.^2)),dim(1),1);
C_0         =     abs(randn(dim(3),Re));

else
% Unconstrained case
A_0         =     randn(dim(1),Re);
K_0         =     randperm(d,Re);
A_0         =     A_0.*repmat(1./sqrt(sum(A_0.^2)),dim(1),1);
C_0         =     randn(dim(3),Re);
end

%-------------------------------------------------------------------------%

%--------------------------   Plain  NALS  -------------------------------%

%N-Way toolbox
Oldload{1} = A_0;Oldload{2} = D(:,K_0);Oldload{3} = C_0; Options(1) = 10^-16; Options(5)=100; Options(6)=iter;
tic
if strcmp(input{5},'NN')
[Fact,~, ~,~,~,err_als]         =     parafac(Y,Re,Options,[2 2 2]);%,Oldload); % nonnegative
else
[Fact,~, ~,~,~,err_als]         =     parafac(Y,Re,Options,[0 0 0]);%,Oldload);
end
t_als(i)=toc;
A_anls = Fact{1}; B_anls=Fact{2}; C_anls=Fact{3};
C_anls  = col_norm(C_anls);
A_anls  = col_norm(A_anls);
B_anls  = col_norm(B_anls);

%--------------------------Dictionary NALS updates----------------------------%

opts{1} = input{5}; opts{2}= 'true';
                                                           
%for i=1:N %for grid on delta

% MPNALS 
tic
[A_mpals,B_mpals,C_mpals,K_mpals,err_mpals_vec]     =     D_mpnals(Y,iter,A_anls,B_anls,C_anls,D,opts);
t_mpals(i)=toc;

% MPNALS randomized initialization
tic
[A_mpals_r,B_mpals_r,C_mpals_r,K_mpals_r,err_mpals_vec_r]     =     D_mpnals(Y,5*iter,A_0,D(:,K_0),C_0,D,opts);
t_rmpals(i)=toc;

% Test Flex MPNALS
%[A_ppals,B_ppals,C_ppals,S_ppals,P_est,~,err_ppals_vec]     =     D_P2MP(Y,iter,A_0,S_0,C_0,D);
%[A_ppals,B_ppals,C_ppals,S_ppals,P_est,~,err_ppals_vec]     =     D_P2mpnals(Y,iter,A_anls,MatchP(B_anls,D),C_anls,D);                                                            
tic
[A_ppals,B_ppals,C_ppals,K_ppals,err_ppals_vec]     =     D_mpnals_flex_2(Y,iter,A_anls,B_anls,C_anls,D,0.04,1.1,opts);                                                            
t_flex_mpals(i)=toc;

% Smooth MPALS
tic
[A_smpals,B_smpals,C_smpals,K_smpals,B_flex,err_smpals]     =     D_smpnals(Y,iter,A_anls,B_anls,C_anls,D,1.1,opts);                                                            
t_smpals(i)=toc;


% Fast Gradient
tic
[A_fgrad,B_fgrad,C_fgrad,K_fgrad,err_fgrad]   =  D_fastgrad(Y,2*iter,A_anls,B_anls,C_anls,D,opts);
t_fgrad(i)=toc;

%------------ Postprocessing -------------%

% Overfactoring ?

norm_anls(i,:)             =     sort(sum(C_anls.^2));
norm_ppals(i,:)            =     sort(sum(C_ppals.^2));
norm_mpals(i,:)            =     sort(sum(C_mpals.^2));
norm_mpals_r(i,:)            =     sort(sum(C_mpals_r.^2));
norm_fgrad(i,:)            =     sort(sum(C_fgrad.^2));
norm_smpals(i,:)           =     sort(sum(C_smpals.^2));

C_anls    =     col_norm(C_anls); C_ppals    =     col_norm(C_ppals); C_mpals    =     col_norm(C_mpals);
C_mpals_r = col_norm(C_mpals_r); C_fgrad=col_norm(C_fgrad); C_smpals = col_norm(C_smpals);


% Fractional Classification error

Rmax = max(Re,R);

classif_anls(i)  =  classif_err(Ktrue,MatchP(B_anls,D))/Rmax;
classif_mpals(i) =  classif_err(Ktrue,K_mpals)/Rmax;
classif_mpals_r(i) =  classif_err(Ktrue,K_mpals_r)/Rmax;
classif_ppals(i) =  classif_err(Ktrue,K_ppals)/Rmax;
classif_fgrad(i) =  classif_err(Ktrue,K_fgrad)/Rmax;
classif_smpals(i) =  classif_err(Ktrue,K_smpals)/Rmax;

classif_g_anls(i)  =  classif_err_group(Ktrue,MatchP(B_anls,D),20)/Rmax;
classif_g_mpals(i) =  classif_err_group(Ktrue,K_mpals,20)/Rmax;
classif_g_mpals_r(i) =  classif_err_group(Ktrue,K_mpals_r,20)/Rmax;
classif_g_ppals(i) =  classif_err_group(Ktrue,K_ppals,20)/Rmax;
classif_g_fgrad(i) =  classif_err_group(Ktrue,K_fgrad,20)/Rmax;
classif_g_smpals(i) =  classif_err_group(Ktrue,K_smpals,20)/Rmax;

% Factor error (From this point, external information is provided)

[A_anls, B_anls, C_anls]   =     amb_correct(A_anls, B_anls ,C_anls, A,B,C);
[A_ppals,B_ppals,C_ppals]  =     amb_correct(A_ppals,B_ppals,C_ppals,A,B,C);
[A_mpals,B_mpals,C_mpals]  =     amb_correct(A_mpals,B_mpals,C_mpals,A,B,C);
[A_mpals_r,B_mpals_r,C_mpals_r]  =     amb_correct(A_mpals_r,B_mpals_r,C_mpals_r,A,B,C);
[A_fgrad,B_fgrad,C_fgrad]  =     amb_correct(A_fgrad,B_fgrad,C_fgrad,A,B,C);
[A_smpals,B_smpals,C_smpals]  =     amb_correct(A_smpals,B_smpals,C_smpals,A,B,C);

Stock_anls(i,:)            =     [sqFnorm(A-A_anls)/sqFnorm(A), sqFnorm(B-B_anls)/sqFnorm(B), sqFnorm(C-C_anls)/sqFnorm(C)];
Stock_anls_proj(i,:)       =     [sqFnorm(A-A_anls)/sqFnorm(A), sqFnorm(B-D(:,MatchP(B_anls,D)))/sqFnorm(B), sqFnorm(C-C_anls)/sqFnorm(C)];
Stock_ppals(i,:)           =     [sqFnorm(A-A_ppals)/sqFnorm(A),sqFnorm(B-B_ppals)/sqFnorm(B),sqFnorm(C-C_ppals)/sqFnorm(C)];
Stock_mpals(i,:)           =     [sqFnorm(A-A_mpals)/sqFnorm(A),sqFnorm(B-B_mpals)/sqFnorm(B),sqFnorm(C-C_mpals)/sqFnorm(C)];
Stock_mpals_r(i,:)         =     [sqFnorm(A-A_mpals_r)/sqFnorm(A),sqFnorm(B-B_mpals_r)/sqFnorm(B),sqFnorm(C-C_mpals_r)/sqFnorm(C)];
Stock_fgrad(i,:)           =     [sqFnorm(A-A_fgrad)/sqFnorm(A),sqFnorm(B-B_fgrad)/sqFnorm(B),sqFnorm(C-C_fgrad)/sqFnorm(C)];
Stock_smpals(i,:)          =     [sqFnorm(A-A_smpals)/sqFnorm(A),sqFnorm(B-B_smpals)/sqFnorm(B),sqFnorm(C-C_smpals)/sqFnorm(C)];
% note : better to compare normalized factors



%-------------------------------------------------------------------------%

end
 
% classif_anls  = 100-100*sum(classif_anls)/N;
% classif_mpals = 100-100*sum(classif_mpals)/N;
% classif_mpals_r = 100-100*sum(classif_mpals_r)/N;
% classif_ppals = 100-100*sum(classif_ppals)/N;
% classif_fgrad = 100-100*sum(classif_fgrad)/N;
% classif_smpals = 100-100*sum(classif_smpals)/N;

% OUTS

out.labels  = {'ALS-ALS_projected','Flexible-MPALS','MPALS','random-MPALS','FG','Smooth-MPALS'};
out.stock   = {Stock_anls,Stock_anls_proj,Stock_ppals,Stock_mpals,Stock_mpals_r,Stock_fgrad,Stock_smpals};
out.classif = {classif_anls,classif_ppals,classif_mpals,classif_mpals_r,classif_fgrad,classif_smpals};
out.classif_g = {classif_g_anls,classif_g_ppals,classif_g_mpals,classif_g_mpals_r,classif_g_fgrad,classif_g_smpals};
out.time    = {t_als,t_flex_mpals,t_mpals,t_rmpals,t_fgrad,t_smpals};
out.SNR     = SNR;

