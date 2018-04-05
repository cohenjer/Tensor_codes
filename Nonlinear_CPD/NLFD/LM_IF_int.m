function [Lf,res,r,Lp]=LM_IF_int(T,R,S,beta,niter,col,ex,em,Lv);

% Compute the Non Linear Fluorescence Decomposition (NLFD) of a fluorescence
% tensor containing a set of FEEM. The NLFD takes into account possible 
% inner filter effects.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T: Fluorescence data tensor. First mode is concentrations mode, second
% mode is for excitations and third mode is for emissions.
%
% R: Number of fluorophores/chromatophores
%
% S: Shift between excitation and emission: 
% S=1+(lambda_em_min-lambda_ex_min)/delta_ex where delta_ex is the
% excitation sampling step
%
% beta is the ratio between the excitation sampling step and the emission
% sampling step: beta=delta_ex/delta_em.
%
% niter: Maximum number of iterations of the algorithm (stopping
% criterion).
%
% Plotting parameters: col (color of the error plot). ex: vector of
% sampled exitation wavelengths, em: vector of the sampled emission 
% wavelengths of the FEEM.
%
% Lv is an optional parameter. If they are known, actual loading matrices 
% can be given in the cell Lv and used in order to solve permutation and 
% scaling ambiguity. Usefull to test the algorithm on known data sets.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lf is a cell array containing the three NLFD estimated loading matrices.
%
% res is the relative root mean squared reconstruction error term 
% (squared root of the objective function to be minimized).

% If Lv is given then vector r contains the relative root mean squared error
% computed between estimated and actual loading matrices.
%
% Lp is a cell array containing the three unconstrained CPD (PARAFAC) 
% estimated loading matrices DIAG algorithm is used).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VERY IMPORTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The algorithm will work only if
% beta is an INTEGER greater than 0 and if lambda_em_min coincides with a
% sampled value of lambda_ex. 
% For instance, we can have the following discrete excitation range: 
% [280 : 20 : 500] nm and discrete emission range: [340 : 5 : 600] nm 
% since here beta=4 (20/5) and lambda_em_min (340 nm) coincides 
% with the fourth wavelength of the excitation range.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=T/norm(T(:));
r=[];
rp=r;
Lp=dec_diagN(T,R,0,1e-5,30); %PARAFAC/CP decomposition of T
for i=1:3
    Lp1{i}=abs(Lp{i});
end
Le=Lp1; % Initial guess of the NLFD
[Hess,Jac,grad,f_old]=construct_jacobian_IF_int(Le,T,S,beta);
mu= 2;
nu=2;
it=0;
f_vec=[f_old; zeros(niter,1)];
f_rel_vec=zeros(niter,1);
flag=0;
% loop for estimating the NLFD loading matrices ( using the Levenberg Marquardt
% method)
while   it<niter %&& flag==0   %&& err_new >=1e-20  %&& err_rel>1e-15 %&& mu> 1e-8  && lda> 1e-8 &&flag==0
    it=it+1 ;
    [Le,f_new,mu,nu,Hess,grad,Jac]=boucle_lm_IF_int(Le,T,S,mu,nu,Hess,grad,Jac,beta);
    clc;home;
    disp(['Cost function: ',num2str(2*f_new/norm(T(:)))]);
    disp(['Iterations: ',num2str(it)]);
    crit=abs(f_new-f_old)/abs(f_old);
    if crit<eps && crit~=0
        flag=1;
    end
    f_old=f_new;
    f_rel_vec(it)= crit;
    f_vec(it+1)=f_new;
end
res=sqrt(f_new/(T(:)'*T(:)))
% error plot (evolution of the reconstruction error according the
% iterations)
figure
semilogy(f_vec,col)
Ap=Lp{1};
Bp=Lp{2};
Cp=Lp{3};
Ae=Le{1};
Be=Le{2};
Ce=Le{3};
%  NLFD (blue) and CP/parafac (green) loading plots 
figure
if nargin==nargin('LM_IF_int');
    Av=Lv{1};
    Bv=Lv{2};
    Cv=Lv{3};   
    [Ae,Be,Ce]=corrige(Ae,Be,Ce,Lv{1},Lv{2},Lv{3});
    for i=1:R
        rA(i)=sum((Ae(:,i)-Av(:,i)).^2)/sum((Av(:,i).^2));
        rB(i)=sum((Be(:,i)-Bv(:,i)).^2)/sum((Bv(:,i).^2));
        rC(i)=sum((Ce(:,i)-Cv(:,i)).^2)/sum((Cv(:,i).^2));
    end
    rA=mean(rA);
    rB=mean(rB);
    rC=mean(rC);
    r=[rA rB rC];
    [Ap,Bp,Cp]=corrige(Ap,Bp,Cp,Lv{1},Lv{2},Lv{3});
    for i=1:R
        rA(i)=sum((Ap(:,i)-Av(:,i)).^2)/sum((Av(:,i).^2));
        rB(i)=sum((Bp(:,i)-Bv(:,i)).^2)/sum((Bv(:,i).^2));
        rC(i)=sum((Cp(:,i)-Cv(:,i)).^2)/sum((Cv(:,i).^2));
    end
    rA=mean(rA);
    rB=mean(rB);
    rC=mean(rC);
    rp=[rA rB rC];
else
    [Ap,Bp,Cp]=corrige(Ap,Bp,Cp,Ae,Be,Ce);
end
for k=1:R
    subplot(3,R,k)
    plot(Ae(:,k),'LineWidth',2)
    hold on
    plot(Ap(:,k),'g','LineWidth',2)
    if nargin==nargin('LM_IF_int');
        plot(Av(:,k),'r','LineWidth',2)
        set(gca, 'YLim', [0 max(Av(:,k))+0.1*max(Av(:,k))])
    end
    
    title(['component ' num2str(k)]);
    xlabel('Sample number')
    ylabel('Relative concentration')
end
for k=1:R
    subplot(3,R,k+R)
    plot(ex,Be(:,k),'LineWidth',2)
    hold on
    plot(ex,Bp(:,k),'g','LineWidth',2)
    if nargin==nargin('LM_IF_int');
        plot(ex,Bv(:,k),'r','LineWidth',2)
        %r(2)=norm(Be(:)-Bv(:))./norm(Bv(:));
        set(gca, 'YLim', [0 max(Bv(:,k))+0.1*max(Bv(:,k))])
    end
    xlabel('Excitation wavelength')
    ylabel('Intensity (a.u.)')
end
for k=1:R
    subplot(3,R,k+2*R)
    plot(em,Ce(:,k),'LineWidth',2)
    hold on
    plot(em,Cp(:,k),'g','LineWidth',2)
    if nargin==nargin('LM_IF_int');
        plot(em,Cv(:,k),'r','LineWidth',2)
        set(gca, 'YLim', [0 max(Cv(:,k))+0.1*max(Cv(:,k))])
    end
    xlabel('Emission wavelength')
    ylabel('Intensity (a.u.)')
end
if nargin==nargin('LM_IF_int');
    legend('NLFD','CPD','Reference')
else
    legend('NLFD','CPD')
end
rp
r
Lf=Le;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)




%
% Lp=Lini;
% % %solve permutation and scaling ambiguity
%    % Lf=permscal2(Le,Lv);
%     %Lp=permscal2(Lp,Lv);
%      %figure
%      k1=1;
%     for i=1:3
%        % rf(i)=norm(Lf{i}-Lv{i},2)/norm(Lf{i},2);
%        % rp(i)=norm(Lp{i}-Lv{i},2)/norm(Lf{i},2);
%        figure
%         k=num2str(k1);
%         %subplot(3,3,k1)
%         %plot(Lv{i}(:,1))
%         hold on
%         plot(Lf{i}(:,1),'r')
%        % plot(Lp{i}(:,1),'g')
%         %k1=k1+1;
%         k=num2str(k1);
%        % subplot(3,3,k1)
%         %plot(Lv{i}(:,2))
%         figure
%         hold on
%         plot(Lf{i}(:,2),'b')
%         %plot(Lp{i}(:,2),'g')
%         k1=k1+1;
%        % subplot(3,3,k1)
%         %plot(Lv{i}(:,2))
%         figure
%         hold on
%         plot(Lf{i}(:,3),'g')
%         %plot(Lp{i}(:,2),'g')
%         k=num2str(k1);
%     end
%
% %fit plot
%
% %plot(err_vec(1:it),col)
% hold on
%  [Te,Tl,Tif]=model_IF(Le,1,dec);
% % res=sum((T(:)-Te(:)).^2)
% % resn=sum((T(:)-Te(:)).^2)/sum(T(:).^2)
% % figure
% % imagesc(squeeze(T(1,:,:)))
% % colorbar
% % figure
% % imagesc(squeeze(Te(1,:,:)))
% % colorbar

