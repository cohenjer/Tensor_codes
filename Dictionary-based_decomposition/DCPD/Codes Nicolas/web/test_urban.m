% test on the Urban data set (hyperspectral image, 162 bands, 307x307 pixels)

load('Urban2.mat'); 
% The image can be downloaded from
% https://sites.google.com/site/nicolasgillis/code 
% It is a 307x307x162 tensor 
for i = 1 : 162
    X(i,:) = reshape(R(:,:,i),1,307^2); 
end
clear R; 
r = 6; 

%  SPA 
Ks = FastSepNMF(X,r); 
Uspa = X(:,Ks); 
Vspa = nnlsHALSupdt(X,Uspa,[],500);  
espa = norm(X-Uspa*Vspa,'fro')/norm(X,'fro')*100; 
% d-SPA
[K,V,U] = NMFdico(X,X,r,30,Ks); 
Udic = X(:,K); 
Vdic = nnlsHALSupdt(X,Udic,V,500);  
edic = norm(X-Udic*Vdic,'fro')/norm(X,'fro')*100; 
% Display results  
fprintf('-----------------------------------\n')
fprintf('Relative error with SPA   = %2.2f %%\n', espa)
fprintf('Relative error with d-SPA = %2.2f %%\n', edic)
fprintf('-----------------------------------\n')
figure;
subplot(1,2,1); plot(Uspa); title('SPA spectral signatures'); 
subplot(1,2,2); plot(Udic); title('d-SPA spectral signatures'); 

affichage(Vspa',3,307,307); title('Abundance maps for SPA'); 
affichage(Vdic',3,307,307); title('Abundance maps for d-SPA'); 