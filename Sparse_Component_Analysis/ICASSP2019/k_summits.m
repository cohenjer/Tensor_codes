function [A,B,aff,err] = k_summits(M,k,iter,A0,B0)
% [A,B,B_aux,err] = k_summits(M,k,iter,A0,B0)
% This function computes the k-sparse matrix factorization of
% data matrix M:
%               M = AB, ||b_i||_0 < k+1
% K-sparisty is imposed on the columns of matrix B, so that
% at most k components are active at each sample along the second mode. 
%
% a.k.a. NOLRAK
%
% Implemented with a low nonnegative rank impothesis. The algorithm alternates
% between the cluster affectation (choosing the k-closest subspaces for each
% points) with a nonnegative matrix factorization where the k-sparsity of
% coefficients B is fixed. This is reminescent of k-planes, but the major
% differences are:
% 1/ We allow data points to lie on multiple hyperplanes.
% 2/ We use the fact that r hyperplanes in our problem are defined by only r
% vectors in matrix factor D.
% 3/ We account for nonnegativity constraints, and replace the PCA step in
% k-planes by a full NMF.
%
% Inputs :
% - M : n x m data matrix
% - k : number of nonzero coefficients in columns of B
% - iter : maximal number of iterations (optional)
% - A0 : initial factor matrix A (random if not provided)
% - B0 : initial factor matrix B (random if not provided)
%
% Outputs :
% - A : estimated factor A
% - B : estimated k-sparse factor B
% - aff : clustering results for columns of M
% - err : reconstruction error (relative)
%
% Updates :
% 07/10/2018 J. E.Cohen
% Creation of the file

A = A0;
B = B0;
normM = norm(M,'fro');
err(1) = norm(M-A*B,'fro')^2/normM;  err(2) = inf;
tol = 10^-5; % minimal relative reconstruction error
t = 0; % loop counter
iter_nn = 50; % Number of HALS iterations for nnls
iter_nmf = 10; % Number of A-B alternance for the partials NMFs
diff_aff = 1; % anything > 0 is OK
locus = ones(size(A,2)-k,size(M,2));

% Main loop
while t<iter && abs((err(end)-err(end-1))/err(end-1))>tol ...
        && err(end) > 10^-15 && diff_aff > 0
        
        
t = t+1;
clc
fprintf('K-summits running...\n')
fprintf('Completion... ')
fprintf('%d%%  \n',round(100*t/iter))
fprintf('Error : %d  ',err(end))

% Affectation
[aff,distances] = affect(M,A);
locus_old = locus; % for affectation difference spotting
locus = aff(:,1:(size(A,2)-k))';
for pos =1:size(A,2) % TODO: improve implementation
    for pos2 = 1:(size(A,2)-k)
        B(pos,locus(pos2,:)==pos) = 0;
    end
end
diff_aff = sum( abs(locus_old(:) - locus(:)));

    % NMF
    for l=1:iter_nmf %only few iterations of alternating nnls in the inside loop
    
        % B estimate with nnls, but known zero positions
        B = nnlsHALSupdt_fixedsupp(M,A,B,iter_nn,locus);
        % TODO fix this function ?

        % A estimate using nonnegative least squares
        A = nnlsHALSupdt(M',B',A',iter_nn)';
        %A = A.*repmat(1./sqrt(sum(A.^2)),size(A,1),1);
        [A,B] = col_norm(A,B);

    end

% error computation
err(t+1) = norm(M-A*B,'fro')/normM;


end

if diff_aff == 0 % Redo NMF steps if zeros have been fixed
        % NMF
    for l=1:100 %do 20 iterations of alternating nnls at the end with fixed support
    %TODO: while condition on error
    
        % B estimate with nnls, but known zero positions
        B = nnlsHALSupdt_fixedsupp(M,A,B,iter_nn,locus);
        % TODO improve this function

        % A estimate using nonnegative least squares
        A = nnlsHALSupdt(M',B',A',iter_nn)';
        %A = A.*repmat(1./sqrt(sum(A.^2)),size(A,1),1);
        [A,B] = col_norm(A,B);

    end
end

% error computation
err(t+2) = norm(M-A*B,'fro')/normM;

clc
fprintf('K-summits running...\n')
fprintf('Completion... ')
fprintf('100%%  ')

function [aff,d] = affect(M,A)
d = zeros(size(M,2),size(A,2));
    for p=1:size(A,2)
        Aj = A(:,setdiff(1:size(A,2),p));
        
        U = orth(Aj);
        d(:,p) = sum((M - U*U'*M).^2);
        
        %U = Aj*((Aj'*Aj)\(Aj'));
        %d(:,p) = sum((M-U*M).^2);
        %d(:,p) = sum((M-Aj*(Aj\M)).^2);
    end
[~,aff] = sort(d,2);

end


end
