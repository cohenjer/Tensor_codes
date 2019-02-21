function [alloc,indices,X_hat] = MatchP_multiD(X,D,d,fun)
% Provide D as a structure, d is a row vector of same size as D containing
% the number of atoms to be selected from each library in D.
% alloc provides the dictionary index for each column of X, and indices
% the according atom in the dictionary.
% fun is a string, either 'eucl' (default), 'mrsa' or 'sam'. It fixes the
% distance metric used for the assignement problem.

if nargin<4
    fun = 'eucl';
end

[~,k] = size(d);
p = sum(d);

if k~=size(D,2)
    fprintf('\n Wrong number of dictionaries or wrong choice of d \n')
end
if sum(d<1)>0
    fprintf('\n WARNING: 0 in d will result in wrong outputs \n')
end
if p<size(X,2)
    %fprintf('\n WARNING: Different number of atoms to be matched and number of components')
    fprintf('\n WARNING: Under-estimation of the number of atoms to be seleted \n')
end

% Dictionnary normalization for MatchP (normalize all dicos similarly so that scores are comparable)
Dn = cell(1,size(D,2)); % initialization
for i=1:size(D,2)
normD     =     sqrt(sum(D{i}.^2));
Dn{i} = D{i}./(repmat(normD,size(D{i},1),1)+1e-6);
end

% Computing the cost
Gram = cell(1,k);
scores = zeros(size(X,2),k);
indices_1 = zeros(size(X,2),k);
for i=1:k
    % Gram contains the distances of all
    % atoms factors for each dictionary.
    %    Gram{i}     =    Dn{i}'*col_norm(X); 
    Gram{i}  =  feval(fun,X,Dn{i});
    if size(Gram{i},1)>1
        [scores(:,i),indices_1(:,i)] = max(Gram{i},[],2); 
    else
        scores(:,i) = Gram{i};
        indices_1(:,i) = ones(size(X,2),1);
    end
end

% The following lines account for multiplicity in d by duplicating the
% dictionaries in the scores.
mark = 1;
compteur = 0;
table = zeros(1,p);
table(1,1:k) = linspace(1,k,k); % original dictionary ref
for j=d
    for l=1:(j-1)
    scores(:,k+l+compteur) = scores(:,mark);
    indices_1(:,k+l+compteur) = indices_1(:,mark);
    %table = [table,mark];
    table(k+l+compteur) = mark;
    end
    compteur = compteur + d(mark) -1;
    mark = mark+1;
end
    

% Hungarian algorithm
[alloc,~] = munkres(-scores);
%[alloc,~] = munkres(scores); % if little is better

% Allocating the right atoms to each column of X
alloc = table(alloc);

% Storing the best atom for each X knowing alloc
%indices = zeros(1,p);
indices = zeros(1,size(X,2));
X_hat   = zeros(size(X));
for i=1:size(X,2)
indices(1,i) = indices_1(i,alloc(i));
% ------    % munkres version for no multiples
X_hat(:,i) = D{alloc(i)}(:,indices(i));
end
