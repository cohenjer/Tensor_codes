function [U,t, nu_t] = prox_ml1(X,l,tol)
% U is the prox of max_i \|X\|_1, t the value of the max.
% Input X is a matrix, l the lambda parameter. tol is the desired
% precision, and the stopping criterion will be tmax*tol.


[n,m] = size(X);

% Precomputation
    % Summing, sorting
Xsort = sort(abs(X),'descend');
Xsum = sum(Xsort,1);
Xcumsum = cumsum(Xsort);
    % maximal value of t (no thresholding)
tmax = max(Xsum);
delta = tol;
tmin = 0; 
t = tmax/2;
%tol = 1e-8; % tolerance on estimation of sum(nu_i)
    % Order of visited columns of X
[sorted_sums, order] = sort(Xsum,'descend');
%t = sum( abs( ST( X(:,order(1) ),l) ) ); %first value
    % Tools for computations later on
ps = linspace(1,n,n)';

% Remove zero solutions (bugged otherwise)
if l>=sum(Xsort(1,:))
    U = zeros(n,m);
    t = 0;
    nu_t = Xsort(1,:)/l;
    return
end

% Bisection
while tmax - tmin > delta 
    %fprintf('looping...')
    % 1/ Computing the current active set
    I = order(sorted_sums>t);
    i = length(I); %todo improve
    
    % cummulative sum of x and preprocessing
    Xcumsum_m_t = (Xcumsum(:,order(1:i))-t)/l;
    
    % Compute the candidate nu values (warning: already sorted)
    Ps = repmat(ps, 1, i);
    Nu = Xcumsum_m_t./Ps;
    
    % Test to keep the only good one
    % Check which nu is compatible with the Soft Thresholding number of nonzero
    % coeffs p.
    nu_t = zeros(1,m);
    for j=1:i % O(nm) in worst case :(
        % TODO: There has to be a better way since Xsort is sorted !!
        for p=n:-1:1
            % does xsort(p) < \lambda\nu(p) < xsort(p+1)
            if p==n % This is often the result is early stages where all is clipped
                if Xsort(p,order(j)) >= l*Nu(p,j)
                    nu_t(order(j)) = Nu(p,j);
                    break
                end
            else
                if Xsort(p,order(j)) >= l*Nu(p,j) && l*Nu(p,j) > Xsort(p+1,order(j))
                    nu_t(order(j)) = Nu(p,j);
                    break
                end %end if
            end %end if p==n
        end %end for p
    end %end for j
    
    % 3/ Checking dual condition, 1< or 1>
    if sum(nu_t) < 1 
        % t must be decreased
        tmax = t;
        t = (t + tmin)/2;
    else
        % t must be increased
        tmin = t;
        t = (t + tmax)/2;
    end
end

% Final step, thresholding vectors that need to be
U = X;
for j=order(1:i)
    U(:,j) = ST(X(:,j), l*nu_t(j));
end

end
