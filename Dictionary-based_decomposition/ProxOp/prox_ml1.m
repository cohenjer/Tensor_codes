function [U,t, nu_t] = prox_ml1(X,l,tol)
% U is the prox of max_i \|X\|_1, t the value of the max.
% Input X is a matrix, l the lambda parameter. tol is the desired precision.


[n,m] = size(X);

% Precomputation
    % Summing, sorting
Xsort = sort(abs(X),'descend');
Xsum = sum(Xsort,1);
Xcumsum = cumsum(Xsort);
    % maximal value of t (no thresholding)
tmax = max(Xsum);
tmin = 0;
%tol = 1e-8; % tolerance on estimation of sum(nu_i)
    % Order of visited columns of X
[Xsumsort,order] = sort(Xsum,'descend');
t = sum( abs( ST( X(:,order(1) ),l) ) ); %first value
    % Tools for computations later on
ps = linspace(1,n,n)';
Ps = repmat(ps, 1, m);

% Remove zero solutions (bugged otherwise)
if l>=sum(Xsort(1,:))
    U = zeros(n,m);
    t = 0;
    nu_t = Xsort(1,:)/l;
    return
end
% Loop on number of considered columns
for i=1:m
   % Selecting columns
   %Xs = X(:,order(1:i)); % Thresholded columns
   
   % Computing the value of max_{j<=i} \|U_j\|_1 
   % with U_j = S_{\lambda nu}(X_j)
   loop = true;
   while loop 
        %fprintf('looping...')
        % 1/ Computing the dual variables
               
            % cummulative sum of x and preprocessing
            Xcumsum_m_t = (Xcumsum-t)/l;

            % Compute the candidate nu values
            Nu = Xcumsum_m_t./Ps;
           
            % Test to keep the only good one
            % Check which nu is compatible with the Soft Thresholding number of nonzero
            % coeffs p.
            nu_t = zeros(1,m);
            for j=1:m % O(nm) in worst case :(
                % TODO: There has to be a better way since Xsort is sorted !!
                for p=n:-1:1
                    % does xsort(p) < \lambda\nu(p) < xsort(p+1)
                    if p==n % This is often the result is early stages where all is clipped
                        if Xsort(p,j) >= l*Nu(p,j)
                            nu_t(j) = Nu(p,j);
                            break
                        end
                    else 
                        if Xsort(p,j) >= l*Nu(p,j) && l*Nu(p,j) > Xsort(p+1,j)
                            nu_t(j) = Nu(p,j);
                            break
                        end %end if
                    end %end if p==n
                end %end for p
            end %end for j
       
        
        % 2/ Checking dual feasibility
        if abs(sum(nu_t)-1)<tol
           loop = false;
        end
        % 3/ Checking dual condition, 1< or 1>
        if sum(nu_t) < 1 && loop
            % t must be decreased
            tmax = t;
            t = (t + tmin)/2;
        elseif loop
            % t must be increased
            tmin = t;            
            t = (t + tmax)/2;
        end
   end
   
   % Checking primal feasibility
   if i==m || t>Xsumsort(i+1) 
       % The right solution has been found, now we only need to clip these
       % vectors
       break
   end
   
   % Resetting the tmax variable for next loop
   % Indeed tmin will only increase from the t value
   tmax = max(Xsum);
   tmin = t; 
   %tmin = 0;
   
end % loop on thresholded vectors

% Final step, thresholding vectors that need to be
U = X;
for j=order(1:i)
    U(:,j) = ST(X(:,j), l*nu_t(j));
end

end
