function theta = sample_dirichlet(alpha, N)
% theta = sample_dirichlet(alpha, N)
%
% Code downladed from 
%
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/22703
%
% 
% SAMPLE_DIRICHLET Sample N vectors from Dir(alpha(1), ..., alpha(k))
% theta = sample_dirichlet(alpha, N)
% theta(i,j) = i'th sample of theta_j, where theta ~ Dir
% 
% We use the method from p. 482 of "Bayesian Data Analysis", Gelman et al.

k = length(alpha);
theta = zeros(N, k);
scale = 1; 
for i=1:k
  theta(:,i) = gamrnd(alpha(i), scale, N, 1);
end
S = sum(theta,2); 
theta = theta ./ repmat(S, 1, k);

end % of function sample_dirichlet