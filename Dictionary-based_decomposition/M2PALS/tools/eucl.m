function [Gram] = eucl(X,D)
% Euclidean distance function

Gram = col_norm(X)'*D;

end
