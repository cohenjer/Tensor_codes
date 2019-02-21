function [Gram] = sam(X,D)
% spectral angle distance function

Gram = 100/pi*real(acos(col_norm(X)'*D));

end
