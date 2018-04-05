function [Gram] = mrsa(X,D)
% Mean-removed spectral angle distance function

Dm = mean(D,1);
Xm = mean(X,1);

Dr = D-Dm; Dr = col_norm(Dr);
Xr = X-Xm; Xr = col_norm(Xr);

Gram =100/pi*real(acos(Xr'*Dr));

end
