function [ d ] = classif_err( K, Ke )
%Returns the number of missmatch between selection matrices S (ground
%truth) and Se (estimated).


% problem with repeated values in Ke
% if length(K)>length(Ke)
%     Ke = [Ke,zeros(1,length(K)-length(Ke))];
%     d = length(setdiff(Ke,K));
% elseif length(K)<length(Ke)
%     K = [K,zeros(1,length(Ke)-length(K))];
%     d = length(setdiff(K,Ke));
% else
%     d = length(setdiff(K,Ke));
% end

% Finding joint support (no repetitions in K so it's ok) then taking
% opposite.
d = max(size(Ke,2),size(K,2)) - length(intersect(K,Ke));

end

