function [ d ] = classif_err_group( K, Ke, group_interv )
%Returns the number of missmatch between true K (ground
%truth) and the estimated indexes Ke (estimated). group_interv is the
%number of atoms in each group (all the same).

% Finding joint group support (no repetitions in K so it's ok) then taking
% opposite.
G  = floor((K-1)/group_interv)+1;
Ge = floor((Ke-1)/group_interv)+1;
g_joint_support = intersect(G,Ge);
d = max(size(Ke,2),size(K,2)) - length(g_joint_support);

end

