function X = parproj(Z, w, ub, diag_idx, verbose)
% Projection on the Polytope W

if issparse(Z) && ~isa(Z, 'int64')
    diag_idx = int64(diag_idx);
elseif ~isa(Z, 'int32')
    diag_idx = int32(diag_idx);
end

X = parproj_priv(Z, w, ub, diag_idx - 1, int32(verbose));

end