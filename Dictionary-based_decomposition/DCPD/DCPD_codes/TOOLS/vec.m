% Implementation of the vec(x) operator
% obs: For a K x L x M tensor, the vec(x) operator
% stacks according to the kronecker product order

function a = vec(G)

[K,L,M]   =     size(G);
a   =     reshape(permute(G,[3,2,1]),K*L*M,1,1);

end