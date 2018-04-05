
function T=tenseur(L)
% create tensor T of dimensions dim_T and rank rang_T from the matrices
% stored in the cell array L
ordre_T=length(L);
rang_T=size(L{1},2);
for i=1:ordre_T
    dim_T(i)=size(L{i},1);
end
T=zeros(dim_T);
for f=1:rang_T
    TT=ones(dim_T);
    for d=1:ordre_T
        v1=ones(1,ordre_T);
        v1(d)=dim_T(d);
        v2=dim_T;
        v2(d)=1;
        TT=repmat(reshape(L{d}(:,f),v1),v2).*TT;
    end
    T=TT+T;
end