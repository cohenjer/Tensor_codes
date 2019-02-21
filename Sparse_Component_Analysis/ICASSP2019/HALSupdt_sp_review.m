function V = HALSupdt_sp_review(V,lambda,UtU,UtM,itermax)
[r,~] = size(V); 
eps0 = 0;
eps = 1;
cnt = 1;
delta = 1e-16;
while cnt <= itermax && eps >= (delta)^2*eps0
    nodelta = 0; 
        for k = 1 : r
            Vkold = V(k,:);
            V(k,:) = max(((UtM(k,:)-UtU(k,:)*V)-lambda)/UtU(k,k)+V(k,:),0); % update verified
            %V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + norm(Vkold - V(k,:),'fro')^2; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta;
    end
    eps = nodelta; cnt = cnt+1; 
end
