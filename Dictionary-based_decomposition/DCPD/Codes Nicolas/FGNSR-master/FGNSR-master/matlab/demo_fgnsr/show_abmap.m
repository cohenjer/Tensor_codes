function [f_handle, Vaff] = show_abmap(V, lig, Li, Co, bw)
% f_handle = show_abmap(V,lig,Li,Co)
%
% Interpret given data as abundance data and display it as an image.
%
% Input:
%   V              : (m x r) matrix whose colums contains vectorized images
%  lig             : number of images per row in the display
%   (Co,Li)        : dimensions of images
%   bw             : if bw=1, reverse gray level
%
% Output:
%   f_handle -- Figure handle
%   Vaff     -- image raw data

V = max(V,0);
[~,r] = size(V);

for i = 1 : r
    if max(V(:,i)) > 0
        V(:,i) = V(:,i)/max(V(:,i));
    end
end

Vaff = []; 
for i = 1 : r
        ligne = floor((i-1)/lig)+1; col = i - (ligne-1)*lig;
        Vaff((ligne-1)*Li+1:ligne*Li,(col-1)*Co+1:col*Co) = reshape(V(:,i),Li,Co);
end

[m,n] = size(Vaff);
for i = 1 : n/Co-1
        Vaff = [Vaff(:,1:Co*i+i-1) 0.5*ones(m,1) Vaff(:,Co*i+i:end)];
end

[m,n] = size(Vaff);
for i = 1 : m/Li-1
        Vaff = [Vaff(1:Li*i+i-1,:); 0.5*ones(1,n); Vaff(Li*i+i:end,:)];
end

f_handle = figure;
if nargin == 5 && bw == 1
    imshow(Vaff,[0 1]); 
else
    imshow(1-Vaff,[0 1]); 
end
colormap(gray); 

end % of show_abmapx