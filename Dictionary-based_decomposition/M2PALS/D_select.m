function [ D,d,pos ] = D_select( DATA,k,indices )
% This function makes you manually select the portions of image h to be
% used as containing pure pixels. d is the number of expected pure pixels
% per region. DATA is the 3-way hyperspectral data (wavelengths 1st mode).

if nargin<3
    indices = 'none';
end

imdata = sum(DATA,1);
imdata = permute(imdata,[2,3,1]);
imdata = imdata(:,:,1);
figure
imagesc(imdata)
pos = [];

if ~strcmp(indices,'none')
    hold on
    plot(indices(1,:),indices(2,:),'+red')
end


for i=1:k
    fprintf('\n Selection of rectangle number (double click after selection) %d \n',i)
    a = imrect;
    %prompt= 'Are you satisfied with your selection (1==y or 0==n) ? ';
    %reponse = input(prompt);
    %if reponse ==0
    %    fprintf('Second selection (there will be no 3d trial ;) ) \n')
    %    a = getrect;
    %end
    a = round(wait(a));
    D{i} = DATA(:,a(2):(a(2)+a(4)),a(1):(a(1)+a(3))); %D{i} = DATA(:,a(1):(a(1)+a(3)),a(2):(a(2)+a(4)));
    D{i} = reshape(D{i},size(D{i},1),size(D{i},2)*size(D{i},3));
    prompt ='How many pure pixels in this rectangle ? (type desired integer) ';
    d(i) = input(prompt);
    pos  = [pos;a];
end

% ajouter aperçu final, demande de réajustement

close all
