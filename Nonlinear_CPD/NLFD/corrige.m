function [Ae_new,Be_new,Ce_new]=corrige(Ae,Be,Ce,A,B,C)
% Corrects the ambiguities if possible

N=size(A,2);

for i= 1:N
    An(:,i)=sign(A(:,i)).*(A(:,i))/norm(A(:,i));
    Aen(:,i)=sign(Ae(:,i)).*(Ae(:,i))/norm(Ae(:,i));
    Bn(:,i)=sign(B(:,i)).*(B(:,i))/norm(B(:,i));
    Ben(:,i)=sign(Be(:,i)).*(Be(:,i))/norm(Be(:,i));
    Cn(:,i)=sign(C(:,i)).*(C(:,i))/norm(C(:,i));
    Cen(:,i)=sign(Ce(:,i)).*(Ce(:,i))/norm(Ce(:,i));
end

F1=[Aen;Ben;Cen];
F2=[An;Bn;Cn];

for i= 1:N
    F1(:,i)=(F1(:,i))/norm(F1(:,i));
    F2(:,i)=(F2(:,i))/norm(F2(:,i));
end
%%% calcul du critere de mesure de distance
D = F1'*F2;
[l,c] = size(D);
D = ones(l,c) - abs(D);
[Crit ind] = min(D);
Ae_new=Ae(:,ind);
Be_new=Be(:,ind);
Ce_new=Ce(:,ind);
for i= 1:N
    Ae_new(:,i)=Ae_new(:,i).*sign(A(:,i)).*sign(Ae_new(:,i))*norm(A(:,i))/norm(Ae_new(:,i));
    Be_new(:,i)=Be_new(:,i).*sign(B(:,i)).*sign(Be_new(:,i))*norm(B(:,i))/norm(Be_new(:,i));
    Ce_new(:,i)=Ce_new(:,i).*sign(C(:,i)).*sign(Ce_new(:,i))*norm(C(:,i))/norm(Ce_new(:,i));
end

