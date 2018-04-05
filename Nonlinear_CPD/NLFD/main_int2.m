function main_int2
% Generates the simulations for the article in CILS 2015-2016

cas=2;
beta=2;
tab=2:6;


switch cas
    case 2
        load cas_OK_2_compo.mat
        Lv{2}=Lv{2}(4:beta:end,:);
        
        T=T(:,4:beta:end,:);
        T=T(tab,:,:);
        Lv{1}=Lv{1}(tab,:);
        ex=[290:5*beta:500];
        em=[350:5:600];
        
        S=(em(1)-ex(1))/(5*beta)+1
        size(T)
        
        [Lf,res,r]=LM_IF_int(T,2,S,beta,100,'b',ex,em);%,Lv);
        
        
    case 3
        [TMc,Trefc,Tlin,Cref,EEM_M_lin,EEM_N_lin]=tenseur_meef;
        Cref([10,16],:)=[];
        Trefc([10,16],:,:)=[];
        Tec1=Trefc;
        Trefc(9,:,:)=Tec1(14,:,:);
        Trefc(14,:,:)=Tec1(9,:,:);
        Tec1=Cref;
        Cref(9,:)=Tec1(14,:);
        Cref(14,:)=Tec1(9,:);
        a=1:14;
        Te=Trefc(tab,:,:);
        Cref=Cref(tab,:);
        [a(tab)' Cref]
        % return
        ex=270:5*beta:550;
        em=270:5:550;
        
        Te=Te(:,1:beta:end,:);
        size(Te)
        S=(em(1)-ex(1))/(5*beta)+1
        
        [Lf,res,r,Lp]=LM_IF_int(Te,3,S,beta,100,'b',ex,em);
        respara2(Lf,Cref,1,'b',beta);
        respara2(Lp,Cref,0,'g',beta);
end

