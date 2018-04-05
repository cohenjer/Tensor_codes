function [ G ] = MultProd( Y,U,V,W )
%This function calculates the coordinates of a given tensor in a new basis
%(U,V,W)*Y (unfolding matrices as in Kolda 2008)

[K,L,M]   =     size(Y);

[u1,~]    =     size(U); %U of size u1 x K
[v1,~]    =     size(V);
[w1,~]    =     size(W);


T1  =     reshape( Y , K , L*M );
G   =     reshape( U*T1 , u1 , L , M ); 

T2  =     reshape( permute( G , [2,1,3] ) , L , u1*M );
G   =     permute( reshape( V*T2 , v1 , u1 , M ), [2,1,3] );

T3  =     reshape( permute( G , [3,1,2] ) , M , u1*v1 );
G   =     permute( reshape( W*T3 , w1 , u1 , v1 ) , [2,3,1] );

end

