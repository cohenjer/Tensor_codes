%% Script for testing the proximal operator of the matrix induced l1 norm

% Dummy test
%toto = [ 3, 2, 1, 4.1 ; 2, 1, 0 , 3; 1, 0, 0, 2 ];
toto = randn(100,10);
tic
[U,t, nu]= prox_ml1(toto,10,10^-8);


%%
clear out

for n=1:10
    100*n
    for m=1:10
        temp = zeros(1,5);
        for i=1:5
            % Parameters
            lamb = 10;
            % Data
            X = randn(100*n,10*m);
            % Prox vs sum
            clock_alg=tic;
            prox_ml1(X,lamb,10^(-8));
            t=toc(clock_alg);
            % output
            temp(i) = t/(10*m*100*n*log2(100*n));
            clear clock_alg
            clear clock_ref
        end
        out(n,m) = mean(temp);
    end
end

surf(100:100:1000,10:10:100,log10(out)')
axis([100 1000 10 100])
xlabel('n')
ylabel('m')
zlabel('log_{10}(t/nm log_2(n))')
title('')
view(0,90)
c = colorbar;
c.Label.String = 'log_{10}(t/nm log_2(n))';
%title('Effective computation time of prox of l11')