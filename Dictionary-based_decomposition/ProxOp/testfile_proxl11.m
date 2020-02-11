%% Script for testing the proximal operator of the matrix induced l1 norm

for n=1:100
    n
    for m=1:100
        temp = zeros(1,5);
        for i=1:5
            % Parameters
            lamb = 1;
            % Data
            X = randn(10*n,m);
            % Prox vs sum
            clock_alg=tic;
            prox_ml1(X,lamb,10^(-8));
            t=toc(clock_alg);
            % output
            temp(i) = t/(m^2*10*n);
            clear clock_alg
            clear clock_ref
        end
        out(n,m) = mean(temp);
    end
end

surf(10:10:1000,1:100,log10(out)')
xlabel('n')
ylabel('m')
zlabel('log_{10}(t/nm^2)')
title('')
view(0,90)
c = colorbar;
c.Label.String = 'log_{10}(t/nm^2)';
%title('Effective computation time of prox of l11')