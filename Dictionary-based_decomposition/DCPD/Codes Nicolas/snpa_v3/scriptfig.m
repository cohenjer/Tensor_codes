colors{1} = 'b' ;   % SPA
colors{2} = 'r-.' ; % SNPA
colors{3} = 'mo-' ; % XRAY 
colors{4} = 'ko-' ; % NMFdico 

figure; 
for i = 1 : nalgo 
    semilogx(delta, results(i,:), colors{i}, 'LineWidth', 2); hold on;
end
xlabel('\delta'); ylabel('Fraction of columns of W exctracted'); 
legend('SPA', 'SNPA', 'XRAY', 'NMFdico'); 

figure; 
for i = 1 : nalgo 
    semilogx(delta, resultsree(i,:), colors{i}, 'LineWidth', 2); hold on;
end
xlabel('\delta'); ylabel('Fraction of columns of W exctracted'); 
legend('SPA', 'SNPA', 'XRAY', 'NMFdico'); 