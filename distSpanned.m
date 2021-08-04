

% Distance spanned by N agents at terminal time
function [L] = distSpanned(velRatio, epsilon, N)
    q = ((1/velRatio) + 1)/((1/velRatio) - 1);
    
    if mod(N, 2) == 1
        n = floor(N / 2);
        
        L = 2*epsilon - 2*epsilon*((1/velRatio) + 1)*(1 - q^n);
    else
        n = N / 2;
        
        L = 2*epsilon/velRatio * (q^n - 1);
    end

end