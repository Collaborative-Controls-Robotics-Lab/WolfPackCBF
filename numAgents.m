
% Number of agents required to span a length
function [N] = numAgents(velRatio, epsilon, L)
    q = ((1/velRatio) + 1)/((1/velRatio) - 1);
    N = ceil(2*log(L/(2*epsilon/velRatio) + 1)/log(q));

%    if distSpanned(velRatio, epsilon, N) > L
    if mod(N, 2) ~= 0
        N =  ceil(2 *log(1 + (L - 2*epsilon)/(2*epsilon*(1/velRatio + 1)))/log(q)+ 1) ;
    end

end
