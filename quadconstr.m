function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
% Quadratic constraint functions
    jj = length(H);
    y = zeros(1,jj);
    
    for i = 1:jj
        y(i) = 1/2*x.'*H{i}*x + k{i}'*x + d{i};
    end
    
    yeq = [];

    if nargout > 2
        grady = zeros(length(x),jj);
        for i = 1:jj
            grady(:,i) = H{i}*x + k{i};
        end
    end
    
    gradyeq = [];
end