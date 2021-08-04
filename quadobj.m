function [y,grady] = quadobj(x,Q,f)
    % Quadratic objective function for optimization
    y = 1/2*x.'*Q*x + f.'*x;
    if nargout > 1
        grady = Q*x + f;
    end
end