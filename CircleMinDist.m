function [point, distance] = CircleMinDist(x0, y0, r, xstar, ystar)
% Find closest point on circle to a fixed point
    % If x0 \neq xstar
    m = (y0 - ystar)/(x0 - xstar);
    y = @(x) m*(x - x0) + y0;
    circ = @(x) (x - x0)^2 + (y(x) - y0)^2 - r^2;
    xcrit = fzero(circ, x0 - r);
    ycrit = m*(xcrit - x0) + y0;
    point = [xcrit; ycrit];
    distance = sqrt((xstar - xcrit)^2 + (ystar - ycrit)^2);
end