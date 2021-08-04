function [x,y,lambda] = closestEllipse(point, focus1, focus2, r)
% Find closest point on the boundary of the ellipse using Lagrange
% multipliers
ox = point(1);
oy = point(2);
[A,B,C,D,E,F] = ellipseData(focus1,focus2,r);
x0 = [ox; oy; 0.1];
for ii = 1:20
    Fx = [2*(x0(1) - ox) + 2*x0(3)*A*x0(1) + x0(3)*B*x0(2) + x0(3)*D; ...
        2*(x0(2) - oy) + 2*x0(3)*C*x0(2) + x0(3)*B*x0(1) + x0(3)*E; ...
        A*x0(1)^2 + B*x0(1)*x0(2) + C*x0(2)^2 + D*x0(1) + E*x0(2) + F];
    J = [2 + 2*x0(3)*A, x0(3)*B, 2*A*x0(1) + B*x0(2) + D; x0(3)*B, 2 + 2*x0(3)*C, 2*C*x0(2) + B*x0(1) + E; ...
        2*A*x0(1) + B*x0(2) + D, 2*C*x0(2) + B*x0(1) + E, 0];
    x0 = x0 - J\Fx;
end
x = x0(1);
y = x0(2);
lambda = x0(3);
end