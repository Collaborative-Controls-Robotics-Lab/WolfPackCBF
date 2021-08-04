function [xs, ys] = CircleIntersection(x0, y0, r0, x1, y1, r1, guess)
% Intersection of two circles
A = (-2*x0 + 2*x1);
B = (-2*y0 + 2*y1);
C = r0^2 - r1^2 - x0^2 + x1^2 - y0^2 + y1^2;
y = @(x) C/B - A/B*x;
circ = @(x) (x - x0)^2 + (y(x) - y0)^2 - r0^2;
xs = fzero(circ, guess);
ys = y(xs);
end