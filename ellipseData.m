function [A,B,C,D,E,F] = ellipseData(focus1, focus2, r)
a = focus1(1);
b = focus1(2);
c = focus2(1);
d = focus2(2);
A = (a - c)^2 - r^2;
B = 2*(a - c)*(b - d);
C = (b - d)^2 - r^2;
D = r^2*(a + c) - (a - c)*(a^2 + b^2 - c^2 - d^2);
E = r^2*(b + d) - (b - d)*(a^2 + b^2 - c^2 - d^2);
F = 1/4*(r^4 - 2*r^2*(a^2 + b^2 + c^2 + d^2) + (a^2 + b^2 - c^2 - d^2)^2);
%fun = @(x,y) A*x.^2 + B*x.*y + C*y.^2 + D*x + E*y + F;
end