function partial = dxstardxi(A,B,C,D,E,F,x,y,lambda,sigma)
% vector of unknowns = [dxdxix; dydxix; dlambdadxix];
% [dxdxiy; dydxiy; dlambdadxiy];
% First system of equations
csigma = 1/(1 - sigma^2);
A1 = zeros(3,3);
b1 = zeros(3,1);
A1(1,:) = [2 + 2*lambda*A; B*lambda; 2*A*x + B*y + D].';
b1(1) = 2*csigma;
A1(2,:) = [B*lambda; 2 + 2*C*lambda; 2*C*y + B*x + E].';
b1(2) = 0;
A1(3,:) = [2*A*x + B*y + D; B*x + 2*C*y + E; 0].';
b1(3) = 0;
first = A1\b1;

b2 = [0; 2*csigma; 0];
second = A1\b2;
partial = [first(1:2) second(1:2)];

end