function partial = dxstardt(A,B,C,D,E,F,x,y,lambda,umax,p,xe,r)
xex = xe(1); xey = xe(2); px = p(1); py = p(2);

% vector of unknowns = [dxdxix; dydxix; dlambdadxix];
% [dxdxiy; dydxiy; dlambdadxiy];
% First system of equations
dAdt = 2*r*umax; dBdt = 0;
dCdt = 2*r*umax; dDdt = -2*r*umax*(xex + px);
dEdt = -2*r*umax*(xey + py); dFdt = -r^3*umax + r*umax*(xex^2 + xey^2 - px^2 - py^2);

A1 = zeros(3,3);
b1 = zeros(3,1);
A1(1,:) = [2 + 2*lambda*A; B*lambda; 2*A*x + B*y + D].';
A1(2,:) = [B*lambda; 2 + 2*C*lambda; 2*C*y + B*x + E].';
A1(3,:) = [2*A*x + B*y + D; B*x + 2*C*y + E; 0].';
b1(1) = -2*lambda*x*dAdt - lambda*y*dBdt - lambda*dDdt;
b1(2) = -2*lambda*y*dCdt - lambda*x*dBdt - lambda*dEdt;
b1(3) = -x^2*dAdt - x*y*dBdt - y^2*dCdt - x*dDdt - y*dEdt - dFdt;
first = A1\b1;

partial = first(1:2);

end