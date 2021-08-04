function partial = dxstardxe(A,B,C,D,E,F,x,y,lambda,sigma,p,xe,r)
xex = xe(1); xey = xe(2); px = p(1); py = p(2);
csigma = 1/(1 - sigma^2);
% vector of unknowns = [dxdxix; dydxix; dlambdadxix];
% [dxdxiy; dydxiy; dlambdadxiy];
% First system of equations
dAdxex = 2*(xex - px); dAdxey = 0;
dBdxex = 2*(xey - py); dBdxey = 2*(xex - px);
dCdxex = 0; dCdxey = 2*(xey - py);
dDdxex = r^2 - (xex^2 + xey^2 - px^2 - py^2) - (xex - px)*(2*xex);
dDdxey = -2*(xex - px)*xey; dEdxex = -2*(xey - py)*xex;
dEdxey = r^2 - (xex^2 + xey^2 - px^2 - py^2) - (xey - py)*(2*xey);
dFdxex = -r^2*xex + (xex^2 + xey^2 - px^2 - py^2)*(xex);
dFdxey = -r^2*xey + (xex^2 + xey^2 - px^2 - py^2)*(xey);

A1 = zeros(3,3);
b1 = zeros(3,1);
A1(1,:) = [2 + 2*lambda*A; B*lambda; 2*A*x + B*y + D].';
A1(2,:) = [B*lambda; 2 + 2*C*lambda; 2*C*y + B*x + E].';
A1(3,:) = [2*A*x + B*y + D; B*x + 2*C*y + E; 0].';
b1(1) = -sigma^2*csigma - 2*lambda*x*dAdxex - lambda*y*dBdxex - lambda*dDdxex;
b1(2) = -2*lambda*y*dCdxex - lambda*x*dBdxex - lambda*dEdxex;
b1(3) = -x^2*dAdxex - x*y*dBdxex - y^2*dCdxex - x*dDdxex - y*dEdxex - dFdxex;
first = A1\b1;

b2 = zeros(3,1);
b2(1) = -2*lambda*x*dAdxey - lambda*y*dBdxey - lambda*dDdxey;
b2(2) = -sigma^2*csigma -2*lambda*y*dCdxey - lambda*x*dBdxey - lambda*dEdxey;
b2(3) = -x^2*dAdxey - x*y*dBdxey - y^2*dCdxey - x*dDdxey - y*dEdxey - dFdxey;
second = A1\b2;
partial = [first(1:2) second(1:2)];

end