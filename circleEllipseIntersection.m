
function [ref1 ref2] = circleEllipseIntersection(cx, cy, r, ellipsefunc)
    
    psi1 = fzero(@(theta) ellipsefunc(r*cos(theta) + cx, r*sin(theta) + cy), 0.75*pi);
    psi2 = fzero(@(theta) ellipsefunc(r*cos(theta) + cx, r*sin(theta) + cy), -0.75*pi);
    
    ref1 = [r*cos(psi1) + cx, r*sin(psi1) + cy];
    ref2 = [r*cos(psi2) + cx, r*sin(psi2) + cy];
end