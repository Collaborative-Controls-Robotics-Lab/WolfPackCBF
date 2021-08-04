function hess = quadhess(x,lambda,Q,H)
% Quadratic hessian matrices for inequalities
hess = Q;
jj = length(H);
for i = 1:jj
    hess = hess + lambda.ineqnonlin(i)*H{i};
end
end