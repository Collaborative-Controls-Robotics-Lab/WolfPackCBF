function [circVal] = LargestSubCircle(circVals, x0)
    cX = circVals(1,:);
    cY = circVals(2,:);
    cR = circVals(3,:);

    m = length(radii);
    Q = zeros(3,3);
    f = [0; 0; -1];
    qo = @(x) quadobj(x, Q, f);

    H = cell(m);
    k = cell(m);
    d = cell(m);
    for ii = 1:m
        H{ii} = [2 0 0;
                 0 2 0;
                 0 0 -2];
        k{ii} = [-2*cX(ii);
                 -2*cY(ii);
                 2*cR(ii)];
        d{ii} = cX(ii)^2 + cY(ii)^2 - cR(ii)^2;
    end
    qc = @(x) quadconstr(x, H, k, d);

    qh = @(x, lambda) quadhess(x, lambda, Q, H);

    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',qh);
    [x,fval,eflag,output,lambda] = fmincon(qo,x0,[],[],[],[],[],[],qc,options);

    circVal = x;
end