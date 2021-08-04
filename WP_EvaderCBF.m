%% An attempt to reproduce WolfPACK Results by Alexander Davydov
% Modified by Mikhail Khrenov

% March 2021 - UMD CCRL
% Using unicycle dynamics

%%

vid = false;
rng(140)

% Number of evaders
M = 1;
% Integration time
dt = 0.05;

xpur = nan(T,2*N);
xeva = nan(T,2);


% Capture distance
epsilon = 0.02;
k = 1;
cvxflag = 0;

% Domain
Dxmin = -0.3;
Dxmax = 0.9;
Dymin = -0.5;
Dymax = 0.5;

% Target set
Pxmin = 0.525;
Pxmax = 0.575;
Pymin = -0.025;
Pymax = 0.025;
Pxavg = 0.5*(Pxmin + Pxmax);
Pyavg = 0.5*(Pymin + Pymax);

% Maximum control inputs
vmax = 0.01;
wmax = 2*pi;
closest = 100000;
velRatio = 0.5;

% Final time
tf = 60;
% Number of iterations
T = tf/dt;


% Set initial evader positions and headings
% xe = [-0.4*ones(1,M); .0*ones(1, M); zeros(1,M)];
xe = [-0.16;
   -0.0; 0];

% Terminal time
t_terminal = norm(xe(1:2, 1) - [Pxavg; Pyavg])/(vmax*(1 + 1/velRatio))
% Dist to span at terminal time
L = ((vmax/velRatio*(tf - t_terminal))^2 - (vmax*t_terminal)^2)/(vmax/velRatio)/(tf - t_terminal)

% Number of pursuers
N = numAgents(velRatio, epsilon, L)
L_pred = distSpanned(velRatio, epsilon, N)

% State Vectors
% Pursuers
xp = zeros(3,N);
adj = zeros(N, N);


% Pursuer-Evader Map
pem = zeros(N,M);
pem = [1, 0; 1, 0; 1, 0; 1, 0; 1, 0; 1, 0; 1, 0; 1, 0;];

% set initial pursuer position
%(-vmax*t_terminal + Pxavg)
xp = [Pxmin*ones(1,N); linspace(-0.25, 0.25, N); 0*ones(1,N)];
%xp = [0.4*ones(1,N) + 0.1*rand(1,N); linspace(-0.13, 0.15, N); pi*ones(1,N)];

% Virtual pursuers on the CRS boundary
xr = zeros(2, 2*M);
nOnes = ones(N,1);
adj = diag(nOnes(1:N-1), -1) + diag(nOnes(1:N-1), 1);
%xp = [0.4 0.4; 0.1 -0.3; pi pi]

% Initialize projections onto defense surface
p = zeros(1,N);
pCoord = zeros(2,N);
c = zeros(1,N);
cCoord = zeros(2,N);
pVoronoi = zeros(2,N);
sigmaC = nan(1,T);
h = nan(1,T);
umax = nan(1,T);
u = zeros(2*N,1);
unorm = zeros(N,1);
Aij = zeros(N-1, 2*N);
Aie = zeros(N+1, 2);
bie = zeros(N+1, 1);
bij = zeros(1, N-1);
hij = zeros(1, N-1);
flag = 0;
midpt = zeros(2, N-1);
He = cell(N+2,1);
ke = cell(N+2,1);
de = cell(N+2,1);

neigh = zeros(N, 2);

if vid == true
    %nFrames = 20;
    vidObj = VideoWriter('test', 'MPEG-4');
    vidObj.Quality = 75;
    vidObj.FrameRate = 50;
    open(vidObj);
end

% Plot everything
figure(1); hold on; axis equal

% Target set
patch([Pxmax Pxmax Pxmin Pxmin], [Pymax Pymin Pymin Pymax], 'green','HandleVisibility','off')
xlim([Dxmin Dxmax]);
ylim([Dymin Dymax]);

% Evader
[xs,ys] = getpatches(xe,0.05);
h1 = patch(xs, ys, 'red', 'DisplayName', 'Evader');

% Pursuers
[xs,ys] = getpatches(xp,0.05);
h2 = patch(xs,ys, 'blue','DisplayName', 'Pursuers');
CircValues = zeros(3, N+1);
IntersectionPoints = zeros(2, 2);
minPts = zeros(2,N);
minDst = zeros(1,N);
CRS = cell(2,1);


for jj = 1:M
    CRS{jj} = ellipse(xe(1:2,jj), [Pxavg; Pyavg], T*dt*vmax/velRatio);

    h3(jj) = fimplicit(CRS{jj}, 'r', 'LineWidth', 1.5,'DisplayName','Reachable Set');
%     h4 = plot(0,0, 'm', 'LineWidth', 1.5);
    %h6 = plot(0,0, 'gd', 'LineWidth', 2.5);
    %h7 = plot(0,0, 'b*', 'LineWidth', 2.5);
end

h5 = plot(xr(1,:),xr(2,:), 'md', 'LineWidth', 2.5,'DisplayName','Reference Point');

Apollonius = cell(N,1);
EpsilonCirc = cell(N, 1);
DefenseSurface = cell(N,1);

ngon_N = 8;
ngon_theta = linspace(0, 2*pi, ngon_N + 1).';
ngon_theta(end) = [];

ngon_H = [cos(ngon_theta), sin(ngon_theta)];
ngon_k = v_max * cos(pi/ngon_N) * ones(ngon_N, 1);

for jj = 1:M
    for ii = 1:N
        rs = norm([xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
        os = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1,jj); xe(2,jj)])/(1 - velRatio^2);
        ApCircle = circle(os(1), os(2), rs);
        CircValues(:,ii) = [os(1); os(2); rs];
        Apollonius{ii,jj} = fimplicit(ApCircle, 'b', 'LineWidth', 1.5,'HandleVisibility','off');
        % Find point on circles for minimum distance
        
        % plot epsilon circle
        epcirc = circle(xp(1,ii), xp(2,ii), epsilon);
        EpsilonCirc{ii} = fimplicit(epcirc, 'g', 'LineWidth', 1,'HandleVisibility','off');
    end
end
count = 1;

% Evader distance at terminal time
r_terminal = vmax * t_terminal;
nom_circle = circle(Pxavg, Pyavg, r_terminal);
fimplicit(nom_circle, 'g', 'LineWidth', 1.5,'HandleVisibility','off');

% Algorithm loop
for qq = 1:T
    u = zeros(2*N,1);
    for jj = 1:M
        CRS{jj} = ellipse(xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
        set(h3(jj), 'Function', CRS{jj});
    end
    
    % Construct nominal control points
%     [r2, r1] = circleEllipseIntersection(Pxavg, Pyavg, r_terminal, CRS{1});
    slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2, jj));
    poin = (-r_terminal * ([Pxavg; Pyavg] - xe(1:2, jj))/norm([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
    r1 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, 1);
    r2 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, N);
    
    poin = (-r_terminal * ([Pxavg; Pyavg] - xe(1:2, jj))/norm([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
    if norm([Pxavg; Pyavg] - xe(1:2, jj)) > norm([Pxavg; Pyavg] - poin.')
        r1 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, 1);
        r2 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, N);
    else
        poin = (-0.9 * ([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
        r1 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, 1);
        r2 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, N);
    end
    
    % Construct Apollonius Circles for each pursuer/evader pair
    for ii = 1:N
        for jj = 1:M
            CircValues(3,ii) = norm([xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
            CircValues(1:2,ii) = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1,jj); xe(2,jj)])/(1 - velRatio^2);
            ApCircle = circle(CircValues(1,ii), CircValues(2,ii), CircValues(3,ii));
            
            set(Apollonius{ii,jj}, 'Function', ApCircle);
            
            epcirc = circle(xp(1,ii), xp(2,ii), epsilon);
            set(EpsilonCirc{ii}, 'Function', epcirc);
            
            % Find point on circles for minimum distance
%             [pt, minD] = CircleMinDist(CircValues(1,ii), CircValues(2,ii), CircValues(3,ii), xe(1,jj), xe(2,jj));
%             minPts(:,ii) = pt;
%             minDst(ii) = minD;
        end
    end
    tic
    
    % CBF QCQP for evader
    for ii = 1:N
        for jj = 1:M
            Aie(ii,:) = -(xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii));
            He{ii} = zeros(2);
            ke{ii} = Aie(ii,:).';
            bie(ii) = -vmax + (norm(xe(1:2,jj) - xp(1:2,ii)) - 1.05*epsilon);
            de{ii} = -bie(ii);
        end
    end
    
    uvmax = vmax/velRatio;
    He{N+1} = zeros(2);
    Aie(N+1,:) = (xe(1:2, 1) - [Pxavg; Pyavg]).'/norm(xe(1:2, 1) - [Pxavg; Pyavg]);
    ke{N+1} = Aie(N+1,:).';
    bie(N+1) = -uvmax + 1e2*(tf - (qq)*dt)*uvmax - norm(xe(1:2, 1) - [Pxavg; Pyavg])^3;
    de{N+1} = -bie(N+1);
    He{N+2} = 2*eye(2);
    ke{N+2} = zeros(2,1);
    de{N+2} = -uvmax^2;
    
    
    % Solve QCQP
    % Compute alpha 1 angle for each pair
    angles = zeros(N+1,1);
    for ii = 1:N+1
        for jj = 1:M
            if ii == 1
                [x,y,~] = closestEllipse(CircValues(1:2,ii), xe(1:2, jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                xstar = [x; y];
                xei = norm(xe(1:2, jj) - xp(1:2,ii));
                xestar = norm(xe(1:2, jj) - xstar);
                xistar = norm(xp(1:2,ii) - xstar);
                % Angle via law of cosines
                angles(ii) = acos((xistar^2 - xestar^2 - xei^2)/(-2*xei*xestar));
            elseif ii == N+1
                [x,y,~] = closestEllipse(CircValues(1:2,ii-1), xe(1:2, jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                xstar = [x; y];
                xei = norm(xe(1:2, jj) - xp(1:2,N));
                xestar = norm(xe(1:2, jj) - xstar);
                xistar = norm(xp(1:2,N) - xstar);
                % Angle via law of cosines
                angles(ii) = acos((xistar^2 - xestar^2 - xei^2)/(-2*xei*xestar));
            else
                xei = norm(xe(1:2, jj) - xp(1:2,ii));
                xej = norm(xe(1:2, jj) - xp(1:2,ii-1));
                xij = norm(xp(1:2,ii) - xp(1:2,ii-1));
                angles(ii) = acos((xij^2 - xej^2 - xei^2)/(-2*xei*xej));
            end
        end
    end
    
    idx = find(angles == max(max(angles)));
    
    if idx == 1
        [x,y,lambda] = closestEllipse(CircValues(1:2,idx), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
        xstar = [x; y];
        [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
        num = xp(1:2,idx) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
        xie = xp(1:2,idx) - xe(1:2,jj);
        unom = velRatio*xie/norm(xie) + ((velRatio^2 - 1)*dxe - velRatio^2*eye(2)).'*num/norm(num);
        unom = unom/norm(unom)*uvmax;
    elseif idx == N+1
        [x,y,lambda] = closestEllipse(CircValues(1:2,idx-1), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
        xstar = [x; y];
        [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
        num = xp(1:2,N) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
        xie = xp(1:2,N) - xe(1:2,jj);
        unom = velRatio*xie/norm(xie) + ((velRatio^2 - 1)*dxe - velRatio^2*eye(2)).'*num/norm(num);
        unom = unom/norm(unom)*uvmax;
    else
        xie = xp(1:2,idx-1) - xe(1:2,jj);
        xje = xp(1:2,idx) - xe(1:2,jj);
        unom = xie/norm(xie) + xje/norm(xje);
        unom = unom/norm(unom)*uvmax;
    end
    % Go straight to the goal
    %unom = ([Pxavg; Pyavg] - xe(1:2,jj));
    % Separate agents 1 and 2
    %unom = velRatio*((xp(1:2,1) - xe(1:2,jj))/norm((xp(1:2,1) - xe(1:2,jj))) + (xp(1:2,2) - xe(1:2,jj))/norm(xp(1:2,2) - xe(1:2,jj)));
    % Separate agents 2 and 3
    %unom = ((xp(1:2,2) - xe(1:2,jj))/norm((xp(1:2,2) - xe(1:2,jj))) + (xp(1:2,3) - xe(1:2,jj))/norm(xp(1:2,3) - xe(1:2,jj)));
    %unom = unom/norm(unom)*uvmax;

    for jj = 1:M
        
        if (tf - (qq)*dt)*uvmax <= 1.08*norm(xe(1:2,jj) - [Pxavg; Pyavg])
            ue(:,jj) = ([Pxavg; Pyavg] - xe(1:2,jj))/norm(([Pxavg; Pyavg] - xe(1:2,jj)))*uvmax;
            flag = 1;
        else
%             cvx_begin quiet
%                variable ueop(2)
%                minimize( (ueop - unom).' * (ueop - unom))
%                subject to
%                    Aie * ueop <= bie;
%                    ueop.'*ueop <= uvmax^2;
%             cvx_end
%             
%             ue(:,jj) = ueop;
            nonlconstr = @(x)quadconstr(x,He,ke,de);
            ue(:,jj) = ([Pxavg; Pyavg] - xe(1:2,jj))/norm(([Pxavg; Pyavg] - xe(1:2,jj)))*uvmax;

            if any(nonlconstr(ue) > 0)
                %Fmincon options
                options = optimoptions(@fmincon, 'Algorithm', 'interior-point','Display', 'off',...
                        'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                        'HessianFcn', @(x,lambda)quadhess(x,lambda,eye(2),He), ...
                                    'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);
                fun = @(x)quadobj(x,zeros(2),xe(1:2,jj) - [Pxavg; Pyavg]);
                [ue(:,jj),fval,eflag,output,lambda] = fmincon(fun,[0;0],...
                    [],[],[],[],[],[],nonlconstr,options);
            end
            
        end
    end
    
    for ii = 1:N
        for jj = 1:M
            
            if pem(ii, jj) > 0.5
            
                % Find two nearest pursuers
                if qq == 1
                    f_dist = arrayfun(@(r_idx) norm(xp(1:2,r_idx) - xp(1:2, ii)), 1:N);
                    for iii = 1:N
                        if pem(iii, jj) < 0.5
                            f_dist(iii) = inf;
                        end
                    end
                    [val, indices] = sort(f_dist);
                    if1 = indices(2);
                    if2 = indices(3);
                    f1 = xp(1:2, indices(2));
                    f2 = xp(1:2, indices(3));
                    neigh(ii, 1) = if1;
                    neigh(ii, 2) = if2;
                else
                    if1 = neigh(ii, 1);
                    if2 = neigh(ii, 2);
                    f1 = xp(1:2, if1);
                    f2 = xp(1:2, if2);
                end

                if N == 1
                    slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2, jj));
                    % Find both intersection points
                    [A,B,C,D,E,F] = ellipseData(xe(1:2, jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    p = xp(1:2,ii);
                    if slope(1) == 0 % Vertical line x = const
                        x = p(1);
                        xr(1,(jj*2-1):2*jj) = [x x];
                        a = C;
                        b = B*x + E;
                        c = A*x^2 + D*x + F;
                        % Quadratic formula
                        xr(2,(jj*2-1):2*jj) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
                    else
                        m = slope(2)/slope(1);
                        a = A + m*B + m^2*C;
                        b = -m*B*p(1) + B*p(2) - 2*m^2*C*p(1) + 2*m*C*p(2) + D + m*E;
                        c = m^2*C*p(1)^2 - 2*m*C*p(1)*p(2) + C*p(2)^2 - m*E*p(1) + E*p(2) + F;
                        % Quadratic formula
                        xr(1,(jj*2-1):2*jj) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
                        xr(2,(jj*2-1):2*jj) = m*xr(1,(jj*2-1):2*jj) - m*p(1) + p(2);
                    end
                    ui = xe(1:2, jj) + xr(:,(jj*2-1)) + xr(:,2*jj) - 3*xp(1:2,ii);
                elseif ii == 1
                    %[x,y] = closestEllipse(xp(1:2,ii), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    %[x,y] = closestEllipse(CircValues(1:2,ii+1), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    %xr(:,1) = [x;y];
%                     slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2, jj));

                    xr(:,(jj*2)) = r1;
                    %ui = xr(:,1) + xp(1:2,ii+1) + xe(1:2,jj) - 3*xp(1:2,ii);
                    %ui = xr(:,1) + xp(1:2,ii+1) - 2*xp(1:2,ii);
                    %ui = (1 - velRatio^2)*(xr(:,1) + xp(1:2,ii+1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3));sin(xe(3))];
                    w1 = 1/(CircValues(3,ii) + 0);
                    w2 = 1/(CircValues(3,if1) + 0);
                    a1 = (CircValues(3,ii) + 0);
                    a2 = (CircValues(3,if1) + 0);
                    %ui = 0.5*xr(:,1) + 0.5*w1/(w1 + w2)*xp(1:2,ii) ...
                    %    + 0.5*w2/(w1 + w2)*xp(1:2,ii+1) - xp(1:2,ii);
                    ui = 0.5*xr(:,(jj*2)) + 0.5*w1/(w1 + w2)*CircValues(1:2,ii) ...
                        + 0.5*w2/(w1 + w2)*CircValues(1:2,if1) - CircValues(1:2,ii);
                    ui = (1 - velRatio^2)*ui + velRatio^2*ue(:,jj);
                    Aij = -(velRatio*(xp(1:2,ii) - xe(1:2, jj)).'/norm(xp(1:2,ii) - xe(1:2, jj)) ...
                        - (xp(1:2,ii) - f1).'/norm(xp(1:2,ii) - f1));
                    bij = a1/(a1 + a2)*(velRatio*((xe(1:2, jj) - xp(1:2,ii)).'/norm(xe(1:2, jj) - xp(1:2,ii)) + (xe(1:2, jj) - f1).'/norm(xe(1:2, jj) - f1)) ...
                        * ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2, jj)) + norm(f1 - xe(1:2, jj))) - norm(xp(1:2,ii) - f1) + epsilon)^3);
                    % A1 - CRS constraint
                    [x,y,lambda] = closestEllipse(CircValues(1:2,ii), xe(1:2, jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    [A,B,C,D,E,F] = ellipseData(xe(1:2, jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    xstar = [x;y];
                    dxi = dxstardxi(A,B,C,D,E,F,x,y,lambda,velRatio);
                    dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
                    dxt = dxstardt(A,B,C,D,E,F,x,y,lambda,vmax/velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
                    num = xp(1:2,ii) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
                    xie = xp(1:2,ii) - xe(1:2,jj);
                    dhdxi = velRatio*xie.'/norm(xie) - num.'/norm(num)*(eye(2) + (velRatio^2 - 1)*dxi);
                    dhdxe = -velRatio*xie.'/norm(xie) - num.'/norm(num)*((velRatio^2 - 1)*dxe - velRatio^2*eye(2));
                    dhdxstar = (1 - velRatio^2)*num.'/norm(num);
                    A1 = -dhdxi - dhdxstar*dxi;
                    b1 = dhdxe*ue(:,jj) + dhdxstar*(dxe*ue(:,jj) + dxt) + 1e3*(velRatio*norm(xie) - norm(num) + epsilon)^3;
                    A = [Aij; A1]; 
                    b = [bij; b1];
                    %z = quadprog(eye(2), -ui, A, b, [], [], ones(2,1)*-vmax/sqrt(2), ones(2,1)*vmax/sqrt(2));
                    %cvx_begin quiet
                    %    variable z(2)
                    %    minimize( (ui - z).' * (ui - z))
                    %   subject to
                    %        Aij * z <= bij
                    %        z.'*z <= vmax^2
                    %cvx_end
                    %tic


                    H = {zeros(2), zeros(2), 2*eye(2)};
                    k = {Aij', A1', zeros(2,1)};
                    d = {-bij, -b1, -vmax^2};
                    
                    nonlconstr = @(x)quadconstr(x,H,k,d);
                    
                    if norm(ui) > vmax || any(nonlconstr(ui) > 0)
                        %Fmincon options
                        options = optimoptions(@fmincon, 'Algorithm', 'interior-point','Display', 'off',...
                            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                            'HessianFcn', @(x,lambda)quadhess(x,lambda,eye(2),H), ...
                            'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);
                        fun = @(x)quadobj(x,eye(2),-ui);
                        [z,fval,eflag,output,lambda] = fmincon(fun,[0; 0],...
                            [],[],[],[],[],[],nonlconstr,options);
                    else
                        z = ui;
                    end

                    if norm(z) > vmax && cvxflag == 0
                        disp('using cvx')
                        cvx_begin quiet
                            variable z(2)
                            minimize( (ui - z).' * (ui - z))
                            subject to
                                Aij * z <= bij
                                norm(z) <= vmax
                        cvx_end
                    end
                    if (size(z,1) == 0 || any(isnan(z)))
                        %z = (xe(1:2,jj) - xp(1:2,ii))/norm(xe(1:2,jj) - xp(1:2,ii))*vmax;
                        z = ui/norm(ui)*vmax;
                        cvxflag = 1;
                    end
                    
                    
                    
                    u(1:2) = u(1:2) + pem(ii,jj) * z;
                    %toc
                    hij(1) = velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f1 - xe(1:2,jj))) - norm(xp(1:2,ii) - f1) + 2*epsilon;
                    midpt(:,ii) = w1/(w1 + w2)*xp(1:2,ii) + w2/(w1 + w2)*f1;
                elseif ii == N
                    %[x,y] = closestEllipse(xp(1:2,ii), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    %[x,y] = closestEllipse(CircValues(1:2,ii+1), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    %xr(:,2) = [x;y];
                    slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2,jj));

                    xr(:,2*jj-1) = r2;
                    %ui = xr(:,2) + xp(1:2, ii-1) + xe(1:2,jj) - 3*xp(1:2,ii);
                    %ui = xr(:,2) + xp(1:2, ii-1) - 2*xp(1:2,ii);
                    %ui = (1 - velRatio^2)*(xr(:,2) + xp(1:2,ii-1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3)); sin(xe(3))];
                    w1 = 1/(CircValues(3,if1) + 0);
                    w2 = 1/(CircValues(3,ii) + 0);
                    a1 = (CircValues(3,if1) + 0);
                    a2 = (CircValues(3,ii) + 0);
                    %ui = 0.5*xr(:,2) + 0.5*w2/(w2 + w1)*xp(1:2,ii) ...
                    %    + 0.5*w1/(w2 + w1)*xp(1:2,ii-1) - xp(1:2,ii);
                    ui = 0.5*xr(:,2*jj-1) + 0.5*w2/(w2 + w1)*CircValues(1:2,ii) ...
                        + 0.5*w1/(w2 + w1)*CircValues(1:2,if1) - CircValues(1:2,ii);
                    ui = (1 - velRatio^2)*ui + velRatio^2*ue(:,jj);
                    
                    Aij = -(velRatio*(xp(1:2,ii) - xe(1:2,jj)).'/norm(xp(1:2,ii) - xe(1:2,jj)) ...
                        - (xp(1:2,ii) - f1).'/norm(xp(1:2,ii) - f1));
                    bij = a1/(a1 + a2)*(velRatio*((xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii)) + (xe(1:2,jj) - f1).'/norm(xe(1:2,jj) - f1)) ...
                        *ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f1 - xe(1:2,jj))) - norm(xp(1:2,ii) - f1) + epsilon)^3);
                    [x,y,lambda] = closestEllipse(CircValues(1:2,ii), xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
                    xstar = [x;y];
                    dxi = dxstardxi(A,B,C,D,E,F,x,y,lambda,velRatio);
                    dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
                    dxt = dxstardt(A,B,C,D,E,F,x,y,lambda,vmax/velRatio,[Pxavg;Pyavg],xe(1:2,jj),(T-qq)*dt*vmax/velRatio);
                    num = xp(1:2,ii) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
                    xie = xp(1:2,ii) - xe(1:2,jj);
                    dhdxi = velRatio*xie.'/norm(xie) - num.'/norm(num)*(eye(2) + (velRatio^2 - 1)*dxi);
                    dhdxe = -velRatio*xie.'/norm(xie) - num.'/norm(num)*((velRatio^2 - 1)*dxe - velRatio^2*eye(2));
                    dhdxstar = (1 - velRatio^2)*num.'/norm(num);
                    A1 = -dhdxi - dhdxstar*dxi;
                    b1 = dhdxe*ue(:,jj) + dhdxstar*(dxe*ue(:,jj) + dxt) + 1e3*(velRatio*norm(xie) - norm(num) + epsilon)^3;
                    A = [Aij; A1]; 
                    b = [bij; b1];
                    %z = quadprog(eye(2), -ui, A, b, [], [], ones(2,1)*-vmax/sqrt(2), ones(2,1)*vmax/sqrt(2));
                    %cvx_begin quiet
                    %    variable z(2)
                    %    minimize( (ui - z).' * (ui - z))
                    %    subject to
                    %        Aij * z <= bij
                    %        z.'*z <= vmax^2
                    %cvx_end
                    H = {zeros(2), zeros(2), 2*eye(2)};
                    k = {Aij', A1', zeros(2,1)};
                    d = {-bij, -b1, -vmax^2};
                    
                    nonlconstr = @(x)quadconstr(x,H,k,d);
                    
                    if norm(ui) > vmax  || any(nonlconstr(ui) > 0)
                        %Fmincon options
                        options = optimoptions(@fmincon, 'Algorithm', 'interior-point','Display', 'off',...
                            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                            'HessianFcn', @(x,lambda)quadhess(x,lambda,eye(2),H), ...
                            'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);
                        
                        fun = @(x)quadobj(x,eye(2),-ui);
                        [z,fval,eflag,output,lambda] = fmincon(fun,[0; 0],...
                            [],[],[],[],[],[],nonlconstr,options);
                    else
                        z = ui;
                    end
                    
                    if norm(z) > vmax && cvxflag == 0
                        disp('using cvx')
                        cvx_begin quiet
                            variable z(2)
                            minimize( (ui - z).' * (ui - z))
                            subject to
                                Aij * z <= bij
                                norm(z) <= vmax
                        cvx_end
                    end
                    if (size(z,1) == 0 || any(isnan(z)))
                        %z = (xe(1:2,jj) - xp(1:2,ii))/norm(xe(1:2,jj) - xp(1:2,ii))*vmax;
                        z = ui/norm(ui)*vmax;
                        cvxflag = 1;
                    end
                    
                    
                    u(2*ii-1:2*ii) = u(2*ii-1:2*ii) + pem(ii,jj) * z;
                else
                    %ui = xp(1:2,ii - 1) + xp(1:2,ii + 1) + xe(1:2,jj) - 3*xp(1:2,ii);
                    Aij = zeros(2);
                    bij = zeros(2,1);
                    w0 = 1/(CircValues(3,if1) + 0);
                    w1 = 1/(CircValues(3,ii) + 0);
                    w2 = 1/(CircValues(3,if2) + 0);
                    a0 = (CircValues(3,if1) + 0);
                    a1 = (CircValues(3,ii) + 0);
                    a2 = (CircValues(3,if2) + 0);

                    %ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*xp(1:2,ii) ...
                    %    + 0.5*w0/(w0 + w1)*xp(1:2,ii-1) ...
                    %    + 0.5*w2/(w2 + w1)*xp(1:2,ii+1) - xp(1:2,ii);

                    ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*CircValues(1:2,ii) ...
                        + 0.5*w0/(w0 + w1)*CircValues(1:2,if1) ...
                        + 0.5*w2/(w2 + w1)*CircValues(1:2,if2) - CircValues(1:2,ii);%+ xe(1:2,jj) - 2*CircValues(1:2,ii+1);
                    ui = (1 - velRatio^2)*ui + velRatio^2*ue(:,jj);% + xe(1:2,jj) - xp(1:2,ii);
                    Aij(1,:) = -(velRatio*(xp(1:2,ii) - xe(1:2,jj)).'/norm(xp(1:2,ii) - xe(1:2,jj)) ...
                        - (xp(1:2,ii) - f2).'/norm(xp(1:2,ii) - f2));
                    Aij(2,:) = -(velRatio*(xp(1:2,ii) - xe(1:2,jj)).'/norm(xp(1:2,ii) - xe(1:2,jj)) ...
                        - (xp(1:2,ii) - f1).'/norm(xp(1:2,ii) - f1));
                    bij(1) = a1/(a1 + a2)*(velRatio*((xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii)) + (xe(1:2,jj) - f2).'/norm(xe(1:2,jj) - f2)) ...
                        *ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f2 - xe(1:2,jj))) - norm(xp(1:2,ii) - f2) + epsilon)^3);
                    bij(2) = a1/(a1 + a0)*(velRatio*((xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii)) + (xe(1:2,jj) - f1).'/norm(xe(1:2,jj) - f1)) ...
                        *ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f1 - xe(1:2,jj))) - norm(xp(1:2,ii) - f1) + epsilon)^3);
                    

                    H = {zeros(2), zeros(2), 2*eye(2)};
                    k = {Aij(1,:).', Aij(2,:).', zeros(2,1)};
                    d = {-bij(1), -bij(2), -vmax^2};
                    
                    nonlconstr = @(x)quadconstr(x,H,k,d);
                    
                    if norm(ui) > vmax || any(nonlconstr(ui) > 0)
                        %Fmincon options
                        options = optimoptions(@fmincon, 'Algorithm', 'interior-point','Display', 'off',...
                            'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                            'HessianFcn', @(x,lambda)quadhess(x,lambda,eye(2),H), ...
                                    'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);
                        fun = @(x)quadobj(x,eye(2),-ui);
                        [z,fval,eflag,output,lambda] = fmincon(fun,[0; 0],...
                            [],[],[],[],[],[],nonlconstr,options);
                    else
                        z = ui;
                    end
                    
                    
                    if norm(z) > vmax && cvxflag == 0
                        disp('using cvx')
                        cvx_begin quiet
                            variable z(2)
                            minimize( (ui - z).' * (ui - z))
                            subject to
                                Aij * z <= bij
                                norm(z) <= vmax
                        cvx_end
                    end
                    if (size(z,1) == 0 || any(isnan(z)))
                        %z = (xe(1:2,jj) - xp(1:2,ii))/norm(xe(1:2,jj) - xp(1:2,ii))*vmax;
                        z = ui/norm(ui)*vmax;
                        cvxflag = 1;
                    end
                    
                    
                    u(2*ii-1:2*ii) = u(2*ii-1:2*ii) + pem(ii,jj) * z;
                    hij(ii) = velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f2 - xe(1:2,jj))) - norm(xp(1:2,ii) - f2)+ 2*epsilon;
                    midpt(:,ii) = w1/(w1 + w2)*xp(1:2,ii) + w2/(w1 + w2)*f2;
                end

            %ui = (ui)/norm(ui)*vmax;
            %u(2*ii - 1:2*ii) = ui;
            %xp(1:2,ii) = xp(1:2,ii) + (ui) * dt;
            end
        end
    end
    toc
    for ii = 1:N
        xp(1:2,ii) = xp(1:2,ii) + u(2*ii-1:2*ii)*dt;
        xp(3,ii) = atan2(u(2*ii), u(2*ii-1));
        xpur(qq, 2*ii-1:2*ii) = xp(1:2,ii).';
        %xp(1:2,ii) = xp(1:2,ii) + vmax*(xe(1:2,jj) - xp(1:2,ii))/norm(xe(1:2,jj) - xp(1:2,ii))*dt;
        %xp(3,ii) = atan2(xe(2) - xp(2,ii), xe(1) - xp(1,ii));
    end
    h(qq) = min(hij);
    for ii = 1:N
        unorm(ii) = norm(u(2*ii-1:2*ii));
    end
    umax(qq) = min(unorm);
    
    set(h5, 'XData', xr(1,:), 'YData', xr(2,:));
    
    [xps,yps] = getpatches(xp, 0.05);
    for jj = 1:M
        xeva(qq,:) = xe(1:2,jj).';
        xe(1:2,jj) = xe(1:2,jj) + ue(:,jj)*dt;
        xe(3,jj) = atan2(ue(2,jj), ue(1,jj));
    end
    
    [xes,yes] = getpatches(xe, 0.05);
    set(h1, 'XData', xes, 'YData', yes);
    set(h2, 'XData', xps, 'YData', yps);
    drawnow limitrate
    
    last_closest = closest;
    dx = xp(1:2,1:N) - xe(1:2,1); 
    closest = 100000;
    
    %pause(0.05)
    for jj = 1:M
        for ii = 1:N
            ev_dist = norm(dx(1:2,ii));
            if ev_dist < closest
                closest = ev_dist;
            end
            
            if norm([xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)]) <= epsilon
                flag = 1;
            end
        end
        
        if xe(1,jj) >= Pxmin && xe(1,jj) <= Pxmax && xe(2,jj) >= Pymin && xe(2,jj) <= Pymax || qq == 900%|| closest > last_closest
            %disp('Evader evaded capture');
            flag = -1;
        elseif flag == 1
            disp('Capture Successful');
        end
    end
    
    
    
    if flag == 1 || flag == -1
        break
    end
    
    if vid == true
        writeVideo(vidObj, getframe(gcf));
    end
end

if vid == true
    close(vidObj);
end
%figure(2);
%plot([1:qq]*dt, h(~isnan(h)));
%ylabel('min(h_i_j)');
%xlabel('Time');

plot(xpur(:,1), xpur(:,2), 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '--','DisplayName','Pursuer Trajectories');
plot(xpur(:,3:2:end), xpur(:,4:2:end), 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '--','HandleVisibility','off');
plot(xeva(:,1), xeva(:,2), 'LineWidth', 2, 'Color', 'red', 'LineStyle', '--','DisplayName','Evader Trajectory');

legend('Location', 'northwest')
%figure(3);
%plot([1:qq]*dt, umax(~isnan(umax)));
%ylabel('min(u_i)');
%xlabel('Time');
%title('Min control input for fmincon - Pure Coverage')
%grid on