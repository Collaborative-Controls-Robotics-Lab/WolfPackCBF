%% An attempt to reproduce WolfPACK Results by Alexander Davydov

% June 2020 - JHU APL
% Using unicycle dynamics

%%
rng(140)
% Number of pursuers
N = 3;
% Integration time
dt = 0.05;
% Capture distance
epsilon = 0.05;
k = 1;

% Domain
Dxmin = -0.6;
Dxmax = 0.8;
Dymin = -0.6;
Dymax = 0.6;

% Target set
Pxmin = 0.5;
Pxmax = 0.6;
Pymin = -0.05;
Pymax = 0.05;
Pxavg = 0.5*(Pxmin + Pxmax);
Pyavg = 0.5*(Pymin + Pymax);

% Maximum control inputs
vmax = 0.01;
wmax = 2*pi;
velRatio = 0.8;

% Final time
tf = 100;
% Number of iterations
T = tf/dt;

% State Vectors
% Pursuers
xp = zeros(3,N);
% Evader
xe = zeros(3,1);

% Set positions and headings
xe = [-0.45; 0; 0];
xp = [0.4 0.4 0.4; 0.3 -0.06 -0.3; pi pi pi+0.2];
adj = [0 1 0; 1 0 1; 0 1 0];
%xp = [0.4 0.4; 0.1 -0.3; pi pi]
% Initialize projections onto defense surface
p = zeros(1,N);
pCoord = zeros(2,N);
c = zeros(1,N);
cCoord = zeros(2,N);
pVoronoi = zeros(2,N);

% Plot everything
figure; hold on
% Target set
patch([Pxmax Pxmax Pxmin Pxmin], [Pymax Pymin Pymin Pymax], 'green')
xlim([Dxmin Dxmax]);
ylim([Dymin Dymax]);
% Evader
[xs,ys] = getpatches(xe,epsilon);
h1 = patch(xs, ys, 'red');
%h1 = plot(xe(1), xe(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
% Pursuers
[xs,ys] = getpatches(xp,epsilon);
%h2 = plot(xp(1,:), xp(2,:), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
h2 = patch(xs,ys, 'blue');
CircValues = zeros(3, N+1);
IntersectionPoints = zeros(2, 2);
minPts = zeros(2,N);
minDst = zeros(1,N);

CRSCenter = 0.5*([xe(1);xe(2)] + [Pxavg; Pyavg]);
ECircle = circle(CRSCenter(1),CRSCenter(2), norm(CRSCenter - [Pxavg; Pyavg] + 0.05));
CircValues(:,1) = [CRSCenter(1);CRSCenter(2); norm(CRSCenter - [Pxavg; Pyavg])+ 0.05];
h3 = fimplicit(ECircle, 'r', 'LineWidth', 1.5);
h4 = plot(0,0, 'm', 'LineWidth', 1.5);
h5 = plot(0,0, 'md', 'LineWidth', 2.5);
h6 = plot(0,0, 'gd', 'LineWidth', 2.5);
h7 = plot(0,0, 'b*', 'LineWidth', 2.5);
Apollonius = cell(N,1);
DefenseSurface = cell(N,1);

for ii = 1:N
    rs = norm([xe(1); xe(2)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
    os = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1); xe(2)])/(1 - velRatio^2);
    ApCircle = circle(os(1), os(2), rs);
    CircValues(:,ii+1) = [os(1); os(2); rs];
    Apollonius{ii} = fimplicit(ApCircle, 'b', 'LineWidth', 1.5);
    % Find point on circles for minimum distance
end



% Algorithm loop
for qq = 1:T
    % Preallocate some vectors
    % Radii of Apollonius Circles
    % rs = zeros(1,N);
    % Centers of Apollonius Circles
    % os = zeros(2,N);
    flag = 0;
    % Construct constrained reachable set (estimate)
    % NOTE: TO BE UPDATED
    % Approximated by a rough circle
    CRSCenter = 0.5*([xe(1);xe(2)] + [Pxavg; Pyavg]);
    ECircle = circle(CRSCenter(1),CRSCenter(2), norm(CRSCenter - [Pxavg; Pyavg]) + 0.05);
    CircValues(:,1) = [CRSCenter(1);CRSCenter(2); norm(CRSCenter - [Pxavg; Pyavg])+ 0.05];
    set(h3, 'Function', ECircle);
    
    % Construct Apollonius Circles for each pursuer
    for ii = 1:N
        CircValues(3,ii+1) = norm([xe(1); xe(2)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
        CircValues(1:2,ii+1) = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1); xe(2)])/(1 - velRatio^2);
        ApCircle = circle(CircValues(1,ii+1), CircValues(2,ii+1), CircValues(3,ii+1));
        %xA = [xA xs];
        %yA = [yA ys];
        %rii = atan2((xe(2) - xp(2,ii)),(xe(1) - xp(1,ii)));
        %xp(1,ii) = xp(1,ii) + vmax*cos(xp(3,ii))*dt;
        %xp(2,ii) = xp(2,ii) + vmax*sin(xp(3,ii))*dt;
        %xp(3,ii) = xp(3,ii) + k*(sin(rii - xp(3,ii)))*dt;
        [xps,yps] = getpatches(xp, epsilon);
        set(Apollonius{ii}, 'Function', ApCircle);
        % Find point on circles for minimum distance
        [pt, minD] = CircleMinDist(CircValues(1,ii+1), CircValues(2,ii+1), CircValues(3,ii+1), xe(1), xe(2));
        minPts(:,ii) = pt;
        minDst(ii) = minD;
    end
    % Get index of min distance
    idx = find(minDst == min(min(minDst)));
    minPoint = minPts(:,idx);
    % Projection
    r = [Pxavg; Pyavg] - xe(1:2);
    n = [0 1; -1 0]*r/norm(r);
    if n(1) ~= 0
        DefSurfX = @(t) t;
        DefSurfY = @(t) n(2)/n(1)*(t - minPoint(1)) + minPoint(2);
        DefSurf = @(x) n(2)/n(1)*(x - minPoint(1)) + minPoint(2);
    else
        DefSurfX = @(t) minPoint(1) + 0*t;
        DefSurfY = @(t) t;
        DefSurf = minPoint(1); % constant x-value
    end
    % Find intersection between defense surface and CRS
    if n(1) == 0
        B = -2*CircValues(2,1);
        C = CircValues(2,1)^2 + (minPoint(1) - CircValues(1,1))^2 - CircValues(3,1)^2;
        IntersectionPoints(:,1) = [minPoint(1); 0.5*(-B + sqrt(B^2 - 4*C))];
        IntersectionPoints(:,2) = [minPoint(1); 0.5*(-B - sqrt(B^2 - 4*C))];
    else
        m = n(2)/n(1);
        b = minPoint(2) - m*minPoint(1);
        A = 1 + m^2;
        B = -2*CircValues(1,1) + 2*m*(b - CircValues(2,1));
        C = CircValues(1,1)^2 + (b - CircValues(2,1))^2 - CircValues(3,1)^2;
        IntersectionPoints(:,1) = [(-B + sqrt(B^2 - 4*A*C))/(2*A); DefSurf((-B + sqrt(B^2 - 4*A*C))/(2*A))];
        IntersectionPoints(:,2) = [(-B - sqrt(B^2 - 4*A*C))/(2*A); DefSurf((-B - sqrt(B^2 - 4*A*C))/(2*A))];
    end
    oc = 0.5*(IntersectionPoints(:,1) + IntersectionPoints(:,2));
    % Get endpoints of the domain
    VoronoiStart = min([n.'*(IntersectionPoints(:,1) - oc), n.'*(IntersectionPoints(:,2) - oc)]);
    VoronoiEnd = max([n.'*(IntersectionPoints(:,1) - oc), n.'*(IntersectionPoints(:,2) - oc)]);
    % Get agents projected positions and tessellation
    for ii = 1:N
        p(ii) = n.'*(CircValues(1:2, ii+1) - oc);
        pCoord(:,ii) = oc + n*p(ii);
        if ii == 1
            pVoronoi(:,ii) = [VoronoiStart; 0.5*(p(ii) + p(ii + 1))];
        elseif ii == N
            pVoronoi(:,ii) = [0.5*(p(ii-1) + p(ii)); VoronoiEnd];
        else
            pVoronoi(:,ii) = [0.5*(p(ii-1) + p(ii)); 0.5*(p(ii) + p(ii + 1))];
        end
        c(ii) = 0.5*(pVoronoi(1,ii) + pVoronoi(2,ii));
        cCoord(:,ii) = oc + n*c(ii);
        
%         ui = [0;0];
%         for jj = 1:N
%             if adj(jj,ii) == 1
%                 xij = xp(1:2,ii) - xp(1:2,jj);
%                 nij = norm(xij);
%                 xie = xp(1:2,ii) - xe(1:2);
%                 nie = norm(xie);
%                 xje = xp(1:2,jj) - xe(1:2);
%                 nje = norm(xje);
%                 Delta = velRatio*(nie + nje);
%                 dij = 1/2*(epsilon + Delta);
%                 hi = 5e-4*(nij^2 - dij^2)^2;
%                 low1 = (Delta^2 - nij^2)^4;
%                 low2 = (nij^2 - epsilon^2)^2;
%                 dhi = 5e-4*2*(nij^2 - dij^2)*2*(xij);
%                 dDeltaSquared = 2*velRatio^2*xie + 2*velRatio^2*xie*nje/nie;
%                 dlow1 = 4*(Delta^2 - nij^2)^3*(dDeltaSquared - 2*xij);
%                 dlow2 = 2*(nij^2 - epsilon^2)*2*xij;
%                 dlow = low1*dlow2 + dlow1*low2;
%                 gradient = ((low1*low2)*dhi - hi*dlow)/(low1*low2)^2;
%                 ui = ui - gradient;
%                 wij = 8e-1*(nij - dij)*(nij*((Delta - dij) - (dij - epsilon)) + dij*(Delta + epsilon) - 2*epsilon*Delta)...
%                     /(nij*(Delta - nij^2)*(nij - epsilon)^2);
%                 ui = ui + wij*(xp(1:2,jj) - xp(1:2,ii));
%             end
%         end
%         inputs = [cos(xp(3,ii)) sin(xp(3,ii)); -2*sin(xp(3,ii)) 2*cos(xp(3,ii))]*ui;
        alpha = n(2)*vmax/(1 - velRatio^2);
        beta = n(1)*vmax/(1 - velRatio^2);
        track = c(ii) - p(ii) + n.'*velRatio^2*[vmax/velRatio*cos(xe(3)); vmax/velRatio*sin(xe(3))]/(1 - velRatio^2);
        rii = fminbnd(@(th) abs(alpha*sin(th) + beta*cos(th) - track), 0, 2*pi);
        %ui = inputs + [vmax; 4*(rii - xp(3,ii))];
        xp(1,ii) = xp(1,ii) + vmax*cos(xp(3,ii))*dt;
        xp(2,ii) = xp(2,ii) + vmax*sin(xp(3,ii))*dt;
        if abs(4*(rii - xp(3,ii))) > wmax
            xp(3,ii) = xp(3,ii) + sign(4*(rii - xp(3,ii)))*wmax*dt;    
        else
            xp(3,ii) = xp(3,ii) + 4*(rii - xp(3,ii))*dt;
        end
        
        %if abs(ui(1)) > vmax
        %    ui(1) = vmax*sign(ui(1));
        %end
        %if abs(ui(2)) > wmax
        %    ui(2) = wmax*sign(ui(2));
        %end
        %xp(1,ii) = xp(1,ii) + ui(1)*cos(xp(3,ii))*dt;
        %(2,ii) = xp(2,ii) + ui(1)*sin(xp(3,ii))*dt;
        %xp(3,ii) = xp(3,ii) + ui(2)*dt;
    end
    
    
    %set(h4, 'XData', xA, 'YData', yA);
    set(h4, 'XData', DefSurfX(-1:0.1:1), 'YData', DefSurfY(-1:0.1:1))
    set(h5, 'XData', pCoord(1,:), 'YData', pCoord(2,:));
    set(h6, 'XData', cCoord(1,:), 'YData', cCoord(2,:));
    set(h7, 'XData', IntersectionPoints(1,:), 'YData', IntersectionPoints(2,:));
    xe(1) = xe(1) + vmax/velRatio*cos(xe(3))*dt;
    xe(2) = xe(2) + vmax/velRatio*sin(xe(3))*dt;
    xe(3) = xe(3) + 0.078*randn*wmax/velRatio*dt;
    [xes,yes] = getpatches(xe, epsilon);
    set(h1, 'XData', xes, 'YData', yes);
    set(h2, 'XData', xps, 'YData', yps);
    drawnow limitrate
    %pause(0.05)
    for ii = 1:N
        if norm([xe(1); xe(2)] - [xp(1,ii); xp(2,ii)]) <= epsilon
            flag = 1;
        end
    end
    if xe(1) >= Pxmin && xe(1) <= Pxmax && xe(2) >= Pymin && xe(2) <= Pymax
        disp('Evader evaded capture');
        break
    elseif flag == 1
        disp('Capture Successful');
        break
    end
end
