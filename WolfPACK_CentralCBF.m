%% An attempt to reproduce WolfPACK Results by Alexander Davydov

% June 2020 - JHU APL
% Using unicycle dynamics

%%

vid = false;
rng(140)
% Number of pursuers
N = 2;
% Integration time
dt = 0.05;
% Capture distance
epsilon = 0.02;
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
velRatio = 0.75;

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
xe = [-0.45; -0; 0];
xp = [0.4*ones(1,N); linspace(-0.3, 0.35, N); pi*ones(1,N)];
% Virtual pursuers on the CRS boundary
xr = zeros(2);
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
bij = zeros(1, N-1);
hij = zeros(1, N-1);
flag = 0;

if vid == true
    %nFrames = 20;
    vidObj = VideoWriter('test', 'MPEG-4');
    vidObj.Quality = 75;
    vidObj.FrameRate = 50;
    open(vidObj);
end

% Plot everything
figure(1); hold on
% Target set
patch([Pxmax Pxmax Pxmin Pxmin], [Pymax Pymin Pymin Pymax], 'green')
xlim([Dxmin Dxmax]);
ylim([Dymin Dymax]);
% Evader
[xs,ys] = getpatches(xe,0.05);
h1 = patch(xs, ys, 'red');
%h1 = plot(xe(1), xe(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
% Pursuers
[xs,ys] = getpatches(xp,0.05);
%h2 = plot(xp(1,:), xp(2,:), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
h2 = patch(xs,ys, 'blue');
CircValues = zeros(3, N+1);
IntersectionPoints = zeros(2, 2);
minPts = zeros(2,N);
minDst = zeros(1,N);

%CRSCenter = 0.5*([xe(1);xe(2)] + [Pxavg; Pyavg]);
%ECircle = circle(CRSCenter(1),CRSCenter(2), norm(CRSCenter - [Pxavg; Pyavg] + 0.05));
%CircValues(:,1) = [CRSCenter(1);CRSCenter(2); norm(CRSCenter - [Pxavg; Pyavg])+ 0.05];
CRS = ellipse(xe(1:2), [Pxavg; Pyavg], T*dt*vmax/velRatio);
h3 = fimplicit(CRS, 'r', 'LineWidth', 1.5);
h4 = plot(0,0, 'm', 'LineWidth', 1.5);
h5 = plot(0,0, 'md', 'LineWidth', 2.5);
%h6 = plot(0,0, 'gd', 'LineWidth', 2.5);
%h7 = plot(0,0, 'b*', 'LineWidth', 2.5);
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
count = 1;

% Algorithm loop
for qq = 1:T
    CRS = ellipse(xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
    set(h3, 'Function', CRS);
    
    
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
        %[xps,yps] = getpatches(xp, epsilon);
        set(Apollonius{ii}, 'Function', ApCircle);
        % Find point on circles for minimum distance
        [pt, minD] = CircleMinDist(CircValues(1,ii+1), CircValues(2,ii+1), CircValues(3,ii+1), xe(1), xe(2));
        minPts(:,ii) = pt;
        minDst(ii) = minD;
    end
     for ii = 1:N
        if N == 1
            slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2));
            % Find both intersection points
            [A,B,C,D,E,F] = ellipseData(xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            p = xp(1:2,ii);
            if slope(1) == 0 % Vertical line x = const
                x = p(1);
                xr(1,:) = [x x];
                a = C;
                b = B*x + E;
                c = A*x^2 + D*x + F;
                % Quadratic formula
                xr(2,:) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
            else
                m = slope(2)/slope(1);
                a = A + m*B + m^2*C;
                b = -m*B*p(1) + B*p(2) - 2*m^2*C*p(1) + 2*m*C*p(2) + D + m*E;
                c = m^2*C*p(1)^2 - 2*m*C*p(1)*p(2) + C*p(2)^2 - m*E*p(1) + E*p(2) + F;
                % Quadratic formula
                xr(1,:) = [(-b + sqrt(b^2 - 4*a*c))/(2*a) (-b - sqrt(b^2 - 4*a*c))/(2*a)];
                xr(2,:) = m*xr(1,:) - m*p(1) + p(2);
            end
            ui = xe(1:2) + xr(:,1) + xr(:,2) - 3*xp(1:2,ii);
        elseif ii == 1
            %[x,y] = closestEllipse(xp(1:2,ii), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %[x,y] = closestEllipse(CircValues(1:2,ii+1), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %xr(:,1) = [x;y];
            slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2));
            xr(:,1) = lineEllipse(xp(1:2,ii), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio,1);
            %ui = xr(:,1) + xp(1:2,ii+1) + xe(1:2) - 3*xp(1:2,ii);
            %ui = xr(:,1) + xp(1:2,ii+1) - 2*xp(1:2,ii);
            %ui = (1 - velRatio^2)*(xr(:,1) + xp(1:2,ii+1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3));sin(xe(3))];
            w1 = 1;%/CircValues(3,ii+1);
            w2 = 1;%/CircValues(3,ii+2);
            %ui = 0.5*xr(:,1) + 0.5*w1/(w1 + w2)*xp(1:2,ii) ...
            %    + 0.5*w2/(w1 + w2)*xp(1:2,ii+1) - xp(1:2,ii);
            ui = 0.5*xr(:,1) + 0.5*w1/(w1 + w2)*CircValues(1:2,ii+1) ...
                + 0.5*w2/(w1 + w2)*CircValues(1:2,ii+2) - CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio*vmax*[cos(xe(3)); sin(xe(3))];
            Aij(1,1:2) = -(velRatio*(xp(1:2,ii) - xe(1:2)).'/norm(xp(1:2,ii) - xe(1:2)) ...
                - (xp(1:2,ii) - xp(1:2,ii+1)).'/norm(xp(1:2,ii) - xp(1:2,ii+1)));
            Aij(1,3:4) = -(velRatio*(xp(1:2,ii+1) - xe(1:2)).'/norm(xp(1:2,ii+1) - xe(1:2)) ...
                - (xp(1:2,ii+1) - xp(1:2,ii)).'/norm(xp(1:2,ii) - xp(1:2,ii+1)));
            bij(1) = velRatio*((xe(1:2) - xp(1:2,ii)).'/norm(xe(1:2) - xp(1:2,ii)) + (xe(1:2) - xp(1:2,ii+1)).'/norm(xe(1:2) - xp(1:2,ii+1))) ...
                * [cos(xe(3)); sin(xe(3))] + 1e6*(velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1)) + 2*epsilon)^3;
            hij(1) = velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1)) + 2*epsilon;
        elseif ii == N
            %[x,y] = closestEllipse(xp(1:2,ii), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %[x,y] = closestEllipse(CircValues(1:2,ii+1), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %xr(:,2) = [x;y];
            slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2));
            xr(:,2) = lineEllipse(xp(1:2,ii), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio,N);
            %ui = xr(:,2) + xp(1:2, ii-1) + xe(1:2) - 3*xp(1:2,ii);
            %ui = xr(:,2) + xp(1:2, ii-1) - 2*xp(1:2,ii);
            %ui = (1 - velRatio^2)*(xr(:,2) + xp(1:2,ii-1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3)); sin(xe(3))];
            w1 = 1;%/CircValues(3,ii);
            w2 = 1;%/CircValues(3,ii+1);
            %ui = 0.5*xr(:,2) + 0.5*w2/(w2 + w1)*xp(1:2,ii) ...
            %    + 0.5*w1/(w2 + w1)*xp(1:2,ii-1) - xp(1:2,ii);
            ui = 0.5*xr(:,2) + 0.5*w2/(w2 + w1)*CircValues(1:2,ii+1) ...
                + 0.5*w1/(w2 + w1)*CircValues(1:2,ii) - CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio*vmax*[cos(xe(3)); sin(xe(3))];

        else
            %ui = xp(1:2,ii - 1) + xp(1:2,ii + 1) + xe(1:2) - 3*xp(1:2,ii);
            w0 = 1;%/CircValues(3,ii);
            w1 = 1;%/CircValues(3,ii+1);
            w2 = 1;%/CircValues(3,ii+2);
            %ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*xp(1:2,ii) ...
            %    + 0.5*w0/(w0 + w1)*xp(1:2,ii-1) ...
            %    + 0.5*w2/(w2 + w1)*xp(1:2,ii+1) + xe(1:2) - 2*xp(1:2,ii);
            ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*CircValues(1:2,ii+1) ...
                + 0.5*w0/(w0 + w1)*CircValues(1:2,ii) ...
                + 0.5*w2/(w2 + w1)*CircValues(1:2,ii+2) - CircValues(1:2,ii+1);%+ xe(1:2) - 2*CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio*vmax*[cos(xe(3)); sin(xe(3))] + xe(1:2) - xp(1:2,ii);
            Aij(ii,2*ii-1:2*ii) = -(velRatio*(xp(1:2,ii) - xe(1:2)).'/norm(xp(1:2,ii) - xe(1:2)) ...
                - (xp(1:2,ii) - xp(1:2,ii+1)).'/norm(xp(1:2,ii) - xp(1:2,ii+1)));
            Aij(ii,2*ii+1:2*ii+2) = -(velRatio*(xp(1:2,ii+1) - xe(1:2)).'/norm(xp(1:2,ii+1) - xe(1:2)) ...
                - (xp(1:2,ii+1) - xp(1:2,ii)).'/norm(xp(1:2,ii) - xp(1:2,ii+1)));
            bij(ii) = velRatio*((xe(1:2) - xp(1:2,ii)).'/norm(xe(1:2) - xp(1:2,ii)) + (xe(1:2) - xp(1:2,ii+1)).'/norm(xe(1:2) - xp(1:2,ii+1))) ...
                * [cos(xe(3)); sin(xe(3))] + 1e6*(velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1)) + 2*epsilon)^3;
            hij(ii) = velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1))+ 2*epsilon;
        end
        
        %ui = (ui)/norm(ui)*vmax;
        u(2*ii - 1:2*ii) = ui;
        %xp(1:2,ii) = xp(1:2,ii) + (ui) * dt;
     end
    uOpt = quadprog(eye(2*N), -u, Aij, bij);%, [], [], ones(2*N,1)*-vmax/sqrt(2), ones(2*N,1)*vmax/sqrt(2));
    if norm(uOpt - u) >= 10^-6
        count = 3;
    end
    for ii = 1:N
        xp(1:2,ii) = xp(1:2,ii) + uOpt(2*ii-1:2*ii)*dt;
    end
    h(qq) = min(hij);
    %umax(qq) = max(uOpt);
    for ii = 1:N
        unorm(ii) = norm(uOpt(2*ii-1:2*ii));
    end
    umax(qq) = max(unorm);
    
    %set(h4, 'XData', xA, 'YData', yA);
    %set(h4, 'XData', DefSurfX(-1:0.1:1), 'YData', DefSurfY(-1:0.1:1))
    set(h5, 'XData', xr(1,:), 'YData', xr(2,:));
    %set(h6, 'XData', cCoord(1,:), 'YData', cCoord(2,:));
    %set(h7, 'XData', IntersectionPoints(1,:), 'YData', IntersectionPoints(2,:));
    [xps,yps] = getpatches(xp, 0.05);
    xe(1) = xe(1) + vmax/velRatio*cos(xe(3))*dt;
    xe(2) = xe(2) + vmax/velRatio*sin(xe(3))*dt;
    xe(3) = xe(3) + 0.07*randn*wmax/velRatio*dt;
    %xe(3) = xe(3) + wmax/velRatio*dt;
    [xes,yes] = getpatches(xe, 0.05);
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
    if vid == true
        writeVideo(vidObj, getframe(gcf));
    end
end
if vid == true
    close(vidObj);
end
figure(2);
plot([1:qq]*dt, h(~isnan(h)));
ylabel('min(h_i_j)');
xlabel('Time');

figure(3);
plot([1:qq]*dt, umax(~isnan(umax)));
ylabel('max(u_i)');
xlabel('Time');