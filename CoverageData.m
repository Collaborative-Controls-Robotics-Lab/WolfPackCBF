%% An attempt to reproduce WolfPACK Results by Alexander Davydov

% June 2020 - JHU APL
% Using unicycle dynamics

%%
conversion = 50000;
totalmisses = ones(100,1);
for kk = 1
vid = false;
rng(kk)
% Number of pursuers
N = 4;
% Integration time
dt = 0.05;
% Capture distance
epsilon = 0.05;
k = 1;
cvxflag = 0;

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
velRatio = 0.5;

% Final time
tf = 65;
% Number of iterations
T = tf/dt;

% State Vectors
% Pursuers
xp = zeros(3,N);
% Evader
xe = zeros(3,1);

% Set positions and headings
%xe = [-0.45; 0.05; 0];
% Fix initial evader position
xe = [-0.45; 0; 0];
% Sample from uniform distribution in half ellipse
halfpoint = 1/2*(xe(1:2) + [Pxavg; Pyavg]);
xrange = 0.5*(T*dt*vmax/velRatio);
yrange = 2*sqrt(xrange^2 - norm(xe(1:2) - halfpoint)^2);
for ii = 1:N
    initflag = 0;
    while initflag == 0
        xp(:,ii) = [halfpoint(1) + xrange*rand(1); halfpoint(2) + yrange*rand(1) - yrange/2; pi];
        if norm(xe(1:2) - xp(1:2,ii)) + norm(xp(1:2,ii) - [Pxavg; Pyavg]) <= T*dt*vmax/velRatio
            break
        end
    end
end
% Sort pursuers by y-position
xp = sortrows(xp.', 2).';
%xp = [0.4*ones(1,N) + 0.05*randn(1,N); linspace(-0.2, 0.2, N) + 0.075*randn(1,N); pi*ones(1,N)];
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
Aie = zeros(N+1, 2);
bie = zeros(N+1, 1);
bij = zeros(1, N-1);
hij = zeros(1, N-1);
flag = 0;
midpt = zeros(2, N-1);
He = cell(N+2,1);
ke = cell(N+2,1);
de = cell(N+2,1);

if vid == true
    %nFrames = 20;
    vidObj = VideoWriter('test', 'MPEG-4');
    vidObj.Quality = 75;
    vidObj.FrameRate = 50;
    open(vidObj);
end

% Plot everything
%figure(1); hold on
% Target set
%patch([Pxmax Pxmax Pxmin Pxmin], [Pymax Pymin Pymin Pymax], 'green')
th = 0:pi/50:2*pi;
xunit = epsilon*cos(th) + Pxavg;
yunit = epsilon*sin(th) + Pyavg;
%fill(xunit, yunit, 'g');
%xlim([Dxmin Dxmax]);
%ylim([Dymin Dymax]);
% Evader
%[xs,ys] = getpatches(xe,0.05);
%h1 = patch(xs, ys, 'red');
%h1 = plot(xe(1), xe(2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
% Pursuers
%[xs,ys] = getpatches(xp,0.05);
%h2 = plot(xp(1,:), xp(2,:), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
%h2 = patch(xs,ys, 'blue');
CircValues = zeros(3, N+1);
IntersectionPoints = zeros(2, 2);
minPts = zeros(2,N);
minDst = zeros(1,N);
miss = 100*ones(T, N);

%CRSCenter = 0.5*([xe(1);xe(2)] + [Pxavg; Pyavg]);
%ECircle = circle(CRSCenter(1),CRSCenter(2), norm(CRSCenter - [Pxavg; Pyavg] + 0.05));
%CircValues(:,1) = [CRSCenter(1);CRSCenter(2); norm(CRSCenter - [Pxavg; Pyavg])+ 0.05];
CRS = ellipse(xe(1:2), [Pxavg; Pyavg], T*dt*vmax/velRatio);
%h3 = fimplicit(CRS, 'r', 'LineWidth', 1.5);
%h4 = plot(0,0, 'm', 'LineWidth', 1.5);
%h5 = plot(0,0, 'md', 'LineWidth', 2.5);
%h6 = plot(0,0, 'gd', 'LineWidth', 2.5);
%h7 = plot(0,0, 'b*', 'LineWidth', 2.5);
Apollonius = cell(N,1);
DefenseSurface = cell(N,1);

for ii = 1:N
    rs = norm([xe(1); xe(2)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
    os = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1); xe(2)])/(1 - velRatio^2);
    ApCircle = circle(os(1), os(2), rs);
    CircValues(:,ii+1) = [os(1); os(2); rs];
    %Apollonius{ii} = fimplicit(ApCircle, 'b', 'LineWidth', 1.5);
    % Find point on circles for minimum distance
end
%count = 1;

% Algorithm loop
for qq = 1:T
    CRS = ellipse(xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
    %set(h3, 'Function', CRS);
    
    
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
        %set(Apollonius{ii}, 'Function', ApCircle);
        % Find point on circles for minimum distance
        %[pt, minD] = CircleMinDist(CircValues(1,ii+1), CircValues(2,ii+1), CircValues(3,ii+1), xe(1), xe(2));
        %minPts(:,ii) = pt;
        %minDst(ii) = minD;
    end
    tic
    
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
            %[xr(:,1),~] = lineEllipse(0.5*(xe(1:2) + [Pxavg; Pyavg]), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            xr(:,1) = lineEllipse(xp(1:2,ii), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio, 1);
            
            %ui = xr(:,1) + xp(1:2,ii+1) + xe(1:2) - 3*xp(1:2,ii);
            %ui = xr(:,1) + xp(1:2,ii+1) - 2*xp(1:2,ii);
            %ui = (1 - velRatio^2)*(xr(:,1) + xp(1:2,ii+1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3));sin(xe(3))];
            w1 = 1;%/CircValues(3,ii+1);
            w2 = 1;%/CircValues(3,ii+2);
            a1 = CircValues(3,ii+1);
            a2 = CircValues(3,ii+2);
            %ui = 0.5*xr(:,1) + 0.5*w1/(w1 + w2)*xp(1:2,ii) ...
            %    + 0.5*w2/(w1 + w2)*xp(1:2,ii+1) - xp(1:2,ii);
            if isreal(xr(:,1))
            ui = 0.5*xr(:,1) + 0.5*w1/(w1 + w2)*CircValues(1:2,ii+1) ...
                + 0.5*w2/(w1 + w2)*CircValues(1:2,ii+2) - CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio^2*ue;
            else
                ui = vmax*(xe(1:2) - xp(1:2,ii))/norm(xe(1:2) - xp(1:2,ii));
            end
            
            z = ui/norm(ui)*vmax;
            
            u(1:2) = z;
            %toc
            %hij(1) = velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1)) + 2*epsilon;
            %midpt(:,ii) = w1/(w1 + w2)*xp(1:2,ii) + w2/(w1 + w2)*xp(1:2,ii+1);
        elseif ii == N
            %[x,y] = closestEllipse(xp(1:2,ii), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %[x,y] = closestEllipse(CircValues(1:2,ii+1), xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            %xr(:,2) = [x;y];
            slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2));
            %[~,xr(:,2)] = lineEllipse(0.5*(xe(1:2) + [Pxavg; Pyavg]), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio);
            xr(:,2) = lineEllipse(xp(1:2,ii), slope, xe(1:2), [Pxavg; Pyavg], (T-qq)*dt*vmax/velRatio,2);
            %ui = xr(:,2) + xp(1:2, ii-1) + xe(1:2) - 3*xp(1:2,ii);
            %ui = xr(:,2) + xp(1:2, ii-1) - 2*xp(1:2,ii);
            %ui = (1 - velRatio^2)*(xr(:,2) + xp(1:2,ii-1) - 2*CircValues(1:2,ii+1)) + velRatio*vmax*[cos(xe(3)); sin(xe(3))];
            w1 = 1;%/CircValues(3,ii);
            w2 = 1;%/CircValues(3,ii+1);
            a1 = CircValues(3,ii);
            a2 = CircValues(3,ii+1);
            %ui = 0.5*xr(:,2) + 0.5*w2/(w2 + w1)*xp(1:2,ii) ...
            %    + 0.5*w1/(w2 + w1)*xp(1:2,ii-1) - xp(1:2,ii);
            if isreal(xr(:,2))
            ui = 0.5*xr(:,2) + 0.5*w2/(w2 + w1)*CircValues(1:2,ii+1) ...
                + 0.5*w1/(w2 + w1)*CircValues(1:2,ii) - CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio^2*ue;
            else
                ui = vmax*(xe(1:2) - xp(1:2,ii))/norm(xe(1:2) - xp(1:2,ii));
            end
            
            z = ui/norm(ui)*vmax;
            
            u(2*ii-1:2*ii) = z;
        else
            %ui = xp(1:2,ii - 1) + xp(1:2,ii + 1) + xe(1:2) - 3*xp(1:2,ii);
            w0 = 1;%/CircValues(3,ii);
            w1 = 1;%/CircValues(3,ii+1);
            w2 = 1;%/CircValues(3,ii+2);
            a0 = CircValues(3,ii);
            a1 = CircValues(3,ii+1);
            a2 = CircValues(3,ii+2);
            %ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*xp(1:2,ii) ...
            %    + 0.5*w0/(w0 + w1)*xp(1:2,ii-1) ...
            %    + 0.5*w2/(w2 + w1)*xp(1:2,ii+1) - xp(1:2,ii);
            ui = 0.5*(w1/(w0 + w1) + w1/(w1 + w2))*CircValues(1:2,ii+1) ...
                + 0.5*w0/(w0 + w1)*CircValues(1:2,ii) ...
                + 0.5*w2/(w2 + w1)*CircValues(1:2,ii+2) - CircValues(1:2,ii+1);%+ xe(1:2) - 2*CircValues(1:2,ii+1);
            ui = (1 - velRatio^2)*ui + velRatio^2*ue;% + xe(1:2) - xp(1:2,ii);
            
            z = ui/norm(ui)*vmax;
            
            u(2*ii-1:2*ii) = z;
            %hij(ii) = velRatio*(norm(xp(1:2,ii) - xe(1:2)) + norm(xp(1:2,ii+1) - xe(1:2))) - norm(xp(1:2,ii) - xp(1:2,ii+1))+ 2*epsilon;
            %midpt(:,ii) = w1/(w1 + w2)*xp(1:2,ii) + w2/(w1 + w2)*xp(1:2,ii+1);
        end

        %ui = (ui)/norm(ui)*vmax;
        %u(2*ii - 1:2*ii) = ui;
        %xp(1:2,ii) = xp(1:2,ii) + (ui) * dt;
    end
    % CBF QCQP for evader
    for ii = 1:N
        Aie(ii,:) = -(xe(1:2) - xp(1:2,ii)).'/norm(xe(1:2) - xp(1:2,ii));
        He{ii} = zeros(2);
        ke{ii} = Aie(ii,:).';
        bie(ii) = -(xe(1:2) - xp(1:2,ii)).'*u(2*ii-1:2*ii)/norm(xe(1:2) - xp(1:2,ii)) + 1e4*(norm(xe(1:2) - xp(1:2,ii)) - 0.05)^3;
        de{ii} = -bie(ii);
    end
    uvmax = vmax/velRatio;
    He{N+1} = zeros(2);
    Aie(N+1,:) = (xe(1:2) - [Pxavg; Pyavg]).'/norm(xe(1:2) - [Pxavg; Pyavg]);
    ke{N+1} = Aie(N+1,:).';
    bie(N+1) = -uvmax + 1e2*(tf - (qq)*dt)*uvmax - norm(xe(1:2) - [Pxavg; Pyavg])^3;
    de{N+1} = -bie(N+1);
    He{N+2} = 2*eye(2);
    ke{N+2} = zeros(2,1);
    de{N+2} = -uvmax^2;
    % Solve QCQP
    
    % Go straight to the goal
    %unom = ([Pxavg; Pyavg] - xe(1:2));
    % Separate agents 1 and 2
    %unom = velRatio*((xp(1:2,1) - xe(1:2))/norm((xp(1:2,1) - xe(1:2))) + (xp(1:2,2) - xe(1:2))/norm(xp(1:2,2) - xe(1:2)));
    % Separate agents 2 and 3
    %unom = ((xp(1:2,2) - xe(1:2))/norm((xp(1:2,2) - xe(1:2))) + (xp(1:2,3) - xe(1:2))/norm(xp(1:2,3) - xe(1:2)));
    %unom = unom/norm(unom)*uvmax;

    if (tf - (qq)*dt)*uvmax <= 1.03*norm(xe(1:2) - [Pxavg; Pyavg])
        ue = ([Pxavg; Pyavg] - xe(1:2))/norm(([Pxavg; Pyavg] - xe(1:2)))*uvmax;
    else
    %cvx_begin quiet
    %    variable ue(2)
    %    minimize( (ue - unom).' * (ue - unom))
    %    subject to
    %        Aie * ue <= bie
    %        ue.'*ue <= uvmax^2
    %cvx_end
    
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'off',...
        'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
        'HessianFcn', @(x,lambda)quadhess(x,lambda,eye(2),He));
    fun = @(x)quadobj(x,zeros(2),xe(1:2) - [Pxavg; Pyavg]);
    nonlconstr = @(x)quadconstr(x,He,ke,de);
    try
    [ue,fval,eflag,output,lambda] = fmincon(fun,[0;0],...
        [],[],[],[],[],[],nonlconstr,options);
    catch
    cvx_begin quiet
        variable ue(2)
        minimize( (xe(1:2) - [Pxavg;Pyavg]).' * (ue) )
        subject to
            Aie * ue <= bie
            ue.'*ue <= uvmax^2
    cvx_end
    end
    end
    toc
    for ii = 1:N
        xp(1:2,ii) = xp(1:2,ii) + u(2*ii-1:2*ii)*dt;
        %xp(3,ii) = atan2(u(2*ii), u(2*ii-1));
        %xp(1:2,ii) = xp(1:2,ii) + vmax*(xe(1:2) - xp(1:2,ii))/norm(xe(1:2) - xp(1:2,ii))*dt;
        %xp(3,ii) = atan2(xe(2) - xp(2,ii), xe(1) - xp(1,ii));
        miss(qq,ii) = norm(xp(1:2,ii) - xe(1:2));
    end
    %h(qq) = min(hij);
    %for ii = 1:N
    %    unorm(ii) = norm(u(2*ii-1:2*ii));
    %end
    %umax(qq) = min(unorm);
    
    %set(h4, 'XData', xA, 'YData', yA);
    %set(h4, 'XData', DefSurfX(-1:0.1:1), 'YData', DefSurfY(-1:0.1:1))
    %set(h5, 'XData', xr(1,:), 'YData', xr(2,:));
    %set(h6, 'XData', midpt(1,:), 'YData', midpt(2,:));
    %set(h7, 'XData', IntersectionPoints(1,:), 'YData', IntersectionPoints(2,:));
    %[xps,yps] = getpatches(xp, 0.05);
    xe(1:2) = xe(1:2) + ue*dt;
    %xe(3) = atan2(ue(2), ue(1));
    %xe(1) = xe(1) + vmax/velRatio*cos(xe(3))*dt;
    %xe(2) = xe(2) + vmax/velRatio*sin(xe(3))*dt;
    %xe(3) = xe(3) + 0.07*randn*wmax/velRatio*dt;
    %xe(3) = xe(3) + wmax/velRatio*dt;
    %[xes,yes] = getpatches(xe, 0.05);
    %set(h1, 'XData', xes, 'YData', yes);
    %set(h2, 'XData', xps, 'YData', yps);
    %drawnow limitrate
    %pause(0.05)
    for ii = 1:N
        if norm([xe(1); xe(2)] - [xp(1,ii); xp(2,ii)]) <= epsilon
            flag = 1;
        end
    end
    if norm(xe(1:2) - [Pxavg; Pyavg]) <= epsilon
        disp('Evader evaded capture');
        break
    %elseif flag == 1
    %    disp('Capture Successful');
    %    break
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

%figure(3);
%plot([1:qq]*dt, umax(~isnan(umax)));
%ylabel('min(u_i)');
%xlabel('Time');
%title('Min control input for fmincon - Pure Coverage')
%grid on
totalmisses(kk) = min(min(miss), [], 2);
end