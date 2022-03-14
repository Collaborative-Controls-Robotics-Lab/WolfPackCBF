%% An attempt to reproduce WolfPACK Results by Alexander Davydov
% Modified by Mikhail Khrenov

% March 2021 - UMD CCRL
% Using unicycle dynamics

%%
close all
vid = false;
isSimulation = true;
rng(140)

% Number of evaders
numEvaders = 2;
% Integration time
dt = 0;

% Capture distance
epsilon = 0.05;
k = 1;
cvxflag = 0;

% Domain
Dxmin = -0.5;
Dxmax = 1.2;
Dymin = -0.6;
Dymax = 0.6;
% Dxmin = -2.5908;
% Dxmax = 2.5908;
% Dymin = -1.8288;
% Dymax = 1.8288;

% Target set
Pxmin = 0.925;
Pxmax = 0.975;
Pymin = -0.025;
Pymax = 0.025;
Pxavg = 0.5*(Pxmin + Pxmax);
Pyavg = 0.5*(Pymin + Pymax);

% Maximum control inputs
vmax = 0.02;
wmax = 2*pi;
closest = 100000;
velRatio = 0.4;%0.5;

% Final time
TF = 30;
t_final = TF*ones(numEvaders);
% Number of iterations
max_iterations = 1e4; %% HACK


% Set initial evader positions and headings
xe = [-0.4*ones(1, numEvaders); linspace(-0.3, 0.3,numEvaders); zeros(1,numEvaders)];
% xe = [-0.3; 0; 0];

% Terminal time
t_terminal = norm(xe(1:2, 1) - [Pxavg; Pyavg])/(vmax*(1 + 1/velRatio));
% Dist to span at terminal time
L = ((vmax/velRatio*(t_final(1) - t_terminal))^2 - (vmax*t_terminal)^2)/(vmax/velRatio)/(t_final(1) - t_terminal)

% Number of pursuers
numPursuers = 5;%numAgents(velRatio, epsilon, L);%4;%
L_pred = distSpanned(velRatio, epsilon, numPursuers)


%% Initialize the System


% State Vectors
% Pursuers
xp = zeros(3,numPursuers);
adj = zeros(numPursuers, numPursuers);


% Pursuer-Evader Map
pem = zeros(numPursuers,numEvaders);
pem = [1, 0; 
       1, 0; 
       1, 0; 1, 0; 1, 0; 1, 0; 1, 0; 1, 0;];

% set initial pursuer position
%(-vmax*t_terminal + Pxavg)
xp = [Pxmin*ones(1,numPursuers) - 0.05; linspace(-0.5, 0.5, numPursuers); pi*ones(1,numPursuers)];
%xp = [0.4*ones(1,N) + 0.1*rand(1,N); linspace(-0.13, 0.15, N); pi*ones(1,N)];

% Virtual pursuers on the CRS boundary
xr = nan(2, 2*numEvaders);
nOnes = ones(numPursuers,1);
adj = diag(nOnes(1:numPursuers-1), -1) + diag(nOnes(1:numPursuers-1), 1);
%xp = [0.4 0.4; 0.1 -0.3; pi pi]

% Initialize projections onto defense surface
p = zeros(1,numPursuers);
pCoord = zeros(2,numPursuers);
c = zeros(1,numPursuers);
cCoord = zeros(2,numPursuers);
pVoronoi = zeros(2,numPursuers);
sigmaC = nan(1,max_iterations);
h = nan(1,max_iterations);
umax = nan(1,max_iterations);
u = zeros(2*numPursuers,1);
unorm = zeros(numPursuers,1);
Aij = zeros(numPursuers-1, 2*numPursuers);
Aie = zeros(numPursuers+1, 2);
bie = zeros(numPursuers+1, 1);
bij = zeros(1, numPursuers-1);
hij = zeros(1, numPursuers-1);
flag = 0;
midpt = zeros(2, numPursuers-1);
He = cell(numPursuers+2,1);
ke = cell(numPursuers+2,1);
de = cell(numPursuers+2,1);

neigh = zeros(numPursuers, 2);


r = CCRL_Robots(numPursuers + numEvaders, isSimulation, Dxmin, Dxmax, Dymin, Dymax, [xp, xe], dt);



if vid == true
    cam = webcam('Logitech BRIO');
    cam.Resolution = cam.AvailableResolutions{end-3}; % Highest resolution
    frameCropY = 161:820; % Depends of camera positioning
    
    frames2write(max_iterations) = struct('cdata',[],'colormap',[]);
    
    vidObj = VideoWriter(['test_',datestr(datetime('now'),'yyyymmdd_HHMMSS')], 'MPEG-4');
    frameCount = 1;
end

% Plot everything
hFig = figure(1);
set(hFig,'color','w');
hold on; 
axis equal off;
% xlim([Dxmin Dxmax]);
% ylim([Dymin Dymax]);
domBounds = r.domainBoundaries;
xlim(domBounds(1:2));
ylim(domBounds(3:4));
%% Visualization Elements
% Load projector calibration matrix
load('projectiveTransform.mat','H')
% Construct an extended display object
extDisp = extendedDisplay(H,hFig);


%%
if ~r.isSimulation
    set(hFig,'Units','normalized','color','w','Position',[1 1/3 2/3 2/3]) % on TV
else
    set(hFig,'color','w','Position',[100 180 700*1.8 360*1.8])
end


% Target set
targetSet = r.sim2rob([Pxmax Pxmax Pxmin Pxmin; Pymax Pymin Pymin Pymax]);
extDisp.patch(targetSet(1,:), targetSet(2,:), 'FaceColor', 'green', 'HandleVisibility','off');


% Evader
[xes,yes] = getpatches(xe,0.05);
hEvaders = gobjects(numEvaders, 1);
xyes = r.sim2rob([xes(:),yes(:)].');
xes = reshape(xyes(1,:),3,[]);
yes = reshape(xyes(2,:),3,[]);

for jj = 1:numEvaders
    hEvaders(jj,:) = patch(extDisp.hAxesTV, 'XData', xes(jj,:), 'YData', yes(jj,:), 'FaceColor', 'red', 'DisplayName', 'Evader');
end

% Pursuers
[xps,yps] = getpatches(xp,0.05);
hPursuers = gobjects(numPursuers, 1);
xyps = r.sim2rob([xps(:),yps(:)].');
xps = reshape(xyps(1,:),3,[]).';
yps = reshape(xyps(2,:),3,[]).';

for ii = 1:numPursuers
    hPursuers(ii,:) = patch(extDisp.hAxesTV, 'XData', xps(ii,:), 'YData', yps(ii,:), 'FaceColor', 'blue', 'DisplayName', 'Pursuers');
end
CircValues = zeros(3, numPursuers, numEvaders);
IntersectionPoints = zeros(2, 2);
minPts = zeros(2,numPursuers);
minDst = zeros(1,numPursuers);
numParameterizationPoints = 100;
CRS = zeros(2,numParameterizationPoints,numEvaders);
thetaParameterEllipse = linspace(0,2*pi,numParameterizationPoints);
% Path recording
xpur = nan(max_iterations,2*numPursuers);
xeva = nan(max_iterations,2);
ue = zeros(2, numEvaders);

hEvaderCRS = gobjects(numEvaders,2);

for jj = 1:numEvaders
    %CRS{jj} = ellipse(xe(1:2,jj), [Pxavg; Pyavg], t_final(jj)*vmax/velRatio);
    %h3(:, jj) = extDisp.fimplicit(CRS{jj}, 'r', 'LineWidth', 1.5, 'DisplayName','Reachable Set');
    [CRS(1,:,jj),CRS(2,:,jj)] = parEllipse(thetaParameterEllipse,xe(1:2,jj), [Pxavg; Pyavg],t_final(jj)*vmax/velRatio);
    CRS(:,:,jj) = r.sim2rob(CRS(:,:,jj));
    hEvaderCRS(jj, :) = extDisp.patch(CRS(1,:,jj),CRS(2,:,jj), 'FaceColor','r', 'EdgeColor','r', 'LineWidth', 1.5, 'FaceAlpha',0.25);
end

xr_plot = r.sim2rob(xr);
hReferencePointsOnEllipse = extDisp.plot(xr_plot(1,:),xr_plot(2,:), 'md', 'LineWidth', 2.5, 'DisplayName','Reference Point');

hApollonius = gobjects(numPursuers,2,numEvaders);
hEpsilonCirc = gobjects(numPursuers, 2);
% DefenseSurface = gobjects(N,2);

ngon_N = 16;
ngon_theta = linspace(0, 2*pi, ngon_N + 1).';
ngon_theta(end) = [];

ngon_H = [cos(ngon_theta), sin(ngon_theta)];
ngon_k = vmax * cos(pi/ngon_N) * ones(ngon_N, 1);

for ii = 1:numPursuers
    for jj = 1:numEvaders
        rs = norm([xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)])*velRatio/(1 - velRatio^2);
        os = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1,jj); xe(2,jj)])/(1 - velRatio^2);
        %ApCircle = circle(os(1), os(2), rs);
        [ApCircleX,ApCircleY] = parCircle(thetaParameterEllipse, os, rs);
        CircValues(:,ii,jj) = [os(1); os(2); rs];
        ApCircle = r.sim2rob([ApCircleX(:),ApCircleY(:)].');
        hApollonius(ii,:,jj) = extDisp.patch(ApCircle(1,:),ApCircle(2,:), 'FaceColor', 'b', 'EdgeColor', 'b', 'LineWidth', 1.5,'FaceAlpha',0.05);
        % Find point on circles for minimum distance
        
    end
    % plot epsilon circle
    % epcirc = circle(xp(1,ii), xp(2,ii), epsilon);
    [epcircX,epcircY] = parCircle(thetaParameterEllipse,xp(1:2,ii), epsilon);
    epcirc = r.sim2rob([epcircX(:),epcircY(:)].');
    hEpsilonCirc(ii,:) = extDisp.patch(epcirc(1,:), epcirc(2,:), 'EdgeColor', 'g', 'LineWidth', 1,'FaceColor','none');
end

% Evader distance at terminal time
r_terminal = vmax * t_terminal;
% nom_circle = circle(Pxavg, Pyavg, r_terminal);
[nom_circleX,nom_circleY] = parCircle(thetaParameterEllipse, [Pxavg;Pyavg], r_terminal);
nom_circle = r.sim2rob([nom_circleX(:),nom_circleY(:)].');
hDefenseCircle = extDisp.patch(nom_circle(1,:), nom_circle(2,:), 'EdgeColor', 'm', 'LineWidth', 1.5,'FaceColor','none');

inpose = [xp(1:2,:).'; xe(1:2,:).'].';
initialPoseX = inpose(1,:);
initialPoseY = inpose(2, :);


%% Coalition Forming Step
pem = zeros(numPursuers, numEvaders);
isendpoint = zeros(numPursuers, numEvaders);
ii = 1;
for jj = 1:numEvaders
    % Terminal time
    t_terminal = norm(xe(1:2, jj) - [Pxavg; Pyavg])/(vmax*(1 + 1/velRatio));
    % Dist to span at terminal time
    L = ((vmax/velRatio*(t_final(jj) - t_terminal))^2 - (vmax*t_terminal)^2)/(vmax/velRatio)/(t_final(jj) - t_terminal);
    
    % Number of pursuers
    numP = numAgents(velRatio, epsilon, L);
    n = 0;

    isendpoint(ii, jj) = 1;
    while n < numP
        pem(ii, jj) = 1;
        n = n + 1;
        ii = ii + 1;
    end
    isendpoint(ii-1,jj) = -1;
end

r1 = zeros(2, numPursuers);
r2 = zeros(2, numPursuers);

%% Algorithm loop
for qq = 1:max_iterations
    % Update position and orientation
    [r, xy_all, theta_all] = r.getXYTheta();
    xp = [xy_all(:, 1:end-numEvaders); theta_all(:, 1:end-numEvaders)];
    xe = [xy_all(:, numPursuers+1:end); theta_all(:, numPursuers+1:end)];
    
    while r.initializing
        [r, xy_all, theta_all] = r.getXYTheta();
        xp = [xy_all(:, 1:end-numEvaders); theta_all(:, 1:end-numEvaders)];
        xe = [xy_all(:, numPursuers+1:end); theta_all(:, numPursuers+1:end)];
        [r, dxu] = r.goToInitialPositions(xy_all, theta_all);
        r = r.set_velocities(dxu); 
        r = r.step();
        [xps,yps] = getpatches(xp, 0.05);
        [xes,yes] = getpatches(xe, 0.05);
        xyes = r.sim2rob([xes(:),yes(:)].');%[xes(:).';yes(:).'];%
        xyps = r.sim2rob([xps(:),yps(:)].');%[xps(:).';yps(:).'];%
        
        for ii = 1:numEvaders 
            set(hEvaders(ii), 'XData', xyes(1,3*(ii-1)+1:3*ii), 'YData', xyes(2,3*(ii-1)+1:3*ii));
        end
%         set(hEvaders, 'XData', xyes(1,:), 'YData', xyes(2,:));
        xps = reshape(xyps(1,:),3,[]).';
        yps = reshape(xyps(2,:),3,[]).';
        for ii = 1:numPursuers 
            set(hPursuers(ii), 'XData', xps(ii,:), 'YData', yps(ii, :));
        end
        drawnow limitrate
    end
    
    t = r.getTime();
    tgo = t_final - t;
    
    
    
    %% Main control loop
    u = zeros(2, numPursuers);
    for jj = 1:numEvaders
        %CRS{jj} = ellipse(xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        [CRS(1,:,jj),CRS(2,:,jj)] = parEllipse(thetaParameterEllipse, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        CRS(:,:,jj) = r.sim2rob(CRS(:,:,jj));
        extDisp.set(hEvaderCRS(jj,:), CRS(1,:,jj),CRS(2,:,jj));
        %set(hEvaderCRS(jj), 'Function', CRS{jj});
    end

    % Construct nominal control points
    for jj = 1:numEvaders
        slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2, jj));
    %     poin = (-r_terminal * ([Pxavg; Pyavg] - xe(1:2, jj))/norm([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
    %     r1 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, 1);
    %     r2 = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, N);
    
        poin = (-r_terminal * ([Pxavg; Pyavg] - xe(1:2, jj))/norm([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
        if norm([Pxavg; Pyavg] - xe(1:2, jj)) > norm([Pxavg; Pyavg] - poin.')
            r1(:,jj) = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, 1);
            r2(:,jj) = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, numPursuers);
        else
            poin = (-0.9 * ([Pxavg; Pyavg] - xe(1:2, jj))).' + [Pxavg Pyavg];
            r1(:,jj) = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, 1);
            r2(:,jj) = lineEllipse(poin, slope, xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio, numPursuers);
        end
    end

    % Construct Apollonius Circles for each pursuer/evader pair
    for ii = 1:numPursuers
        for jj = 1:numEvaders
            diffi = [xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)];
            di = norm(diffi); % distance between pursuer and evader
            gammai = atan2(-diffi(2),-diffi(1)); % angle from pursuer to evader
            CircValues(3,ii) = (di*velRatio + epsilon)/(1 - velRatio^2); % revised to Andrew's method
            CircValues(1:2,ii) = (velRatio*epsilon + di)/(1 - velRatio^2)*[cos(gammai);sin(gammai)] + xe(1:2,jj);
%             CircValues(1:2,ii) = ([xp(1,ii); xp(2,ii)] - velRatio^2*[xe(1,jj); xe(2,jj)])/(1 - velRatio^2);
            %ApCircle = circle(CircValues(1,ii), CircValues(2,ii), CircValues(3,ii));
            [ApCircleX,ApCircleY] = parCircle(thetaParameterEllipse, CircValues(1:2,ii), CircValues(3,ii));
            ApCircle = r.sim2rob([ApCircleX(:),ApCircleY(:)].');
            extDisp.set(hApollonius(ii,:, jj), ApCircle(1,:),ApCircle(2,:));
            
            
        end
        [epcircX,epcircY] = parCircle(thetaParameterEllipse, xp(1:2,ii), epsilon);
        epcirc = r.sim2rob([epcircX(:),epcircY(:)].');
        extDisp.set(hEpsilonCirc(ii,:), epcirc(1,:), epcirc(2,:));
        
    end
    %tic

    % CBF QCQP for evader
    for ii = 1:numPursuers
        for jj = 1:numEvaders
            Aie(ii,:) = -(xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii));
            He{ii} = zeros(2);
            ke{ii} = Aie(ii,:).';
            bie(ii) = -vmax + (norm(xe(1:2,jj) - xp(1:2,ii)) - 1.05*epsilon);
            de{ii} = -bie(ii);
        end
    end

    uvmax = vmax/velRatio;
    He{numPursuers+1} = zeros(2);
    Aie(numPursuers+1,:) = (xe(1:2, 1) - [Pxavg; Pyavg]).'/norm(xe(1:2, 1) - [Pxavg; Pyavg]);
    ke{numPursuers+1} = Aie(numPursuers+1,:).';
    bie(numPursuers+1) = -uvmax + 1e2*tgo(jj)*uvmax - norm(xe(1:2, 1) - [Pxavg; Pyavg])^3;
    de{numPursuers+1} = -bie(numPursuers+1);

    He{numPursuers+2} = 2*eye(2);
    ke{numPursuers+2} = zeros(2,1);
    de{numPursuers+2} = -uvmax^2;


    % Solve QCQP
    % Compute alpha 1 angle for each pair
    angles = zeros(numPursuers+1,1);
    for ii = 1:numPursuers+1
        for jj = 1:numEvaders
            if ii == 1
                [x,y,~] = closestEllipse(CircValues(1:2,ii), xe(1:2, jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                xstar = [x; y];
                xei = norm(xe(1:2, jj) - xp(1:2,ii));
                xestar = norm(xe(1:2, jj) - xstar);
                xistar = norm(xp(1:2,ii) - xstar);
                % Angle via law of cosines
                angles(ii) = acos((xistar^2 - xestar^2 - xei^2)/(-2*xei*xestar));
            elseif ii == numPursuers+1
                [x,y,~] = closestEllipse(CircValues(1:2,ii-1), xe(1:2, jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                xstar = [x; y];
                xei = norm(xe(1:2, jj) - xp(1:2,numPursuers));
                xestar = norm(xe(1:2, jj) - xstar);
                xistar = norm(xp(1:2,numPursuers) - xstar);
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
        [x,y,lambda] = closestEllipse(CircValues(1:2,idx), xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        xstar = [x; y];
        [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
        num = xp(1:2,idx) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
        xie = xp(1:2,idx) - xe(1:2,jj);
        unom = velRatio*xie/norm(xie) + ((velRatio^2 - 1)*dxe - velRatio^2*eye(2)).'*num/norm(num);
        unom = unom/norm(unom)*uvmax;
    elseif idx == numPursuers+1
        [x,y,lambda] = closestEllipse(CircValues(1:2,idx-1), xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        xstar = [x; y];
        [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
        num = xp(1:2,numPursuers) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
        xie = xp(1:2,numPursuers) - xe(1:2,jj);
        unom = velRatio*xie/norm(xie) + ((velRatio^2 - 1)*dxe - velRatio^2*eye(2)).'*num/norm(num);
        unom = unom/norm(unom)*uvmax;
    else
        xie = xp(1:2,idx-1) - xe(1:2,jj);
        xje = xp(1:2,idx) - xe(1:2,jj);
        unom = xie/norm(xie) + xje/norm(xje);
        unom = unom/norm(unom)*uvmax;
    end

    for jj = 1:numEvaders

        if tgo(jj)*uvmax <= 1.08*norm(xe(1:2,jj) - [Pxavg; Pyavg])
            ue(:,jj) = ([Pxavg; Pyavg] - xe(1:2,jj))/norm(([Pxavg; Pyavg] - xe(1:2,jj)))*uvmax;
            [t_final, tgo, xe] = resetEvader(jj, xe, r, numPursuers, t_final, TF, tgo, t);
%             flag = 1; %% COMMENTED OUT FOR "ASTEROIDS"
        else
%             nonlconstr = @(x)quadconstr(x,He,ke,de);
            ue(:,jj) = ([Pxavg; Pyavg] - xe(1:2,jj))/norm(([Pxavg; Pyavg] - xe(1:2,jj)))*uvmax;

            % If nominal solution does not meet our constraints, solve CBF 
            % with linear programming
%             if any(nonlconstr(ue) > 0)

                quad_A = [Aie(1:end-1, :); ngon_H];
                quad_b = [bie(1:end-1, :); ngon_k / velRatio];
                options = optimset('Display', 'off');
                ue(:,jj) = linprog(xe(1:2,jj) - [Pxavg; Pyavg], quad_A, quad_b, [], [], [], [], [0; 0], options);

%             end

        end
    end

    for ii = 1:numPursuers
        for jj = 1:numEvaders

            if pem(ii, jj) > 0.5

                % Find two nearest pursuers
                if qq == 1
                    f_dist = arrayfun(@(r_idx) norm(xp(1:2,r_idx) - xp(1:2, ii)), 1:numPursuers);
                    for iii = 1:numPursuers
                        if pem(iii, jj) < 0.5
                            f_dist(iii) = inf;
                        end
                    end
                    [val, indices] = sort(f_dist);
                    if1 = indices(2);
                    f1 = xp(1:2, indices(2));
                    if numPursuers > 2
                        if2 = indices(3);
                        f2 = xp(1:2, indices(3));
                    else
                        if2 = 1;
                        f2 = zeros(2,1);
                    end
                    
                    neigh(ii, 1) = if1;
                    neigh(ii, 2) = if2;
                else
                    if1 = neigh(ii, 1);
                    if2 = neigh(ii, 2);
                    f1 = xp(1:2, if1);
                    f2 = xp(1:2, if2);
                end

                if numPursuers == 1
                    slope = [0 -1; 1 0]*([Pxavg; Pyavg] - xe(1:2, jj));
                    % Find both intersection points
                    [A,B,C,D,E,F] = ellipseData(xe(1:2, jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
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
                else

                    if isendpoint(ii, jj) == 1
                        xr(:,(jj*2)) = r1(:, jj);

                        w1 = 1/(CircValues(3,ii) + 0);
                        w2 = 1/(CircValues(3,if1) + 0);
                        a1 = (CircValues(3,ii) + 0);
                        a2 = (CircValues(3,if1) + 0);

                        ui = 0.5*xr(:,(jj*2)) + 0.5*w1/(w1 + w2)*CircValues(1:2,ii) ...
                            + 0.5*w2/(w1 + w2)*CircValues(1:2,if1) - CircValues(1:2,ii);
                        ui = (1 - velRatio^2)*ui + velRatio^2*ue(:,jj);
                        Aij = -(velRatio*(xp(1:2,ii) - xe(1:2, jj)).'/norm(xp(1:2,ii) - xe(1:2, jj)) ...
                            - (xp(1:2,ii) - f1).'/norm(xp(1:2,ii) - f1));
                        bij = a1/(a1 + a2)*(velRatio*((xe(1:2, jj) - xp(1:2,ii)).'/norm(xe(1:2, jj) - xp(1:2,ii)) + (xe(1:2, jj) - f1).'/norm(xe(1:2, jj) - f1)) ...
                            * ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2, jj)) + norm(f1 - xe(1:2, jj))) - norm(xp(1:2,ii) - f1) + epsilon)^3);
                        % A1 - CRS constraint
                        [x,y,lambda] = closestEllipse(CircValues(1:2,ii), xe(1:2, jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                        [A,B,C,D,E,F] = ellipseData(xe(1:2, jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                        xstar = [x;y];
                        dxi = dxstardxi(A,B,C,D,E,F,x,y,lambda,velRatio);
                        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
                        dxt = dxstardt(A,B,C,D,E,F,x,y,lambda,vmax/velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
                        num = xp(1:2,ii) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
                        xie = xp(1:2,ii) - xe(1:2,jj);
                        dhdxi = velRatio*xie.'/norm(xie) - num.'/norm(num)*(eye(2) + (velRatio^2 - 1)*dxi);
                        dhdxe = -velRatio*xie.'/norm(xie) - num.'/norm(num)*((velRatio^2 - 1)*dxe - velRatio^2*eye(2));
                        dhdxstar = (1 - velRatio^2)*num.'/norm(num);
                        A1 = -dhdxi - dhdxstar*dxi;
                        b1 = dhdxe*ue(:,jj) + dhdxstar*(dxe*ue(:,jj) + dxt) + 1e3*(velRatio*norm(xie) - norm(num) + epsilon)^3;

                        H = {zeros(2), zeros(2), 2*eye(2)};
                        k = {Aij', A1', zeros(2,1)};
                        d = {-bij, -b1, -vmax^2};   
                    elseif isendpoint(ii, jj) == -1
                        xr(:,2*jj-1) = r2(:,jj);

                        w1 = 1/(CircValues(3,if1) + 0);
                        w2 = 1/(CircValues(3,ii) + 0);
                        a1 = (CircValues(3,if1) + 0);
                        a2 = (CircValues(3,ii) + 0);

                        ui = 0.5*xr(:,2*jj-1) + 0.5*w2/(w2 + w1)*CircValues(1:2,ii) ...
                            + 0.5*w1/(w2 + w1)*CircValues(1:2,if1) - CircValues(1:2,ii);
                        ui = (1 - velRatio^2)*ui + velRatio^2*ue(:,jj);

                        Aij = -(velRatio*(xp(1:2,ii) - xe(1:2,jj)).'/norm(xp(1:2,ii) - xe(1:2,jj)) ...
                            - (xp(1:2,ii) - f1).'/norm(xp(1:2,ii) - f1));
                        bij = a1/(a1 + a2)*(velRatio*((xe(1:2,jj) - xp(1:2,ii)).'/norm(xe(1:2,jj) - xp(1:2,ii)) + (xe(1:2,jj) - f1).'/norm(xe(1:2,jj) - f1)) ...
                            *ue(:,jj) + 1e4*(velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f1 - xe(1:2,jj))) - norm(xp(1:2,ii) - f1) + epsilon)^3);
                        [x,y,lambda] = closestEllipse(CircValues(1:2,ii), xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                        [A,B,C,D,E,F] = ellipseData(xe(1:2,jj), [Pxavg; Pyavg], tgo(jj)*vmax/velRatio);
                        xstar = [x;y];
                        dxi = dxstardxi(A,B,C,D,E,F,x,y,lambda,velRatio);
                        dxe = dxstardxe(A,B,C,D,E,F,x,y,lambda,velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
                        dxt = dxstardt(A,B,C,D,E,F,x,y,lambda,vmax/velRatio,[Pxavg;Pyavg],xe(1:2,jj),tgo(jj)*vmax/velRatio);
                        num = xp(1:2,ii) - xstar + velRatio^2*xstar - velRatio^2*xe(1:2,jj);
                        xie = xp(1:2,ii) - xe(1:2,jj);
                        dhdxi = velRatio*xie.'/norm(xie) - num.'/norm(num)*(eye(2) + (velRatio^2 - 1)*dxi);
                        dhdxe = -velRatio*xie.'/norm(xie) - num.'/norm(num)*((velRatio^2 - 1)*dxe - velRatio^2*eye(2));
                        dhdxstar = (1 - velRatio^2)*num.'/norm(num);
                        A1 = -dhdxi - dhdxstar*dxi;
                        b1 = dhdxe*ue(:,jj) + dhdxstar*(dxe*ue(:,jj) + dxt) + 1e3*(velRatio*norm(xie) - norm(num) + epsilon)^3;

                        H = {zeros(2), zeros(2), 2*eye(2)};
                        k = {Aij', A1', zeros(2,1)};
                        d = {-bij, -b1, -vmax^2};   
                    else
                        Aij = zeros(2);
                        bij = zeros(2,1);
                        w0 = 1/(CircValues(3,if1) + 0);
                        w1 = 1/(CircValues(3,ii) + 0);
                        w2 = 1/(CircValues(3,if2) + 0);
                        a0 = (CircValues(3,if1) + 0);
                        a1 = (CircValues(3,ii) + 0);
                        a2 = (CircValues(3,if2) + 0);

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
                    end                    

                    nonlconstr = @(x)quadconstr(x,H,k,d);

                    % If solution given by weighted consensus does not meet
                    % our constraints, solve CBF with quadratic programming
                    if norm(ui) > vmax || any(nonlconstr(ui) > 0)
                        quad_A = [Aij; A1; ngon_H];
                        quad_b = [bij; b1; ngon_k];
                        options = optimset('Display', 'off');
                        z = quadprog(eye(2), -ui, quad_A, quad_b, [], [], [], [], [0; 0], options);
                    else
                        z = ui;
                    end

                    % If solution is infeasible, fall back to saturating
                    % weighted Consensus result
                    if norm(z) > vmax && cvxflag == 0
                        z = ui/norm(ui)*vmax;
%                         cvxflag = 1;
                    end

                    u(:,ii) = u(:,ii) + pem(ii,jj) * z;
                    hij(ii) = velRatio*(norm(xp(1:2,ii) - xe(1:2,jj)) + norm(f2 - xe(1:2,jj))) - norm(xp(1:2,ii) - f2)+ 2*epsilon;                   
                end
            end
        end
    end

%         if simulate_flag
%             % Updating step and time
%             step_cnt = step_cnt + 1;
%             t = t + dt;
%         else
%             step_cnt = round(t/dtSim);
%             t_last = t;
%             t = toc(timer);
%             dt = t-t_last;
%         end

    dxi = [u ue];
    dxip = u;
    dxie = ue;
    dxip = r.sim2rob_vel(dxip);
    dxie = r.sim2rob_vel(dxie);
    dxup = r.si2uni(dxip,xp);
    dxue = r.si2unie(dxie,xe);
%     dxu = r.si2uni(dxi,[r.robotXY;r.robotTheta]);
    dxu = [dxup dxue]; 

%     toc
    
     %% Send velocities to agents
     % Convert to unicycle model
     

    r = r.set_velocities(dxu); 
    r = r.step();
    
    
    % record pursuer positions
    for ii = 1:numPursuers
        xpur(qq, 2*ii-1:2*ii) = xp(1:2,ii).';
    end
    
    % record evaders positions
    for jj = 1:numEvaders
        xeva(qq,:) = xe(1:2,jj).';
    end
        
    h(qq) = min(hij);
    for ii = 1:numPursuers
        unorm(ii) = norm(u(:,ii));
    end
    umax(qq) = min(unorm);

    xr_plot = r.sim2rob(xr);
    extDisp.set(hReferencePointsOnEllipse,  xr_plot(1,:), xr_plot(2,:));
    [xps,yps] = getpatches(xp, 0.05);
    [xes,yes] = getpatches(xe, 0.05);
    
    xyes = r.sim2rob([xes(:),yes(:)].');%[xes(:).';yes(:).']; %
    xyps = r.sim2rob([xps(:),yps(:)].');%;[xps(:).';yps(:).']; %
    
    xps = reshape(xyps(1,:),3,[]).';
    yps = reshape(xyps(2,:),3,[]).';
    for ii = 1:numPursuers 
        set(hPursuers(ii), 'XData', xps(ii,:), 'YData', yps(ii, :));
    end

    for ii = 1:numEvaders 
        set(hEvaders(ii), 'XData', xyes(1,3*(ii-1)+1:3*ii), 'YData', xyes(2,3*(ii-1)+1:3*ii));
    end
%     set(hEvaders, 'XData', xyes(1,:), 'YData', xyes(2,:));        
    
%     set(h1, 'XData', xes, 'YData', yes);
%     set(h2, 'XData', xps, 'YData', yps);
    drawnow limitrate
    
    last_closest = closest;
    dx = xp(1:2,1:numPursuers) - xe(1:2,1); 
    closest = 100000;
    
%     if ~(r.isSimulation)
%         pause(dt - 0.01)
%     end
    
    for jj = 1:numEvaders
        for ii = 1:numPursuers
            ev_dist = norm(dx(1:2,ii));
            if ev_dist < closest
                closest = ev_dist;
            end
            
            if norm([xe(1,jj); xe(2,jj)] - [xp(1,ii); xp(2,ii)]) <= epsilon
                [t_final, tgo, xe] = resetEvader(jj, xe, r, numPursuers, t_final, TF, tgo, t)
%                 flag = 1; %% Commented out for "Asteroids" implementation
            end
        end
        
        if xe(1,jj) >= Pxmin && xe(1,jj) <= Pxmax && xe(2,jj) >= Pymin && xe(2,jj) <= Pymax %|| qq == 900%|| closest > last_closest
            disp('Evader evaded capture');
            flag = -1;
        elseif flag == 1
            disp('Capture Successful');
        end
    end
    
    
    
    if flag == 1 || flag == -1
        break
    end
    
    if vid == true
        frame = cam.snapshot;
        frames2write(frameCount) = im2frame(frame(frameCropY,:,:));
        %
        frameCount = frameCount + 1;
    end
end

%figure(2);
%plot([1:qq]*dt, h(~isnan(h)));
%ylabel('min(h_i_j)');
%xlabel('Time');

% hFigStopLoop.Clear();
% clear hFigStopLoop;
r = r.stop();

xpur_plot = permute(xpur(~isnan(xpur(:,1)),1:2:end),[3,2,1]);
ypur_plot = permute(xpur(~isnan(xpur(:,1)),2:2:end),[3,2,1]);
xypur_plot = cat(1,xpur_plot,ypur_plot);
xypur_plot = r.sim2rob(xypur_plot);
xeva_plot = r.sim2rob(xeva(~isnan(xeva(:,1)),:).').';
extDisp.plot(squeeze(xypur_plot(1,1,:)),squeeze(xypur_plot(2,1,:)), 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '--','DisplayName','Pursuer Trajectories');
extDisp.plot([permute(xypur_plot(1,2:end,:),[3,2,1]);nan(1,numPursuers-1)], [permute(xypur_plot(2,2:end,:),[3,2,1]);nan(1,numPursuers-1)], 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '--','HandleVisibility','off');
extDisp.plot(xeva_plot(:,1), xeva_plot(:,2), 'LineWidth', 2, 'Color', 'red', 'LineStyle', '--','DisplayName','Evader Trajectory');
extDisp.legend('Location', 'northwest');

if vid == true
    frame = cam.snapshot;
    frames2write(frameCount) = im2frame(frame(frameCropY,:,:));
    vidObj.Quality = 75;
    vidObj.FrameRate = round(frameCount/t);
    frames2write = frames2write(1:frameCount);
    
    open(vidObj);
    writeVideo(vidObj, frames2write);
    close(vidObj);
    clear('cam')
end

%figure(3);
%plot([1:qq]*dt, umax(~isnan(umax)));
%ylabel('min(u_i)');
%xlabel('Time');
%title('Min control input for fmincon - Pure Coverage')
%grid on

function [t_final, tgo, xe] = resetEvader(jj, xe, r, numPursuers, t_final, TF, tgo, t)
    h = pi/2*rand(1) - pi/4;
    xy = [0.5; 0]-0.9*[cos(h); sin(h)];
    xe(:,jj) = [xy; h];
    r.robotXY(:, numPursuers + jj) = xe(1:2, jj);
    t_final(jj) = t + TF;
    tgo = t_final - t;
end