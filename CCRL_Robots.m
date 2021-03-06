classdef CCRL_Robots < handle
    %CCRL_ROBOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numRobots
        isSimulation
        robotXY
        robotTheta
        robotVelLin
        robotVelTheta
        
        initializing = true;
        achievedInitialPosition = false;
        achievedInitialHeading  = false;
        robotDiameterRob
        robotDiameterSim
        
        controlScaleFactor = 1;
        maxIterations
        indxRobots
        domainBoundaries
        dtSim
        barrierSafetyDist
        safety_radius_tight
        si2uniGain
        angularVelLimit
        linearVelLimit
        positionController
        barrierCertificate
        si2uni
        si2unie
        
        sim2rob
        rob2sim
        scaleRatioSim2Rob
        sim2rob_vel
        
        initialPoseX
        initialPoseY
        initialPoseTheta
        robotThetaBias
        
        ViconClient
        viconPose
        robotIndicesVicon
        K4Drv
        looptimer
        previous_time = 0
        dynamicTime
    end
    
    methods
        function obj = CCRL_Robots(numRobots, isSimulation, Dxmin, Dxmax, Dymin, Dymax, initialPoseXYTheta, dt)
            %CCRL_ROBOTS Construct an instance of this class
            %   Detailed explanation goes here
            
            % Initialize parameters and data fields
            obj.looptimer = tic();
            obj.previous_time = obj.getTime;
            
            obj.numRobots = numRobots;
            obj.isSimulation = isSimulation;
            obj.robotXY = zeros(2, numRobots);
            obj.robotTheta = zeros(2, numRobots);
            obj.initialPoseX = initialPoseXYTheta(1,:);
            obj.initialPoseY = initialPoseXYTheta(2,:);
            obj.initialPoseTheta = initialPoseXYTheta(3,:);
            if dt == 0 || ~isSimulation
                obj.dynamicTime = true;
            else
                obj.dynamicTime = false;
            end
            
            % Transform between workspaces
            feet2meter = unitsratio('meter', 'feet');
            simCenter = [(Dxmin + Dxmax)/2;(Dymin + Dymax)/2];
            robWorkspaceHeight = 12*feet2meter;
            robWorkspaceWidth = 17*feet2meter;
            obj.domainBoundaries = 0.5*[-robWorkspaceWidth, robWorkspaceWidth, -robWorkspaceHeight, robWorkspaceWidth];
            obj.scaleRatioSim2Rob = min(robWorkspaceHeight/(Dymax - Dymin),robWorkspaceWidth/(Dxmax - Dxmin));
            obj.sim2rob = @(xx) cat(1,xx(1,:,:)-simCenter(1),xx(2,:,:)-simCenter(2))*obj.scaleRatioSim2Rob;
            obj.rob2sim = @(xx) cat(1,xx(1,:,:)/obj.scaleRatioSim2Rob+simCenter(1),xx(2,:,:)/obj.scaleRatioSim2Rob+simCenter(2));
            obj.sim2rob_vel = @(xx) xx * obj.scaleRatioSim2Rob;
            
            obj.robotDiameterRob = 0.14;
            obj.robotDiameterSim = obj.robotDiameterRob / obj.scaleRatioSim2Rob;
            
            %% Initialize Control Utilities, Parameters and Gains
            % Select the number of iterations for the experiment.  This value is
            % arbitrary
            obj.maxIterations = 1e7;
            % Ordered list of robot id
            obj.indxRobots = 1:numRobots;
            
            
            % Define the domain boundaries
%             obj.domainBoundaries = r.boundaries;
            % Fixed time step for simulation
            obj.dtSim = dt;
            % Safety distance (barrier ceritificates)
            obj.barrierSafetyDist = 0.5 * obj.robotDiameterRob;
            % Termination condition distance
            obj.safety_radius_tight = 0.5 * obj.robotDiameterRob;
            % Gain for single integrator to unicycle conversion
            obj.si2uniGain = 10;
            % Set actuator limits
            % angularVelLimit = 3*controlMaxSpeed;
            obj.angularVelLimit = 2*pi;
            % linearVelLimit = 0.49;
            % linearVelLimit = 0.6*controlMaxSpeed; % 0.8 ORIGINAL
            obj.linearVelLimit = 0.3;%2/3; % 0.8 ORIGINAL
            % Geneterate a go-to-goal controller utility
            obj.positionController = create_si_position_controller;
            % Generate a barrier certificate safety wrap around controller utility
            obj.barrierCertificate = create_si_barrier_certificate('SafetyRadius', obj.barrierSafetyDist);
            % barrierCertificate = @(x,y) x;
            % Generate a single integrator model to unicycle model control conversion
            obj.si2uni = ...
                create_si_to_uni_mapping('ProjectionDistance', 0.5*obj.robotDiameterRob);
            obj.si2unie = ...
                create_si_to_uni_mapping('ProjectionDistance', 0.15*obj.robotDiameterRob);
            
            
            if isSimulation
                obj.robotXY = initialPoseXYTheta(1:2,:);
                obj.robotTheta = initialPoseXYTheta(3,:);
            else
                disp('>>>> Establishing connections with Vicon...')
                [obj.ViconClient,numTrackables] = ViconClientInit;
                disp('>>>> Connections successful!')
                
                % Fetch initial data to create pose object
                obj.viconPose = ViconClientFetch(obj.ViconClient);
                
                % Get the names of the trackable objects
                viconNames = {obj.viconPose.name};
                
                % Determine which objects are named K**
                obj.robotIndicesVicon = contains(viconNames,'K','IgnoreCase',true);
                
                % Determine the IP address of the robots
                robotIP = cell2mat(viconNames(obj.robotIndicesVicon).');
                robotIP = robotIP(:,2:3);
                
                % Determine number of robots
                obj.numRobots = size(robotIP,1);
                disp(['>>>> Number of robots detected: ',num2str(numRobots)])
                disp(['>>>> Evader robot: K', robotIP(numRobots,:)])
                
                % Initialize the robots
                robotDriverVerbosity = 0;
                disp( '>>>> Establishing connections with the robots...')
                obj.K4Drv = KheperaDriverInit(robotIP,[],robotDriverVerbosity);
                disp('>>>> Connections successful!')
                
                % Get initial positions and headings
                robotPose = obj.viconPose(obj.robotIndicesVicon);
                robotPositionXYZ0 = [robotPose.positionXYZ];
                robotAnglesXYZ0 = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                robotXY0 = robotPositionXYZ0(1:2,:);
                robotTheta0 = robotAnglesXYZ0(3,:); % Yaw
                
                % Send a stop command to the robots and wait for a moment
                KheperaSetSpeed(obj.K4Drv, 0, 0)
                pause(0.1)
                disp('>>>> Determining heading bias...')
                
                % Recalibrate to rule out bias. Make the robot move for a short time
                cal_lim_v = obj.linearVelLimit/2;
                KheperaSetSpeed(obj.K4Drv, cal_lim_v, 0)
                pause(0.5)
                KheperaSetSpeed(obj.K4Drv, 0, 0)
                
                % Get updated data
                viconPose = ViconClientFetch(obj.ViconClient,obj.viconPose,1);
                robotPose = viconPose(obj.robotIndicesVicon);
                robotPositionXYZ = [robotPose.positionXYZ];
%               robotAnglesXYZ = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                obj.robotXY = robotPositionXYZ(1:2,:);
                obj.robotTheta = atan2(obj.robotXY(2,:)-robotXY0(2,:),obj.robotXY(1,:)-robotXY0(1,:));
                obj.robotThetaBias = atan2(sin(robotTheta0-obj.robotTheta),cos(robotTheta0-obj.robotTheta));
                disp(['>>>>>>>> Robot vicon heading (degrees):',num2str(robotTheta0*180/pi)])
                disp(['>>>>>>>> Robot actual heading (degrees):',num2str(obj.robotTheta*180/pi)])
                disp(['>>>>>>>> Robot heading bias (degrees):',num2str(obj.robotThetaBias*180/pi)])
            end
            
            
            
            
            
        end
        
        function [obj, XYpos, Theta] = getXYTheta(obj)
            if obj.isSimulation
                XYpos = obj.robotXY;
                Theta = obj.robotTheta;
            else
                % Grabbing new data
                obj.viconPose = ViconClientFetch(obj.ViconClient,obj.viconPose,1);
                robotPose = obj.viconPose(obj.robotIndicesVicon);
                robotPositionXYZ = [robotPose.positionXYZ];
                robotAnglesXYZ = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                obj.robotXY = robotPositionXYZ(1:2,:);
                obj.robotTheta = robotAnglesXYZ(3,:)-obj.robotThetaBias; % Yaw
                
                
                XYpos = obj.rob2sim(obj.robotXY);
                Theta = obj.robotTheta;
            end
        end
        
        function ret = beginVisualization(obj,inputArg)
            % Visualization Elements
            % Load projector calibration matrix
            load('projectiveTransform.mat','H')
            
            % Construct an extended display object
            extDisp = extendedDisplay(H,r.figure_handle);
            
            % Allocate pursuer related handles
%             attDomHandle     = gobjects(2,numPursuers);
%             apolCircHandle   = gobjects(2,numPursuers);
%             apolArcHandle    = gobjects(2,numPursuers);
%             apolIntHandle    = gobjects(2,numPursuers-1);
%             projPosHandle    = gobjects(2,numPursuers);
%             crefPosHandle    = gobjects(2,numPursuers);
%             refHeadingHandle = gobjects(2,numPursuers);
%             % Define a distinct color matrix
%             colorMat = linspecer(numPursuers);
%             % Initialize plot handles
%             reachHandle = extDisp.fill(nan,nan,0.5*[1 1 1],'FaceAlpha',0.05,...
%                 'EdgeColor',0.85*[1 1 1],'EdgeAlpha',0.5,'LineWidth',2);
%             constVelHandle   = extDisp.plot(nan,nan,'--','Color',0.85/2*[1 1 1]);
%             rachSampHandle = extDisp.plot(nan,nan,'bo');
%             for ii = 1:numPursuers
%                 attDomHandle(:,ii) = extDisp.fill(nan,nan,colorMat(ii,:)*0.5,'FaceAlpha',.15,...
%                     'EdgeColor',colorMat(ii,:),'EdgeAlpha',0.5,'LineWidth',2);
% 
%                 apolCircHandle(:,ii) = extDisp.fill(nan,nan,colorMat(ii,:)*0.5,'FaceAlpha',0.15,'EdgeColor',colorMat(ii,:)*0.5,'EdgeAlpha',0.25,'LineWidth',1);
% 
%                 apolArcHandle(:,ii) = extDisp.plot(nan(1,100),nan(1,100),'Color',colorMat(ii,:),'LineWidth',3);
%                 projPosHandle(:,ii) = extDisp.plot(nan,nan,'o','Color',colorMat(ii,:),'MarkerSize',10);
%                 crefPosHandle(:,ii) = extDisp.plot(nan,nan,'x','Color',colorMat(ii,:),'MarkerSize',10);
%                 refHeadingHandle(:,ii) = extDisp.plot(nan,nan,'--','Color',colorMat(ii,:),'LineWidth',2);
% 
%                 if ii < numPursuers
%                     apolIntHandle(:,ii) = extDisp.scatter(nan(1,2),nan(1,2),40,'filled','MarkerEdgeColor',0.66*(colorMat(ii,:)+colorMat(ii+1,:)),'MarkerFaceColor',0.66*(colorMat(ii,:)+colorMat(ii+1,:)),'LineWidth',1.5);
%                 end
% 
%             end
        end
        
        function [obj, dxu] = goToInitialPositions(obj, XYPos, Theta)
            % Initialize robot positions and heading
            if ~obj.achievedInitialPosition
                % Go to goal
                dxi = obj.controlScaleFactor*obj.positionController(XYPos,[obj.initialPoseX; obj.initialPoseY]);
                % Add collision avoidance
                dxi = obj.barrierCertificate(dxi,[XYPos;Theta]);
                % Convert to unicycle model
                dxu = obj.si2uni(dxi,[XYPos;Theta]);
                % Check if they've converged
                robotsDone = hypot(XYPos(1,:)-obj.initialPoseX,XYPos(2,:)-obj.initialPoseY)<obj.robotDiameterSim/2;
                % Set LED's to yellow if they've converged
%                 r.set_left_leds(obj.indxRobots(robotsDone),[255;255;0]*ones(1,nnz(robotsDone)));
%                 r.set_right_leds(obj.indxRobots(robotsDone),[255;255;0]*ones(1,nnz(robotsDone)));
%                 % Set LED's to red if the have not converged
%                 r.set_left_leds(indxRobots(~robotsDone),[255;0;0]*ones(1,nnz(~robotsDone)));
%                 r.set_right_leds(indxRobots(~robotsDone),[255;0;0]*ones(1,nnz(~robotsDone)));

                if all(robotsDone)
                    obj.achievedInitialPosition = true;
                end
            else
                if ~obj.achievedInitialHeading
                    % Turn to your corresponding headings
%                     dxu(2,1:numPursuers) = pi-map_angle(obj.robotTheta(1:numPursuers));
                    dxu(2,:) = sin(obj.initialPoseTheta(:)-Theta(:));
                    % Check if they've converged
                    robotsDone = abs(dxu(2,:))<1*pi/180; % Within 1 degree
                    % Set LED's to green if they've converged
%                     r.set_left_leds(indxRobots(robotsDone),[0;255;0]*ones(1,nnz(robotsDone)));
%                     r.set_right_leds(indxRobots(robotsDone),[0;255;0]*ones(1,nnz(robotsDone)));
%                     % Set LED's to yellow if the have not converged
%                     r.set_left_leds(indxRobots(~robotsDone),[255;255;0]*ones(1,nnz(~robotsDone)));
%                     r.set_right_leds(indxRobots(~robotsDone),[255;255;0]*ones(1,nnz(~robotsDone)));

                    if all(robotsDone) % Within 1 degree
                        obj.achievedInitialHeading = true;
                    end
                else

                    % Done initializing, start a clock using odometry
                    %
                    obj.initializing = false;
                    obj.looptimer = tic();
                    obj.previous_time = 0;
                    dxu = zeros(2, obj.numRobots);

                    % Initializing CRS variables and timer

%                     % Odometry & Reachable set parameters
%                     xf = mean(robotXY(1,1:numPursuers));
%                     yf = mean(robotXY(2,1:numPursuers));
%                     x0 = robotXY(1,numRobots);
% 
%                     tfinal =  1.2 * (xf-x0)/(robotHorizontalSpeed); % Adding 20% to the time required to reach origin w/o maneuvering
%                     % Keeping for now to not screw up reach set calc
%                     %                 tcoverage = tfinal/2; % Assumes same speed from all robots
% 
%                     if simulate_flag
%                         t = 0;
%                         step_cnt = 1;
%                         dt = dtSim;
%                     else
%                         % Create a timer object to determine execution time
%                         timer = tic;
%                         tLast = toc(timer);
%                         % Force rendering and introduce a small delay
%                         drawnow;
%                         % Obtain new time
%                         t = toc(timer);
%                         dt = t-tLast;
%                     end
                end
            end
        end
        
        function obj = set_velocities(obj, dxu)
            % Saturate the control to avoid actuator limits
%             dxu(1,[dxu(1,:)>obj.linearVelLimit false]) = obj.linearVelLimit;
%             dxu(1,[dxu(1,:)<-obj.linearVelLimit false])= -obj.linearVelLimit;
%             dxu(2,dxu(2,:)>obj.angularVelLimit) = obj.angularVelLimit;
%             dxu(2,dxu(2,:)<-obj.angularVelLimit) = -obj.angularVelLimit;
            
            
            if obj.isSimulation
                obj.robotVelLin = dxu(1,:)/obj.scaleRatioSim2Rob;
                obj.robotVelTheta = dxu(2, :);
            else
                % Send velocity updates to robots
                KheperaSetSpeed(obj.K4Drv, dxu(1,:), dxu(2,:))
                % KheperaSetSpeed(K4Drv,0,0)
            end
        
        end
        
        function obj = step(obj)
            if obj.dynamicTime
                t = obj.getTime();
                obj.dtSim = t - obj.previous_time;
                obj.previous_time = t;
            else
                obj.previous_time = obj.previous_time + obj.dtSim;
            end
            
            if obj.isSimulation
                obj.robotTheta = obj.robotTheta + obj.robotVelTheta * obj.dtSim;
                obj.robotXY = obj.robotXY + [cos(obj.robotTheta); sin(obj.robotTheta)] .* obj.robotVelLin * obj.dtSim;
            end
        end
        
        function obj = stop(obj)
            if ~obj.isSimulation
                % Stop all robots at the end
                KheperaSetSpeed(obj.K4Drv, 0, 0)
                % Disconnect the Vicon client
                ViconClientDisconnect(obj.ViconClient)
                % Disconnect the Khepera driver client
                KheperaDriverDisconnect(obj.K4Drv)
                % Close the "stop" window it is still open
            end

        end
        
        function t = getTime(obj)
            if ~obj.dynamicTime
                t = obj.previous_time;
            else
                t = toc(obj.looptimer);
            end
            
        end
        
    end
end

