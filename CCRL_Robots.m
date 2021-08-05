classdef CCRL_Robots
    %CCRL_ROBOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numRobots
        isSimulation
        robotXY
        robotTheta
        
        initializing = true;
        achievedInitialPosition = false;
        achievedInitialHeading  = false;
        robotDiameter = 0.14;
    end
    
    methods
        function obj = CCRL_Robots(numRobots,isSimulation)
            %CCRL_ROBOTS Construct an instance of this class
            %   Detailed explanation goes here
            
            % Initialize parameters and data fields
            obj.numRobots = numRobots;
            obj.isSimulation = isSimulation;
            obj.robotXY = zeros(2, numRobots);
            obj.robotTheta = zeros(2, numRobots);
            
            if isSimulation
                obj.robotXY = diag([robWorkspaceHalfWidth, robWorkspaceHalfHeight]) * (2*rand(2, numRobots) - 1);
                obj.robotTheta = 2*pi*rand(1,numRobots);
            else
                disp('>>>> Establishing connections with Vicon...')
                [ViconClient,numTrackables] = ViconClientInit;
                disp('>>>> Connections successful!')
                
                % Fetch initial data to create pose object
                viconPose = ViconClientFetch(ViconClient);
                
                % Get the names of the trackable objects
                viconNames = {viconPose.name};
                
                % Determine which objects are named K**
                robotIndicesVicon = contains(viconNames,'K','IgnoreCase',true);
                
                % Determine the IP address of the robots
                robotIP = cell2mat(viconNames(robotIndicesVicon).');
                robotIP = robotIP(:,2:3);
                
                % Determine number of robots
                obj.numRobots = size(robotIP,1);
                disp(['>>>> Number of robots detected: ',num2str(numRobots)])
                disp(['>>>> Evader robot: K', robotIP(numRobots,:)])
                
                % Initialize the robots
                robotDriverVerbosity = 0;
                disp( '>>>> Establishing connections with the robots...')
                K4Drv = KheperaDriverInit(robotIP,[],robotDriverVerbosity);
                disp('>>>> Connections successful!')
                
                % Get initial positions and headings
                robotPose = viconPose(robotIndicesVicon);
                robotPositionXYZ0 = [robotPose.positionXYZ];
                robotAnglesXYZ0 = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                robotXY0 = robotPositionXYZ0(1:2,:);
                robotTheta0 = robotAnglesXYZ0(3,:); % Yaw
                
                % Send a stop command to the robots and wait for a moment
                KheperaSetSpeed(K4Drv, 0, 0)
                pause(0.1)
                disp('>>>> Determining heading bias...')
                
                % Recalibrate to rule out bias. Make the robot move for a short time
                KheperaSetSpeed(K4Drv, controlMaxSpeed/2, 0)
                pause(0.5)
                KheperaSetSpeed(K4Drv, 0, 0)
                
                % Get updated data
                viconPose = ViconClientFetch(ViconClient,viconPose,1);
                robotPose = viconPose(robotIndicesVicon);
                robotPositionXYZ = [robotPose.positionXYZ];
%               robotAnglesXYZ = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                obj.robotXY = robotPositionXYZ(1:2,:);
                obj.robotTheta = atan2(obj.robotXY(2,:)-robotXY0(2,:),obj.robotXY(1,:)-robotXY0(1,:));
                robotThetaBias = atan2(sin(robotTheta0-robotTheta),cos(robotTheta0-robotTheta));
                disp(['>>>>>>>> Robot vicon heading (degrees):',num2str(robotTheta0*180/pi)])
                disp(['>>>>>>>> Robot actual heading (degrees):',num2str(robotTheta*180/pi)])
                disp(['>>>>>>>> Robot heading bias (degrees):',num2str(robotThetaBias*180/pi)])
            end
            
        end
        
        function XYpos, Theta = getXYTheta(obj)
            if obj.isSimulation
                XYpos = obj.robotXY;
                Theta = obj.robotTheta;
            else
                % Grabbing new data
                viconPose = ViconClientFetch(ViconClient,viconPose,1);
                robotPose = viconPose(robotIndicesVicon);
                robotPositionXYZ = [robotPose.positionXYZ];
                robotAnglesXYZ = [robotPose.anglesXYZ];
                
                % Extract controllable states for convenience
                obj.robotXY = robotPositionXYZ(1:2,:);
                obj.robotTheta = robotAnglesXYZ(3,:)-robotThetaBias; % Yaw
                
                
                XYpos = obj.robotXY;
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
    end
end

