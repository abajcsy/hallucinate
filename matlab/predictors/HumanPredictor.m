classdef HumanPredictor < handle
    %HUMANPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        human       % (obj) dynamical system representing human model
        xcurr       % (arr) vector of the current [x,y,P(beta=0)] state
        tau         % (arr) time vector representing prediction horizon
        grid        % (obj) compute grid for prediction
        targetR     % (float) radius of initial target set
        
        uMode 
        minWith
        
        % Value function computation information
        preds           % (float arr) Stores most recent prediction (i.e. value fun)
        timeDisc        % (float arr) Discretized time vector          
        schemeData      % (struct) Used to specify dyn, grid, etc. for HJIPDE_solve()
        HJIextraArgs    % (struct) Specifies extra args to HJIPDE_solve()
    end
    
    methods
        %% Constructor
        function obj = HumanPredictor(params)
            obj.human = params.humanModel;
            obj.xcurr = params.x0;
            obj.tau = params.tMin:params.dt:params.tMax;
            obj.grid = params.predGrid;
            obj.targetR = params.targetR;
            obj.uMode = params.uMode;
            obj.minWith = params.minWith;
            
            % Put grid and dynamic systems into schemeData.
            obj.schemeData.grid = obj.grid;
            obj.schemeData.dynSys = obj.human;
            obj.schemeData.accuracy = 'medium'; 
            obj.schemeData.uMode = obj.uMode;
            obj.schemeData.tMode = 'forward';
            obj.schemeData.hamFunc = params.hamFunc;
            obj.schemeData.partialFunc = params.partialFunc;
            
            % Populate extra arguments.
            obj.HJIextraArgs.quiet = params.quiet;
            obj.HJIextraArgs.visualize.valueSet = 1;
            obj.HJIextraArgs.visualize.initialValueSet = 0;
            obj.HJIextraArgs.visualize.figNum = 1; %set figure number
            obj.HJIextraArgs.visualize.deleteLastPlot = true; 
            obj.HJIextraArgs.visualize.viewGrid = false;
            obj.HJIextraArgs.visualize.viewAxis = [params.lowEnv(1) params.upEnv(1) ...
                                                    params.lowEnv(2) params.upEnv(2) ...
                                                    0 1];
            obj.HJIextraArgs.visualize.xTitle = "$p^x$";
            obj.HJIextraArgs.visualize.yTitle = "$p^y$";
            obj.HJIextraArgs.visualize.zTitle = "$P(\beta = 0)$";
            obj.HJIextraArgs.visualize.fontSize = 15;

            % Uncomment if you want to see a 2D slice
            obj.HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
            obj.HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
            obj.HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D
            
            fprintf('------ Human Predictor Setup -------\n');
            fprintf('   human type: %s\n', params.humanType);
            fprintf('   x0: [%d, %d, %d]\n', obj.xcurr(1), obj.xcurr(2), obj.xcurr(3));
            fprintf('   target radius: %d\n', obj.targetR);
            fprintf('   uMode: %s\n', obj.uMode);
            fprintf('--------------------------------------\n');
            
            % Run first prediction given initial state measurement.
            %fprintf('Running first prediction ..........\n');
            %obj.updatePredictions();
        end
        
        %% Updates the current state of the human. 
        %  This includes updating the distribution over beta = 0!
        function updateState(obj, xt1, ut)
            
            % xcurr = [px_t, py_t, pbeta_t]
            if ~iscell(obj.xcurr)
                xt = num2cell(obj.xcurr);
            end
            pbeta_t1 = obj.human.betaPosterior(xt, ut);
            obj.xcurr = [xt1; pbeta_t1];
        end
        
        %% Predicts the FRS of the human.
        function updatePredictions(obj)
            
            % Target set.
            data = shapeSphere(obj.grid, obj.xcurr, obj.targetR);
            
            % Run HJI computation 
            [obj.preds, obj.timeDisc, ~] = HJIPDE_solve(data, obj.tau, ...
                obj.schemeData, obj.minWith, obj.HJIextraArgs);
        end
        
        %% Gets the current predictions and timestamps
        function [preds, times] = getPredictions(obj)
            preds = obj.preds;
            times = obj.timeDisc;
        end
    end
end

