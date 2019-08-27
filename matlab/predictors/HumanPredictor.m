classdef HumanPredictor < handle
    %HUMANPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        human       % (obj) dynamical system representing human model
        zcurr       % (arr) vector of the current [x,y,P(beta=0)] state
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
            obj.zcurr = params.z0;
            obj.tau = params.tMin:params.dt:params.tMax;
            obj.grid = params.predGrid;
            obj.targetR = params.targetRad;
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
            obj.HJIextraArgs.visualize = 0;
            
            % since we have a finite compute grid, we may not want to 
            % trust values near the boundary of grid
            HJIextraArgs.ignoreBoundary = 0; 

            fprintf('------ Human Predictor Setup -------\n');
            fprintf('   human type: %s\n', params.humanType);
            fprintf('   z0: [%d, %d, %d]\n', obj.zcurr(1), obj.zcurr(2), obj.zcurr(3));
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
            
            % zcurr = [px_t, py_t, pbeta_t]
            if ~iscell(obj.zcurr)
                zt = num2cell(obj.zcurr);
            end
            pbeta_t1 = obj.human.betaPosterior(zt, ut);
            obj.zcurr = [xt1; pbeta_t1];
        end
        
        %% Predicts the FRS of the human.
        function updatePredictions(obj)
            
            % Target set.
            data = shapeSphere(obj.grid, obj.zcurr, obj.targetR);
            
            %Run HJI computation 
            [obj.preds, obj.timeDisc, ~] = HJIPDE_solve_pred(data, obj.tau, ...
               obj.schemeData, obj.minWith, obj.HJIextraArgs);
     
        end
        
        %% Gets the current predictions and timestamps
        function [preds, times] = getPredictions(obj)
            preds = obj.preds;
            times = obj.timeDisc;
        end
    end
end

